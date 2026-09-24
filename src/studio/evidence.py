"""How strong a measurement lets a claim be.

Shared deliberately. The reporter assigns a strength to each piece of evidence
and the reviewer checks those assignments, so if the two used different rules
the reviewer would flag a systematic overclaim on every run -- and that finding
would be an artefact of the disagreement, not a fact about the reports. A
built-in mismatch is worse than no check, because the taste stage would then
learn from it and write guidance about a problem nobody has.

With one source, a flagged overclaim means the reporter genuinely reached past
its evidence.
"""

from __future__ import annotations

from typing import Any

from .schemas import Strength

RANK = {
    Strength.ABSENT: 0,
    Strength.WEAK: 1,
    Strength.SUGGESTIVE: 2,
    Strength.STRONG: 3,
    Strength.DIRECT: 4,
}


def strength_for(key: str, metrics: dict[str, Any]) -> Strength | None:
    """The strongest label a measurement justifies for an evidence key.

    Returns ``None`` for keys with no rule, which the reviewer treats as
    "not checkable" rather than as a pass.
    """
    if key == "catalytic_core":
        return Strength.DIRECT if metrics.get("triad_intact") else Strength.ABSENT

    if key == "family_membership":
        z = float(metrics.get("profile_z", 0.0))
        cov = float(metrics.get("core_coverage", 0.0))
        if z >= 20 and cov >= 0.8:
            return Strength.DIRECT
        return Strength.STRONG if z >= 10 else Strength.SUGGESTIVE

    if key == "partner_association":
        p = float(metrics.get("association_neglog_p", 0.0))
        n = int(metrics.get("n_independent_loci", 0))
        conserved = float(metrics.get("partner_conservation", 0.0))
        if p >= 6 and n >= 8 and conserved >= 0.8:
            return Strength.STRONG
        if p >= 3 and n >= 3:
            return Strength.SUGGESTIVE
        return Strength.WEAK

    if key == "noncoding_component":
        q = float(metrics.get("array_quality", 0.0))
        gap = float(metrics.get("array_operon_gap_bp", 10 ** 9))
        if q <= 0.0:
            return Strength.ABSENT
        if gap > 3000:
            # Present, but too far away to be part of the same system.
            return Strength.WEAK
        if q >= 0.7:
            return Strength.STRONG
        return Strength.SUGGESTIVE if q >= 0.4 else Strength.WEAK

    if key == "taxonomic_spread":
        phyla = int(metrics.get("n_phyla", 0))
        if metrics.get("clade_restricted", True):
            return Strength.WEAK
        return Strength.STRONG if phyla >= 3 else Strength.SUGGESTIVE

    if key == "architecture":
        novelty = float(metrics.get("novelty", 0.0))
        return Strength.STRONG if novelty >= 0.6 else Strength.SUGGESTIVE

    if key in ("array_placement", "partner_is_known"):
        return Strength.STRONG

    return None


def overclaims(claims: list[dict[str, Any]], metrics: dict[str, Any]) -> list[str]:
    """Claims asserted more strongly than their measurement supports."""
    out: list[str] = []
    for claim in claims:
        key = claim.get("evidence_key", "")
        supported = strength_for(key, metrics)
        if supported is None:
            continue
        try:
            asserted = Strength(claim.get("strength", "suggestive"))
        except ValueError:
            continue
        if RANK[asserted] > RANK[supported]:
            out.append(
                f"{key}: asserted as {asserted.value} but the measurement "
                f"supports only {supported.value}"
            )
    return out
