"""Campaign configuration and run state.

A campaign is a question ("what does this protein family do that nobody has
described?"); a run is one attempt at it, with its own directory, its own
rubric version and its own manifest. Runs are kept rather than overwritten
because the taste stage learns from the history: the point of keeping a
campaign's set-aside candidates is that they are half the training signal.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from .bio.context import KnownSystem
from .rubric import Rubric, default_rubric
from .schemas import Manifest, StageRecord

_ASSETS = Path(__file__).resolve().parent / "assets"


@dataclass
class CampaignConfig:
    """Everything that defines a campaign, loaded from ``campaign.yaml``."""

    name: str
    family: str = "RT"
    question: str = ""
    source: str = "local:corpus"
    catalogue: str = "builtin"
    backend: str = "rubric"
    max_workers: int = 8

    # Search thresholds. These are campaign policy, not hard-coded constants,
    # because the right cut-off depends on the family's divergence.
    profile_min_z: float = 8.0
    profile_max_evalue: float = 1e-5
    profile_min_core_coverage: float = 0.7
    # Labelling significance. E-value plus coverage, not a raw bit score:
    # bit score scales with alignment length, so a fixed cut-off silently
    # turns every short protein family into an "unknown".
    # Labelling significance: a z-score against each family profile's own
    # shuffled-decoy null, which is comparable across families of different
    # lengths. A raw bit score is not, and a fixed E-value cut-off only
    # partly is.
    label_min_z: float = 6.0
    label_min_coverage: float = 0.20
    sensitive_label_min_z: float = 4.0
    sensitive_label_min_coverage: float = 0.15
    # Anchor recall for the reference path in the methods check.
    sensitive_label_max_evalue: float = 1e-3
    # -log10(E) above which a partner counts as a known family in triage.
    partner_known_neglog_e_cut: float = 3.0
    neighbourhood_bp: int = 12_000
    operon_max_gap: int = 60
    min_novelty: float = 0.35
    max_candidates: int = 60
    background_loci: int = 300
    cluster_identity: float = 0.9

    extras: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def load(cls, path: Path) -> "CampaignConfig":
        doc = yaml.safe_load(Path(path).read_text()) or {}
        known = {f for f in cls.__dataclass_fields__ if f != "extras"}
        kwargs = {k: v for k, v in doc.items() if k in known}
        extras = {k: v for k, v in doc.items() if k not in known}
        return cls(**kwargs, extras=extras)

    def resolved_source(self, base: Path) -> str:
        """Resolve a relative ``local:`` path against the campaign directory."""
        if self.source.startswith("local:"):
            p = Path(self.source[6:])
            if not p.is_absolute():
                p = (base / p).resolve()
            return f"local:{p}"
        return self.source


def load_catalogue(spec: str = "builtin", *, base: Path | None = None) -> tuple[list[KnownSystem], list[dict], dict]:
    """Load the described-systems catalogue.

    Returns ``(systems, reproduction_targets, raw_document)``.
    """
    path = _ASSETS / "known_systems.yaml"
    if spec != "builtin":
        p = Path(spec)
        path = p if p.is_absolute() or base is None else (base / p)
    doc = yaml.safe_load(Path(path).read_text())
    systems = [
        KnownSystem(
            system_id=s["system_id"],
            name=s["name"],
            required_labels=list(s.get("required_labels", [])),
            forbidden_labels=list(s.get("forbidden_labels", [])),
            optional_labels=list(s.get("optional_labels", [])),
            requires_array=bool(s.get("requires_array", False)),
            forbids_array=bool(s.get("forbids_array", False)),
            array_repeat_range=(
                tuple(s["array_repeat_range"]) if s.get("array_repeat_range") else None
            ),
            typical_anchor_aa=(
                tuple(s["typical_anchor_aa"]) if s.get("typical_anchor_aa") else None
            ),
            reference=s.get("reference", ""),
            notes=s.get("notes", ""),
        )
        for s in doc.get("systems", [])
    ]
    return systems, list(doc.get("reproduction_targets", [])), doc


class Run:
    """One execution of a campaign, on disk."""

    def __init__(self, directory: Path, config: CampaignConfig, *, run_id: str | None = None):
        self.config = config
        self.run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        self.dir = Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)
        for sub in ("reports", "bench", "transcripts"):
            (self.dir / sub).mkdir(exist_ok=True)
        self.manifest = self._load_manifest()

    # -- manifest -----------------------------------------------------------

    @property
    def manifest_path(self) -> Path:
        return self.dir / "manifest.json"

    def _load_manifest(self) -> Manifest:
        if self.manifest_path.exists():
            return Manifest.model_validate_json(self.manifest_path.read_text())
        return Manifest(campaign=self.config.name, run_id=self.run_id)

    def save_manifest(self) -> None:
        self.manifest_path.write_text(self.manifest.model_dump_json(indent=1))

    def record_stage(self, record: StageRecord) -> None:
        self.manifest.stages = [s for s in self.manifest.stages if s.stage != record.stage]
        self.manifest.stages.append(record)
        self.save_manifest()

    # -- artifacts ----------------------------------------------------------

    def write_json(self, name: str, payload: Any) -> Path:
        path = self.dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if hasattr(payload, "model_dump_json"):
            path.write_text(payload.model_dump_json(indent=1))
        else:
            path.write_text(json.dumps(payload, indent=1, default=_encode))
        return path

    def read_json(self, name: str) -> Any:
        path = self.dir / name
        if not path.exists():
            return None
        return json.loads(path.read_text())

    def write_text(self, name: str, text: str) -> Path:
        path = self.dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def exists(self, name: str) -> bool:
        return (self.dir / name).exists()

    # -- rubric -------------------------------------------------------------

    @property
    def rubric_path(self) -> Path:
        return self.dir / "rubric.yaml"

    def rubric(self) -> Rubric:
        """The rubric this run uses, pinned to the run directory on first use.

        Pinning matters: the taste stage writes a *new* version rather than
        editing the one a completed run was judged under, so a report always
        remains readable against the standard that produced it.
        """
        if self.rubric_path.exists():
            return Rubric.load(self.rubric_path)
        campaign_rubric = self.dir.parent.parent / "rubric.yaml"
        r = Rubric.load(campaign_rubric) if campaign_rubric.exists() else default_rubric()
        r.save(self.rubric_path)
        return r


def _encode(obj: Any) -> Any:
    if hasattr(obj, "model_dump"):
        return obj.model_dump()
    if hasattr(obj, "__dict__"):
        return {k: v for k, v in vars(obj).items() if not k.startswith("_")}
    return str(obj)


def latest_run(campaign_dir: Path) -> Path | None:
    """Most recent run directory for a campaign, if any."""
    runs = sorted((Path(campaign_dir) / "runs").glob("*/"), key=lambda p: p.name)
    return runs[-1] if runs else None
