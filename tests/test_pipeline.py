"""Tests for the rubric, the harness, the stages and an end-to-end campaign.

The end-to-end test runs a real campaign on a small corpus and asserts the two
properties the whole design exists to produce: the planted novel system is
recovered, and every decoy class is eliminated. Those run against the withheld
ground truth, which no stage can see.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from studio.agents import Harness, RubricBackend, Session, TranscriptStore
from studio.agents.backends import ReplayBackend
from studio.campaign import CampaignConfig, Run, load_catalogue
from studio.corpus.build import build_corpus
from studio.pipeline import run_campaign
from studio.rubric import Rubric, RubricExpressionError, default_rubric, evaluate
from studio.scoring import score_run
from studio.stages.taste import auroc, load_labels


# --- the rubric expression evaluator -----------------------------------------


def test_rubric_expressions_evaluate():
    v = {"a": 3.0, "flag": True, "n": 10}
    assert evaluate("a > 2 and flag", v) is True
    assert evaluate("clamp(a / 10.0)", v) == pytest.approx(0.3)
    assert evaluate("1.0 if flag else 0.0", v) == 1.0
    assert evaluate("min(n, 5)", v) == 5
    assert evaluate("not flag", v) is False


@pytest.mark.parametrize(
    "expr",
    [
        "__import__('os').system('true')",
        "open('/etc/passwd').read()",
        "(lambda: 1)()",
        "a.__class__",
        "[x for x in range(10)]",
        "eval('1+1')",
    ],
)
def test_rubric_expressions_refuse_arbitrary_code(expr):
    """The taste stage rewrites the rubric automatically, so the evaluator is
    an attack surface on its own outputs. It must not be a path to execution."""
    with pytest.raises(RubricExpressionError):
        evaluate(expr, {"a": 1})


def test_unknown_names_raise_rather_than_defaulting():
    """A missing metric is a pipeline fault, not a failed candidate."""
    with pytest.raises(RubricExpressionError):
        evaluate("missing_metric > 1", {"a": 1})


def test_default_rubric_gates_on_a_complete_candidate():
    rubric = default_rubric()
    good = {
        "triad_intact": True, "core_coverage": 0.9, "aligned_fraction": 0.6,
        "low_complexity": 0.05, "evalue": 1e-30, "profile_z": 90.0,
        "novelty": 0.9, "partner_conservation": 0.95,
        "partner_known_neglog_e": 0.0, "n_independent_loci": 12,
        "n_clusters": 12, "n_phyla": 3, "clade_restricted": False,
        "association_neglog_p": 12.0, "array_quality": 0.85,
        "array_operon_gap_bp": 150, "n_predictions": 4,
    }
    failed, _reasons, scores, priority = rubric.apply(good)
    assert failed == []
    assert priority > rubric.promote_threshold
    assert rubric.disposition(failed, priority) == "promote"
    assert set(scores) == {s.score_id for s in rubric.scores}


@pytest.mark.parametrize(
    "override,expected_gate",
    [
        ({"triad_intact": False}, "catalytic_core"),
        ({"core_coverage": 0.3}, "domain_coverage"),
        ({"low_complexity": 0.8}, "not_low_complexity"),
        ({"profile_z": 2.0}, "significance"),
        ({"novelty": 0.05}, "novelty"),
        ({"partner_known_neglog_e": 20.0}, "partner_is_not_known"),
        ({"n_independent_loci": 1}, "reproducible_association"),
        ({"n_clusters": 1}, "not_redundant"),
        (
            {"partner_conservation": 0.0, "array_quality": 0.0},
            "has_testable_component",
        ),
        (
            {"partner_conservation": 0.0, "array_operon_gap_bp": 9000},
            "has_testable_component",
        ),
    ],
)
def test_each_gate_fires_on_its_own_failure_mode(override, expected_gate):
    rubric = default_rubric()
    good = {
        "triad_intact": True, "core_coverage": 0.9, "aligned_fraction": 0.6,
        "low_complexity": 0.05, "evalue": 1e-30, "profile_z": 90.0,
        "novelty": 0.9, "partner_conservation": 0.95,
        "partner_known_neglog_e": 0.0, "n_independent_loci": 12,
        "n_clusters": 12, "n_phyla": 3, "clade_restricted": False,
        "association_neglog_p": 12.0, "array_quality": 0.85,
        "array_operon_gap_bp": 150, "n_predictions": 4,
    }
    failed, reasons, _scores, priority = rubric.apply({**good, **override})
    assert expected_gate in failed
    assert reasons and all(r for r in reasons)
    assert rubric.disposition(failed, priority) == "eliminate"


def test_rubric_round_trips_through_yaml(tmp_path):
    rubric = default_rubric()
    path = tmp_path / "r.yaml"
    rubric.save(path)
    back = Rubric.load(path)
    assert back.version == rubric.version
    assert [g.gate_id for g in back.gates] == [g.gate_id for g in rubric.gates]
    assert [s.expr for s in back.scores] == [s.expr for s in rubric.scores]


def test_rubric_prompt_names_every_gate():
    """The deterministic and the model-backed reviewer must see one standard."""
    rubric = default_rubric()
    text = rubric.as_prompt()
    for gate in rubric.gates:
        assert gate.gate_id in text
        assert gate.expr in text
    for score in rubric.scores:
        assert score.score_id in text


# --- the harness --------------------------------------------------------------


def test_harness_runs_sessions_in_order_and_records_transcripts(tmp_path):
    rubric = default_rubric()
    harness = Harness(RubricBackend(rubric), transcripts=tmp_path, max_workers=4)
    metrics = {
        "triad_intact": True, "core_coverage": 0.9, "aligned_fraction": 0.6,
        "low_complexity": 0.05, "evalue": 1e-30, "profile_z": 90.0,
        "novelty": 0.9, "partner_conservation": 0.9,
        "partner_known_neglog_e": 0.0, "n_independent_loci": 9, "n_clusters": 9,
        "n_phyla": 3, "clade_restricted": False, "association_neglog_p": 9.0,
        "array_quality": 0.8, "array_operon_gap_bp": 100, "n_predictions": 3,
    }
    sessions = [
        Session(
            session_id=f"s{i}", role="reviewer", instructions=rubric.as_prompt(),
            task="review", bundle={"candidate_id": f"c{i}", "metrics": metrics},
        )
        for i in range(6)
    ]
    results = harness.run(sessions)
    assert [r.session_id for r in results] == [s.session_id for s in sessions]
    assert all(r.ok for r in results)
    assert harness.stats.sessions == 6
    assert len(list(Path(tmp_path).glob("*.json"))) == 6


def test_harness_isolates_a_failing_session(tmp_path):
    """One bad session must not take the batch down."""

    class Flaky:
        name = "flaky"

        def run(self, session):
            if session.session_id == "bad":
                raise RuntimeError("boom")
            from studio.agents.session import SessionResult

            return SessionResult(session.session_id, session.role, {"ok": True}, self.name)

    harness = Harness(Flaky(), transcripts=tmp_path)
    # Distinct bundles: identical work is deduplicated by design, so a batch of
    # three identical sessions would be one call and test nothing.
    results = harness.run([
        Session(session_id=i, role="reviewer", instructions="", task="",
                bundle={"candidate_id": i})
        for i in ("good1", "bad", "good2")
    ])
    assert [r.ok for r in results] == [True, False, True]
    assert "boom" in results[1].error
    assert harness.stats.failures == 1


def test_identical_sessions_are_answered_once(tmp_path):
    calls = {"n": 0}

    class Counting:
        name = "counting"

        def run(self, session):
            from studio.agents.session import SessionResult

            calls["n"] += 1
            return SessionResult(session.session_id, session.role, {}, self.name)

    harness = Harness(Counting(), transcripts=tmp_path)
    a = Session(session_id="a", role="reviewer", instructions="i", task="t",
                bundle={"x": 1})
    b = Session(session_id="b", role="reviewer", instructions="i", task="t",
                bundle={"x": 1})
    assert a.fingerprint == b.fingerprint
    harness.run([a, b])
    assert calls["n"] == 1
    assert harness.stats.cached == 1


def test_replay_refuses_a_transcript_recorded_against_different_inputs(tmp_path):
    store = TranscriptStore(tmp_path)
    original = Session(session_id="s", role="reviewer", instructions="i", task="t",
                       bundle={"x": 1})
    from studio.agents.session import SessionResult

    store.write(original, SessionResult("s", "reviewer", {"disposition": "promote"},
                                        "rubric", fingerprint=original.fingerprint))
    changed = Session(session_id="s", role="reviewer", instructions="i", task="t",
                      bundle={"x": 2})
    with pytest.raises(ValueError, match="different inputs"):
        ReplayBackend(store).run(changed)


# --- the taste stage ----------------------------------------------------------


def test_auroc_is_half_at_chance_and_one_when_separable():
    assert auroc([1, 2, 3], [1, 2, 3]) == pytest.approx(0.5)
    assert auroc([4, 5, 6], [1, 2, 3]) == 1.0
    assert auroc([1, 2, 3], [4, 5, 6]) == 0.0


def test_labels_file_is_read_when_present(tmp_path):
    path = tmp_path / "labels.yaml"
    path.write_text(yaml.safe_dump({"labels": {"cand-a": "promote", "cand-b": "set_aside"}}))
    labels = load_labels(path)
    assert labels == {"cand-a": "promote", "cand-b": "set_aside"}
    assert load_labels(tmp_path / "absent.yaml") == {}


# --- the catalogue ------------------------------------------------------------


def test_builtin_catalogue_declares_reproduction_targets():
    systems, targets, doc = load_catalogue("builtin")
    assert len(systems) >= 5
    assert targets, "a methods check with no positive control cannot fail"
    ids = {s.system_id for s in systems}
    for t in targets:
        assert t["system_id"] in ids
        assert 0 < t["min_recall"] <= 1 and 0 < t["min_precision"] <= 1
    # Every described system must state something that would rule it out.
    for s in systems:
        assert s.required_labels
        assert s.forbidden_labels or s.forbids_array or s.requires_array


# --- end to end ---------------------------------------------------------------


@pytest.fixture(scope="module")
def small_campaign(tmp_path_factory):
    """A small campaign, run once and shared by the end-to-end assertions."""
    root = tmp_path_factory.mktemp("campaign")
    corpus_dir = root / "corpus"
    build_corpus(seed=4242, scale=0.35, n_background_contigs=25).write(corpus_dir)

    (root / "campaign.yaml").write_text(yaml.safe_dump({
        "name": "test-rt",
        "family": "RT",
        "source": "local:corpus",
        "backend": "rubric",
        "max_workers": 4,
        "background_loci": 80,
        "max_candidates": 40,
    }))
    outcome = run_campaign(root, run_id="t1")
    return root, corpus_dir, outcome


def test_end_to_end_gate_passes(small_campaign):
    _root, _corpus, outcome = small_campaign
    assert outcome.gate_passed, "the methods check must pass on the reference corpus"
    assert not outcome.stopped_at


def test_end_to_end_recovers_the_planted_novel_system(small_campaign):
    root, corpus, outcome = small_campaign
    report = score_run(root / "runs" / "t1", corpus)
    assert report.novel_recovered > 0, "no planted novel locus reached the candidate list"
    assert report.novel_promoted > 0, "the planted novel system was eliminated"


def test_end_to_end_kills_every_decoy_class(small_campaign):
    """The property the workflow is defined by: what it throws away."""
    root, corpus, outcome = small_campaign
    report = score_run(root / "runs" / "t1", corpus)
    leaked = [d.system_id for d in report.decoys if not d.caught]
    assert not leaked, f"decoys promoted: {leaked}"
    assert report.described_promoted == 0, "a described system was promoted as novel"


def test_end_to_end_eliminates_most_candidates(small_campaign):
    """'Typically most candidates are eliminated at this stage.'"""
    _root, _corpus, outcome = small_campaign
    eliminated = [v for v in outcome.verdicts if v.disposition == "eliminate"]
    assert len(eliminated) > len(outcome.verdicts) / 2
    assert all(v.kill_reason for v in eliminated), "every elimination must name its reason"


def test_every_report_carries_evidence_and_falsifiable_predictions(small_campaign):
    _root, _corpus, outcome = small_campaign
    assert outcome.reports
    for r in outcome.reports:
        assert r.evidence, f"{r.candidate_id} has no evidence"
        assert all(e.method for e in r.evidence), "every claim needs a method"
        assert r.alternatives, f"{r.candidate_id} lists no alternative explanation"
    promoted_ids = {v.candidate_id for v in outcome.promoted}
    for r in outcome.reports:
        if r.candidate_id in promoted_ids:
            assert r.predictions, "a promoted candidate must be falsifiable"
            assert all(p.would_falsify for p in r.predictions)


def test_manifest_records_provenance_and_marks_synthetic_data(small_campaign):
    root, _corpus, _outcome = small_campaign
    manifest = json.loads((root / "runs" / "t1" / "manifest.json").read_text())
    assert manifest["synthetic_data"] is True
    assert "SIMULATED" in manifest["source"]["banner"]
    assert manifest["gate_passed"] is True
    assert {s["stage"] for s in manifest["stages"]} >= {
        "survey", "reproduce", "search", "report", "triage", "bench", "taste"
    }
    for stage in manifest["stages"]:
        assert stage["duration_s"] >= 0
        assert stage["notes"]


def test_synthetic_banner_reaches_every_rendered_report(small_campaign):
    root, _corpus, _outcome = small_campaign
    for path in (root / "runs" / "t1" / "reports").glob("*.md"):
        assert "SIMULATED DATA" in path.read_text().splitlines()[0]


def test_bench_handoff_is_per_system_and_leads_with_falsification(small_campaign):
    _root, _corpus, outcome = small_campaign
    for pkg in outcome.packages:
        assert pkg.falsification, "a handoff with nothing that would kill it is not a test"
        assert pkg.controls
        assert any("mutant" in c.lower() for c in pkg.controls)
        assert pkg.n_instances == len(pkg.supporting_loci) + 1
    # One package per architecture, not one per locus.
    assert len(outcome.packages) <= len(outcome.promoted)


def test_taste_stage_refuses_to_close_the_loop_on_its_own_verdicts(small_campaign):
    """Without external labels, weight changes to existing criteria are held."""
    root, _corpus, outcome = small_campaign
    taste = json.loads((root / "runs" / "t1" / "taste.json").read_text())
    assert taste["label_source"] == "verdicts"
    circular = [p for p in taste["proposals"] if p["basis"] == "circular"]
    assert all(not p["applied"] for p in circular)
    assert any("circular" in n or "same rubric" in n for n in taste["notes"])


def test_transcripts_are_written_for_every_session(small_campaign):
    root, _corpus, outcome = small_campaign
    transcripts = list((root / "runs" / "t1" / "transcripts").glob("*.json"))
    assert len(transcripts) >= len(outcome.reports)
    one = json.loads(transcripts[0].read_text())
    assert one["fingerprint"] and one["instructions"] and "bundle" in one
