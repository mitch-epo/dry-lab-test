"""``studio`` -- the command line for the molecular design studio."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .campaign import CampaignConfig, latest_run
from .pipeline import STAGES, run_campaign
from .rubric import Rubric, default_rubric

app = typer.Typer(
    add_completion=False,
    help="A dry-lab studio: from an uncharacterised protein family to a short "
         "list of lab-testable hypotheses.",
)
console = Console()


def _resolve_run(campaign: Path, run_id: str | None) -> Path:
    if run_id:
        return campaign / "runs" / run_id
    found = latest_run(campaign)
    if found is None:
        raise typer.BadParameter(f"no runs under {campaign / 'runs'}")
    return found


@app.command("build-corpus")
def build_corpus_cmd(
    out: Path = typer.Argument(..., help="directory to write the corpus into"),
    seed: int = typer.Option(20240917, help="generator seed"),
    scale: float = typer.Option(1.0, help="scale every planted instance count"),
    background: int = typer.Option(120, help="number of background contigs"),
):
    """Generate the offline corpus.

    Every sequence is simulated. The corpus exists so the pipeline can run and
    be scored where the sequence databases are unreachable.
    """
    from .corpus.build import build_corpus

    with console.status("generating corpus..."):
        corpus = build_corpus(seed=seed, scale=scale, n_background_contigs=background)
        corpus.write(out)
    console.print(Panel.fit(
        f"[bold]{corpus.meta['n_contigs']}[/] contigs  "
        f"[bold]{corpus.meta['n_genes']}[/] genes  "
        f"[bold]{corpus.meta['n_arrays']}[/] arrays  "
        f"[bold]{corpus.meta['n_planted']}[/] planted loci\n"
        f"[yellow]{corpus.meta['warning']}[/]\n"
        f"ground truth withheld in TRUTH.withheld.json (read only by `studio score`)",
        title=f"corpus written to {out}",
    ))


@app.command()
def run(
    campaign: Path = typer.Argument(..., help="campaign directory containing campaign.yaml"),
    backend: str = typer.Option(None, help="rubric | anthropic[:model] | replay"),
    run_id: str = typer.Option(None, help="reuse or name a run directory"),
    stages: str = typer.Option(None, help=f"comma-separated subset of: {','.join(STAGES)}"),
    force: bool = typer.Option(False, help="continue past a failed reproduce gate"),
):
    """Run a campaign end to end."""
    selected = [s.strip() for s in stages.split(",")] if stages else None
    if selected:
        bad = [s for s in selected if s not in STAGES]
        if bad:
            raise typer.BadParameter(f"unknown stage(s): {', '.join(bad)}")

    outcome = run_campaign(
        campaign, backend=backend, run_id=run_id, stages=selected,
        force=force, log=lambda m: console.print(m),
    )

    console.print()
    if outcome.stopped_at:
        console.print(Panel.fit(
            "[red]The methods check did not pass.[/] No novel candidate was "
            f"proposed.\nSee [bold]{outcome.run.dir / 'reproduce.txt'}[/]",
            title="run stopped",
        ))
        raise typer.Exit(code=2)

    promoted = outcome.promoted
    packages = outcome.packages
    table = Table(
        title=(
            f"{len(packages)} candidate system(s) from {len(promoted)} promoted loci"
            if packages else "nothing promoted"
        )
    )
    table.add_column("lead candidate")
    table.add_column("loci", justify="right")
    table.add_column("architecture")
    table.add_column("proposal")
    for pkg in packages:
        table.add_row(
            pkg.candidate_id, str(pkg.n_instances), pkg.architecture, pkg.title
        )
    console.print(table)
    if not promoted:
        console.print(
            "[yellow]A survey may end with a single candidate worth testing, or "
            "with none. The eliminations are in verdicts.json.[/]"
        )
    console.print(f"\nrun directory: [bold]{outcome.run.dir}[/]")
    console.print(f"summary: [bold]{outcome.run.dir / 'summary.md'}[/]")


@app.command()
def score(
    campaign: Path = typer.Argument(..., help="campaign directory"),
    corpus: Path = typer.Option(None, help="corpus directory (default: the campaign's source)"),
    run_id: str = typer.Option(None, help="which run to score (default: the latest)"),
):
    """Score a run against the corpus's withheld ground truth.

    Only meaningful for the local corpus; a real campaign has no answer key.
    """
    from .scoring import score_run

    run_dir = _resolve_run(campaign, run_id)
    if corpus is None:
        config = CampaignConfig.load(campaign / "campaign.yaml")
        spec = config.resolved_source(campaign)
        if not spec.startswith("local:"):
            raise typer.BadParameter("scoring needs a local corpus")
        corpus = Path(spec[6:])

    report = score_run(run_dir, corpus)
    console.print(report.text())
    raise typer.Exit(code=0 if report.clean else 1)


@app.command()
def show(
    campaign: Path = typer.Argument(...),
    run_id: str = typer.Option(None),
    candidate: str = typer.Option(None, help="print one candidate's report"),
):
    """Print a run's summary, or one candidate report."""
    run_dir = _resolve_run(campaign, run_id)
    if candidate:
        path = run_dir / "reports" / f"{candidate}.md"
        if not path.exists():
            raise typer.BadParameter(f"no report at {path}")
        console.print(path.read_text())
        return
    summary = run_dir / "summary.md"
    console.print(summary.read_text() if summary.exists() else f"no summary in {run_dir}")


@app.command("rubric")
def rubric_cmd(
    campaign: Path = typer.Argument(None, help="campaign directory, or omit for the default"),
    run_id: str = typer.Option(None),
    prompt: bool = typer.Option(False, "--prompt", help="render as the reviewer instruction block"),
    init: bool = typer.Option(False, "--init", help="write the default rubric into the campaign"),
):
    """Show or initialise the rubric."""
    if init:
        if campaign is None:
            raise typer.BadParameter("--init needs a campaign directory")
        path = campaign / "rubric.yaml"
        default_rubric().save(path)
        console.print(f"wrote {path}")
        return

    if campaign is None:
        r = default_rubric()
    else:
        run_dir = _resolve_run(campaign, run_id) if (campaign / "runs").exists() else None
        candidate_paths = [
            run_dir / "rubric.yaml" if run_dir else None,
            campaign / "rubric.yaml",
        ]
        path = next((p for p in candidate_paths if p and p.exists()), None)
        r = Rubric.load(path) if path else default_rubric()

    console.print(r.as_prompt() if prompt else
                  f"rubric v{r.version} ({r.name}): "
                  f"{len(r.gates)} gates, {len(r.scores)} scores, "
                  f"{len(r.guidance)} guidance lines")
    if not prompt:
        for g in r.gates:
            console.print(f"  gate  {g.gate_id:26s} {g.expr}")
        for s in r.scores:
            console.print(f"  score {s.score_id:26s} w={s.weight:<4g} {s.expr}")


@app.command()
def doctor():
    """Report what this environment can reach and which backends are available."""
    table = Table(title="studio doctor")
    table.add_column("check")
    table.add_column("status")
    table.add_column("detail")

    import os

    table.add_row(
        "anthropic backend",
        "[green]ready[/]" if os.environ.get("ANTHROPIC_API_KEY") else "[yellow]unavailable[/]",
        "ANTHROPIC_API_KEY is set" if os.environ.get("ANTHROPIC_API_KEY")
        else "set ANTHROPIC_API_KEY and install the 'live' extra",
    )
    table.add_row("rubric backend", "[green]ready[/]", "deterministic, always available")

    try:
        from .sources.live import probe
        for p in probe(timeout=8.0):
            table.add_row(
                f"live source: {p.endpoint}",
                "[green]reachable[/]" if p.ok else "[red]blocked[/]",
                p.detail or f"HTTP {p.status}" if p.status else p.detail,
            )
    except Exception as exc:
        table.add_row("live source", "[red]error[/]", f"{type(exc).__name__}: {exc}")

    console.print(table)


@app.command("list-runs")
def list_runs(campaign: Path = typer.Argument(...)):
    """List a campaign's runs."""
    import json

    runs_dir = campaign / "runs"
    if not runs_dir.exists():
        console.print(f"no runs under {runs_dir}")
        return
    table = Table(title=f"runs of {campaign.name}")
    for col in ("run", "backend", "rubric", "gate", "candidates", "promoted"):
        table.add_column(col)
    for d in sorted(runs_dir.iterdir()):
        mf = d / "manifest.json"
        if not mf.exists():
            continue
        m = json.loads(mf.read_text())
        gate = m.get("gate_passed")
        table.add_row(
            d.name, m.get("backend", ""), f"v{m.get('rubric_version', '?')}",
            "[green]pass[/]" if gate else ("[red]fail[/]" if gate is False else "-"),
            str(m.get("counts", {}).get("candidates", "-")),
            str(m.get("counts", {}).get("promoted", "-")),
        )
    console.print(table)


if __name__ == "__main__":
    app()
