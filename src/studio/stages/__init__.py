"""The seven stages of the workflow.

    survey -> reproduce -> search -> report -> triage -> bench -> taste
              (gate)                                              (feeds back
                                                                   into survey)

``reproduce`` is a gate, not a step: if the pipeline cannot recover what is
already described, the run stops. ``taste`` closes the loop by rewriting the
rubric that ``report`` and ``triage`` are held to.
"""

from .bench import build_package, package_markdown, run_bench
from .report import run_report
from .reproduce import run_reproduce
from .search import run_search
from .survey import FamilyLabeller, build_reference_index, run_survey
from .taste import run_taste
from .triage import run_triage

__all__ = [
    "FamilyLabeller", "build_package", "build_reference_index", "package_markdown",
    "run_bench", "run_report", "run_reproduce", "run_search", "run_survey",
    "run_taste", "run_triage",
]
