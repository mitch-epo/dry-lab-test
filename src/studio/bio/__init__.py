"""Offline bioinformatics kernel for the molecular design studio.

Everything in this package is deterministic and network-free: given the same
sequences it produces the same numbers whether the records came from a live
UniProt query or a local snapshot. That property is what makes the
``reproduce`` stage meaningful -- a methods check is worthless if the methods
drift between runs.
"""

from .align import Alignment, align, bit_score, evalue, profile_align, quick_score
from .cluster import Cluster, TaxonomicSpread, greedy_cluster, non_redundant, taxonomic_spread
from .context import (
    Gene,
    KnownSystem,
    Locus,
    NoveltyAssessment,
    SystemMatch,
    architecture,
    assess_novelty,
    conservation_across_loci,
    core_architecture,
    extract_locus,
    match_known_system,
    neighbours,
    predicted_operon,
)
from .motifs import MotifReport, MotifSet, MotifSpec, RT_MOTIFS, catalytic_triad_intact
from .profile import Profile, build_profile, profile_from_unaligned
from .repeats import RepeatArray, array_transcript_units, find_repeat_arrays
from .search import Hit, SequenceIndex, reciprocal_best_hits
from .seq import (
    Feature,
    clean_protein,
    find_orfs,
    gc_content,
    low_complexity_fraction,
    parse_fasta,
    revcomp,
    shannon_entropy,
    translate,
    write_fasta,
)
from .stats import FisherResult, benjamini_hochberg, fisher_exact, neighbour_enrichment

__all__ = [
    "Alignment", "align", "bit_score", "evalue", "profile_align", "quick_score",
    "Cluster", "TaxonomicSpread", "greedy_cluster", "non_redundant", "taxonomic_spread",
    "Gene", "KnownSystem", "Locus", "NoveltyAssessment", "SystemMatch",
    "architecture", "assess_novelty", "conservation_across_loci", "core_architecture",
    "extract_locus", "match_known_system", "neighbours", "predicted_operon",
    "MotifReport", "MotifSet", "MotifSpec", "RT_MOTIFS", "catalytic_triad_intact",
    "Profile", "build_profile", "profile_from_unaligned",
    "RepeatArray", "array_transcript_units", "find_repeat_arrays",
    "Hit", "SequenceIndex", "reciprocal_best_hits",
    "Feature", "clean_protein", "find_orfs", "gc_content", "low_complexity_fraction",
    "parse_fasta", "revcomp", "shannon_entropy", "translate", "write_fasta",
    "FisherResult", "benjamini_hochberg", "fisher_exact", "neighbour_enrichment",
]
