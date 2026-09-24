"""Data sources: where a campaign's sequences and literature come from."""

from __future__ import annotations

from pathlib import Path

from .base import DataSource, ProteinRecord, Publication, SourceProvenance
from .local import LocalSource


def open_source(spec: str, **kwargs) -> DataSource:
    """Resolve a source specification.

    ``local:<path>`` serves a generated corpus; ``live`` serves the public
    databases. The spec is recorded in the run manifest, so a report always
    carries the identity of the data behind it.
    """
    if spec == "live" or spec.startswith("live:"):
        from .live import LiveSource
        rest = spec[5:] if spec.startswith("live:") else ""
        if rest:
            kwargs.setdefault("family_query", rest)
        return LiveSource(**kwargs)
    if spec.startswith("local:"):
        return LocalSource(Path(spec[6:]), **kwargs)
    return LocalSource(Path(spec), **kwargs)


__all__ = [
    "DataSource", "LocalSource", "ProteinRecord", "Publication",
    "SourceProvenance", "open_source",
]
