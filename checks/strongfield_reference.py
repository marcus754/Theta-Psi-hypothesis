# -*- coding: utf-8 -*-
"""
Reference-level strong-field check for the Θ–Ψ compact branch.

This module treats the public forward event report as the finished
event-level comparator for the compact branch:
- ring observables for M87* and Sgr A*;
- GW echo-scale as a downstream strong-field observable.

The conservative compact-bound registry remains a separate diagnostic layer.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class StrongFieldEventReference:
    name: str
    predicted: float
    target: float
    sigma: float
    pull_sigma: float
    within_tolerance: bool


@dataclass(frozen=True)
class StrongFieldReferenceSummary:
    n_events: int
    max_abs_pull_sigma: float
    within_tolerance: bool
    events: tuple[StrongFieldEventReference, ...]


def strongfield_public_reference_check(
    *,
    report_json_path: str | Path = "results/public_forward/forward_event_report.json",
    pull_tolerance_sigma: float = 1.0,
) -> StrongFieldReferenceSummary:
    """
    Load the public strong-field event report and summarize whether the
    event-level comparator stays within the chosen tolerance.
    """
    path = Path(report_json_path)
    data = json.loads(path.read_text(encoding="utf-8"))
    items = list(data.get("events", []))
    events: list[StrongFieldEventReference] = []
    max_abs_pull = 0.0
    for e in items:
        pull = float(e["pull_sigma"])
        max_abs_pull = max(max_abs_pull, abs(pull))
        events.append(
            StrongFieldEventReference(
                name=str(e["name"]),
                predicted=float(e["predicted"]),
                target=float(e["target"]),
                sigma=float(e["sigma"]),
                pull_sigma=pull,
                within_tolerance=abs(pull) <= float(pull_tolerance_sigma),
            )
        )
    return StrongFieldReferenceSummary(
        n_events=len(events),
        max_abs_pull_sigma=max_abs_pull,
        within_tolerance=max_abs_pull <= float(pull_tolerance_sigma),
        events=tuple(events),
    )
