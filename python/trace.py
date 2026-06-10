"""
====================================================
TRACE — Debug and Trace Utilities
====================================================

Purpose:
    Provides trace output capabilities for debugging and inspecting
    every stage of the pipeline.

    When trace=True is enabled, every pipeline stage writes its
    inputs, outputs, metrics, and execution time to JSON files
    in the traces/ directory.

    Trace files produced:
        01_raw_centroids.json        — Raw centroids from centroiding
        02_filtered_centroids.json    — After magnitude/brightest filtering
        03_star_ids.json              — Star identification results
        04_attitude.json              — Final attitude solution
        05_pipeline_output.json       — Complete pipeline output summary

Reference:
    Inspired by the C++ LOST diagnostic outputs via --print-*
    and --compare-* CLI options.
"""

from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional


@dataclass
class TraceEntry:
    """A single trace entry recording a pipeline stage execution."""

    stage_name: str
    stage_number: int
    algorithm: str
    input_schema: Dict[str, str]
    output_schema: Dict[str, str]
    execution_time_s: float
    input_metrics: Dict[str, Any] = field(default_factory=dict)
    output_metrics: Dict[str, Any] = field(default_factory=dict)
    config: Dict[str, Any] = field(default_factory=dict)


class TraceSession:
    """Manages trace output for a single pipeline run.

    Collects entries from each stage and writes them to files
    and/or a summary JSON.
    """

    def __init__(self, trace_dir: str = "traces", enabled: bool = True):
        self.enabled = enabled
        self.trace_dir = trace_dir
        self.entries: List[TraceEntry] = []

    def record(self, entry: TraceEntry) -> None:
        if not self.enabled:
            return
        self.entries.append(entry)
        os.makedirs(self.trace_dir, exist_ok=True)
        filename = f"{entry.stage_number:02d}_{entry.stage_name}.json"
        filepath = os.path.join(self.trace_dir, filename)
        with open(filepath, 'w') as f:
            json.dump(asdict(entry), f, indent=2, default=str)

    def write_summary(self) -> None:
        if not self.enabled:
            return
        filepath = os.path.join(self.trace_dir, "00_trace_summary.json")
        with open(filepath, 'w') as f:
            json.dump({
                'num_stages': len(self.entries),
                'total_time_s': sum(e.execution_time_s for e in self.entries),
                'stages': [asdict(e) for e in self.entries],
            }, f, indent=2, default=str)
