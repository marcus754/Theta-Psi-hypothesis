# -*- coding: utf-8 -*-
"""
Data loaders for background-level fits.

For this repository we keep a minimal dimensionless demo dataset for CC (H(z))
to avoid unit conversion complexity. Real applications should use
`background.components.omegas_from_km_s_Mpc` to convert H0 and transform units.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple
import csv
import json
from datetime import date


@dataclass
class CCData:
    z: List[float]
    H: List[float]          # dimensionless H(z) in code units (8πG=1, c=1)
    sigma: List[float]      # uncertainties in same units


def load_cc_dimless(path: str | Path) -> CCData:
    path = Path(path)
    z: List[float] = []
    H: List[float] = []
    s: List[float] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            z.append(float(row["z"]))
            H.append(float(row["H"]))
            s.append(float(row["sigma_H"]))
    return CCData(z=z, H=H, sigma=s)


@dataclass
class SNData:
    z: List[float]
    mu: List[float]         # dimensionless proxy for distance modulus (offset absorbed elsewhere)
    sigma: List[float]


def load_sn_dimless(path: str | Path) -> SNData:
    path = Path(path)
    z: List[float] = []
    mu: List[float] = []
    s: List[float] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            z.append(float(row["z"]))
            mu.append(float(row["mu"]))
            s.append(float(row["sigma_mu"]))
    return SNData(z=z, mu=mu, sigma=s)


@dataclass
class BAOData:
    z: List[float]
    DV: List[float]         # dimensionless DV(z)
    sigma: List[float]


def load_bao_dimless(path: str | Path) -> BAOData:
    path = Path(path)
    z: List[float] = []
    dv: List[float] = []
    s: List[float] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            z.append(float(row["z"]))
            dv.append(float(row["DV"]))
            s.append(float(row["sigma_DV"]))
    return BAOData(z=z, DV=dv, sigma=s)


def load_compact_targets_json(path: str | Path) -> dict:
    """
    Load compact-object target summary from JSON.

    Supported schemas:
    - gaussian (legacy/default):
      z_surface, sigma_z_surface,
      delay_proxy, sigma_delay_proxy,
      n_peak, sigma_n_peak
    - constraints (one-sided soft limits):
      mode="constraints",
      optional:
        z_surface_min, sigma_z_surface_min,
        delay_proxy_min, sigma_delay_proxy_min,
        n_peak_max, sigma_n_peak_max
    """
    p = Path(path)
    data = json.loads(p.read_text(encoding="utf-8"))
    mode = str(data.get("mode", "gaussian")).strip().lower()
    if mode == "gaussian":
        required = (
            "z_surface",
            "sigma_z_surface",
            "delay_proxy",
            "sigma_delay_proxy",
            "n_peak",
            "sigma_n_peak",
        )
        miss = [k for k in required if k not in data]
        if miss:
            raise ValueError(f"Missing compact target keys: {', '.join(miss)}")
        out = {"mode": "gaussian"}
        out.update({k: float(data[k]) for k in required})
        return out
    if mode == "constraints":
        out = {"mode": "constraints"}
        keys = (
            "z_surface_min",
            "sigma_z_surface_min",
            "delay_proxy_min",
            "sigma_delay_proxy_min",
            "n_peak_max",
            "sigma_n_peak_max",
        )
        for k in keys:
            if k in data and data[k] is not None:
                out[k] = float(data[k])
        if not any(k in out for k in ("z_surface_min", "delay_proxy_min", "n_peak_max")):
            raise ValueError("compact constraints mode needs at least one bound (z_surface_min, delay_proxy_min, n_peak_max)")
        return out
    raise ValueError(f"Unknown compact target mode: {mode}")


def load_compact_star_targets_json(path: str | Path) -> dict:
    """
    Load stationary compact-star targets from JSON.

    Expected keys:
      z_surface, sigma_z_surface, delay_proxy, sigma_delay_proxy
    """
    p = Path(path)
    data = json.loads(p.read_text(encoding="utf-8"))
    required = (
        "z_surface",
        "sigma_z_surface",
        "delay_proxy",
        "sigma_delay_proxy",
    )
    miss = [k for k in required if k not in data]
    if miss:
        raise ValueError(f"Missing compact-star target keys: {', '.join(miss)}")
    return {k: float(data[k]) for k in required}


def load_forward_targets_json(path: str | Path) -> dict:
    """
    Load forward-observable targets (EHT proxies; optional legacy echo proxy).

    Expected keys:
      ring_m87_uas, sigma_ring_m87_uas,
      ring_sgra_uas, sigma_ring_sgra_uas
    """
    p = Path(path)
    data = json.loads(p.read_text(encoding="utf-8"))
    required = (
        "ring_m87_uas",
        "sigma_ring_m87_uas",
        "ring_sgra_uas",
        "sigma_ring_sgra_uas",
    )
    miss = [k for k in required if k not in data]
    if miss:
        raise ValueError(f"Missing forward target keys: {', '.join(miss)}")
    out = {k: float(data[k]) for k in required}
    if "echo_delay_ms" in data and "sigma_echo_delay_ms" in data:
        out["echo_delay_ms"] = float(data["echo_delay_ms"])
        out["sigma_echo_delay_ms"] = float(data["sigma_echo_delay_ms"])
    return out


def load_strong_field_bounds_json(path: str | Path) -> dict:
    """
    Load strong-field bounds registry.

    Expected top-level shape:
      {
        "entries": [
          {
            "name": "...",
            "type": "lower" | "upper",
            "observable": "z_surface" | "delay_proxy" | "n_peak",
            "value": <float>,
            "sigma": <float>,
            ...
          },
          ...
        ]
      }
    """
    p = Path(path)
    data = json.loads(p.read_text(encoding="utf-8"))
    if "entries" not in data or not isinstance(data["entries"], list):
        raise ValueError("Strong-field bounds JSON must contain list key 'entries'")
    out_entries = []
    for i, e in enumerate(data["entries"]):
        if not isinstance(e, dict):
            raise ValueError(f"Entry #{i} must be an object")
        miss = [k for k in ("name", "type", "observable", "value", "sigma") if k not in e]
        if miss:
            raise ValueError(f"Entry #{i} missing keys: {', '.join(miss)}")
        typ = str(e["type"]).strip().lower()
        if typ not in ("lower", "upper"):
            raise ValueError(f"Entry #{i} has invalid type: {typ}")
        obs = str(e["observable"]).strip()
        if obs not in ("z_surface", "delay_proxy", "n_peak"):
            raise ValueError(f"Entry #{i} has invalid observable: {obs}")
        observed_on = str(e.get("observed_on", "")).strip()
        if observed_on:
            try:
                date.fromisoformat(observed_on)
            except Exception as ex:
                raise ValueError(f"Entry #{i} has invalid observed_on date: {observed_on}") from ex
        out_entries.append(
            {
                "name": str(e["name"]),
                "type": typ,
                "observable": obs,
                "value": float(e["value"]),
                "sigma": float(e["sigma"]),
                "source": str(e.get("source", "")),
                "source_url": str(e.get("source_url", "")),
                "observed_on": observed_on,
                "note": str(e.get("note", "")),
            }
        )
    return {"entries": out_entries}


def load_forward_event_catalog_json(path: str | Path) -> list[dict]:
    """
    Load event-level forward catalog.

    Expected top-level:
      {"events": [ ... ]}

    Event schemas:
    - ring_uas:
      {name, observable="ring_uas", target, sigma, mass_msun, distance_m}
    - echo_delay_ms:
      {name, observable="echo_delay_ms", target, sigma, echo_mass_msun}
    """
    p = Path(path)
    data = json.loads(p.read_text(encoding="utf-8"))
    if "events" not in data or not isinstance(data["events"], list):
        raise ValueError("forward event catalog must contain list key 'events'")
    out = []
    for i, e in enumerate(data["events"]):
        if not isinstance(e, dict):
            raise ValueError(f"Event #{i} must be object")
        miss = [k for k in ("name", "observable", "target", "sigma") if k not in e]
        if miss:
            raise ValueError(f"Event #{i} missing keys: {', '.join(miss)}")
        obs = str(e["observable"]).strip()
        if obs == "ring_uas":
            for k in ("mass_msun", "distance_m"):
                if k not in e:
                    raise ValueError(f"Event #{i} ring_uas missing key: {k}")
            out.append(
                {
                    "name": str(e["name"]),
                    "observable": obs,
                    "target": float(e["target"]),
                    "sigma": float(e["sigma"]),
                    "mass_msun": float(e["mass_msun"]),
                    "distance_m": float(e["distance_m"]),
                    "source": str(e.get("source", "")),
                }
            )
        elif obs == "echo_delay_ms":
            if "echo_mass_msun" not in e:
                raise ValueError(f"Event #{i} echo_delay_ms missing key: echo_mass_msun")
            out.append(
                {
                    "name": str(e["name"]),
                    "observable": obs,
                    "target": float(e["target"]),
                    "sigma": float(e["sigma"]),
                    "echo_mass_msun": float(e["echo_mass_msun"]),
                    "source": str(e.get("source", "")),
                }
            )
        else:
            raise ValueError(f"Event #{i} invalid observable: {obs}")
    return out
