# -*- coding: utf-8 -*-
"""
Converters from common public cosmology datasets to the repo's simple
dimensionless CSV formats for CC/SN/BAO demos and event catalog JSON for
strong-field forward checks.

Usage examples:

  python -m fitting.convert_public --cc CC_public.csv --H0_km_s_Mpc 70.0 \
      --outdir fitting/data/real_dimless

  python -m fitting.convert_public --sn Pantheon.csv --outdir fitting/data/real_dimless

  python -m fitting.convert_public --bao BAO_DV_public.csv --outdir fitting/data/real_dimless

  python -m fitting.convert_public --forward_events strong_field_events.csv \
      --outdir fitting/data/real_dimless

Notes:
- CC: expects columns 'z','H','sigma_H' in km/s/Mpc; converts to dimensionless
  H by multiplying with reduced Planck time (via H0_from_km_s_Mpc for H0=1 s^-1).
- SN: expects columns 'z','mu','sigma_mu' (Pantheon-like). μ0 offset is left
  implicit and can be marginalized in the fit; here we just reformat.
- BAO: simplest path expects columns 'z','DV','sigma_DV' (dimensionless or in
  consistent length units); the model compares dimensionless DV curves, so if
  DV is given in Mpc, you may rescale by H0/c to get dimensionless.
- Forward events: empirical event-level catalog currently accepts only
  `observable=ring_uas`; echo-like entries are treated as non-empirical
  diagnostics and are not converted here.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import csv
import json

from fitting.core_api import H0_from_km_s_Mpc


def _ensure_outdir(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)


def convert_cc(input_csv: Path, out_csv: Path, *, H0_km_s_Mpc: float):
    # Convert H(z) [km/s/Mpc] to dimensionless by multiplying by reduced Planck time
    # Using H0_from_km_s_Mpc(1 km/s/Mpc) to get 1 s^-1 in code units
    scale = H0_from_km_s_Mpc(1.0)
    with input_csv.open("r", newline="", encoding="utf-8") as fi, out_csv.open("w", newline="", encoding="utf-8") as fo:
        r = csv.DictReader(fi)
        w = csv.DictWriter(fo, fieldnames=["z", "H", "sigma_H"])
        w.writeheader()
        for row in r:
            z = float(row["z"])
            H = float(row["H"]) * scale
            sH = float(row.get("sigma_H", row.get("sigma", 0.0))) * scale
            w.writerow({"z": z, "H": H, "sigma_H": sH})


def convert_sn(input_csv: Path, out_csv: Path):
    with input_csv.open("r", newline="", encoding="utf-8") as fi, out_csv.open("w", newline="", encoding="utf-8") as fo:
        r = csv.DictReader(fi)
        w = csv.DictWriter(fo, fieldnames=["z", "mu", "sigma_mu"])
        w.writeheader()
        for row in r:
            z = float(row["z"]) if "z" in row else float(row["zcmb"])  # Pantheon often uses zcmb
            mu = float(row["mu"]) if "mu" in row else float(row["mures"])  # fallback to residual-like
            s = float(row.get("sigma_mu", row.get("sigmu", row.get("muerr", 0.1))))
            w.writerow({"z": z, "mu": mu, "sigma_mu": s})


def convert_bao_dv(input_csv: Path, out_csv: Path, *, dimless_scale: float | None = None):
    """
    Convert BAO DV table with columns 'z','DV','sigma_DV'. If DV is in Mpc,
    supply dimless_scale (e.g., H0/c in chosen units) to convert; otherwise
    leave as-is for dimensionless comparison.
    """
    with input_csv.open("r", newline="", encoding="utf-8") as fi, out_csv.open("w", newline="", encoding="utf-8") as fo:
        r = csv.DictReader(fi)
        w = csv.DictWriter(fo, fieldnames=["z", "DV", "sigma_DV"])
        w.writeheader()
        for row in r:
            z = float(row["z"]) if "z" in row else float(row["zeff"])  # many BAO tables use zeff
            dv = float(row["DV"]) if "DV" in row else float(row["dv"])  # normalize cases
            sdv = float(row.get("sigma_DV", row.get("sdv", 0.05*dv)))
            if dimless_scale is not None:
                dv *= dimless_scale
                sdv *= dimless_scale
            w.writerow({"z": z, "DV": dv, "sigma_DV": sdv})


def convert_forward_events(input_csv: Path, out_json: Path):
    """
    Convert strong-field event table to forward_event_catalog JSON.

    Input CSV unified schema:
      name,observable,target,sigma,mass_msun,distance_m,source

    Required by observable:
    - ring_uas: mass_msun, distance_m
    """
    events = []
    with input_csv.open("r", newline="", encoding="utf-8") as fi:
        r = csv.DictReader(fi)
        for i, row in enumerate(r):
            name = str(row.get("name", "")).strip()
            obs = str(row.get("observable", "")).strip()
            if not name:
                raise ValueError(f"Row #{i} missing name")
            if obs != "ring_uas":
                raise ValueError(f"Row #{i} invalid observable: {obs}")
            if "target" not in row or "sigma" not in row:
                raise ValueError(f"Row #{i} missing target/sigma")

            e = {
                "name": name,
                "observable": obs,
                "target": float(row["target"]),
                "sigma": float(row["sigma"]),
            }
            src = str(row.get("source", "")).strip()
            if src:
                e["source"] = src

            if not row.get("mass_msun") or not row.get("distance_m"):
                raise ValueError(f"Row #{i} ring_uas requires mass_msun and distance_m")
            e["mass_msun"] = float(row["mass_msun"])
            e["distance_m"] = float(row["distance_m"])

            events.append(e)

    with out_json.open("w", encoding="utf-8") as fo:
        fo.write(json.dumps({"events": events}, indent=2, ensure_ascii=False) + "\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Convert public CC/SN/BAO to dimless CSVs")
    p.add_argument("--cc", type=str, default=None, help="Path to CC table with z,H,sigma_H in km/s/Mpc")
    p.add_argument("--sn", type=str, default=None, help="Path to SN table (Pantheon-like) with z,mu,sigma_mu")
    p.add_argument("--bao", type=str, default=None, help="Path to BAO DV table with z,DV,sigma_DV")
    p.add_argument("--forward_events", type=str, default=None, help="Path to strong-field events CSV for forward_event_catalog conversion")
    p.add_argument("--H0_km_s_Mpc", type=float, default=70.0, help="H0 for converting CC to dimless")
    p.add_argument("--bao_scale_dimless", type=float, default=None, help="If BAO DV in Mpc, multiply by this to make dimensionless (e.g., H0/c)")
    p.add_argument("--outdir", type=str, required=True)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    outdir = Path(args.outdir)
    _ensure_outdir(outdir)
    if args.cc:
        convert_cc(Path(args.cc), outdir / "cc_dimless.csv", H0_km_s_Mpc=args.H0_km_s_Mpc)
    if args.sn:
        convert_sn(Path(args.sn), outdir / "sn_dimless.csv")
    if args.bao:
        convert_bao_dv(Path(args.bao), outdir / "bao_dimless.csv", dimless_scale=args.bao_scale_dimless)
    if args.forward_events:
        convert_forward_events(Path(args.forward_events), outdir / "forward_event_catalog.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
