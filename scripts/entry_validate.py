# -*- coding: utf-8 -*-
"""Single entrypoint for full model validation."""
from __future__ import annotations

import subprocess
import sys


def main() -> int:
    cmd = [sys.executable, "-m", "scripts.run_validation_suite"]
    print("[entry_validate] Starting validation via python -m scripts.run_validation_suite", flush=True)
    proc = subprocess.run(cmd)
    print(f"[entry_validate] Validation finished with rc={proc.returncode}", flush=True)
    return int(proc.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
