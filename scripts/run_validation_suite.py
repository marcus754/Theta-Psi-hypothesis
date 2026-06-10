# -*- coding: utf-8 -*-
"""Thin entrypoint for the observational validation layer."""
from __future__ import annotations

from typing import Sequence

from observational import validation as _impl

DIAGNOSTIC_STEPS = _impl.DIAGNOSTIC_STEPS
EMPIRICAL_STEPS = _impl.EMPIRICAL_STEPS
build_steps = _impl.build_steps
parse_args = _impl.parse_args
_run = _impl._run


def main(argv: Sequence[str] | None = None) -> int:
    prev = _impl._run
    _impl._run = _run
    try:
        return _impl.main(argv)
    finally:
        _impl._run = prev


if __name__ == "__main__":
    raise SystemExit(main())
