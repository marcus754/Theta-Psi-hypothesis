# -*- coding: utf-8 -*-
"""Thin entrypoint for observational head-to-head comparison."""
from __future__ import annotations

from typing import Sequence

from observational.compare_models import main, parse_args


def run(argv: Sequence[str] | None = None) -> int:
    return main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
