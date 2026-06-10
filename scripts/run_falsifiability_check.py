# -*- coding: utf-8 -*-
"""Thin entrypoint for the observational falsifiability layer."""
from __future__ import annotations

from typing import Sequence

from observational import falsifiability as _impl

parse_args = _impl.parse_args


def main(argv: Sequence[str] | None = None) -> int:
    return _impl.main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
