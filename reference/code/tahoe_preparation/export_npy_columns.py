#!/usr/bin/env python3
"""Memory-map a C-order NPY matrix and export selected columns.

This is a mechanical plotting-data loader. It does not transform, normalize,
fit, or otherwise alter the scientific data.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("columns", type=Path,
                        help="Text file containing zero-based column indexes")
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()

    matrix = np.load(args.source, mmap_mode="r")
    columns = np.loadtxt(args.columns, dtype=np.int64, ndmin=1)
    if matrix.ndim != 2:
        raise ValueError(f"Expected a two-dimensional NPY matrix: {args.source}")
    if columns.ndim != 1 or len(np.unique(columns)) != len(columns):
        raise ValueError("Column indexes must be a unique one-dimensional vector")
    if np.any(columns < 0) or np.any(columns >= matrix.shape[1]):
        raise ValueError("A requested column is outside the source matrix")
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.destination, np.asarray(matrix[:, columns], dtype=np.float32))


if __name__ == "__main__":
    main()
