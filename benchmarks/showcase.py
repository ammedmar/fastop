"""Benchmark the compact odd-primary showcase models."""

from __future__ import annotations

import argparse
import gc
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
if (ROOT.parent / "oddp").exists():
    sys.path.insert(0, str(ROOT.parent / "oddp"))

from fastop import spaces  # noqa: E402


def median_time(function, repeats: int):
    samples = []
    output = None
    for _ in range(repeats):
        gc.collect()
        start = time.perf_counter()
        output = function()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples), output


def total_cells(model) -> int:
    if hasattr(model, "f_vector"):
        return sum(model.f_vector())
    return sum(
        len(model.cells(degree))
        for degree in range(model.dimension + 1)
    )


def report(label, build, operations, repeats: int) -> None:
    build_seconds, model = median_time(build, repeats)

    def compute_cohomology():
        cohomology = model.cohomology(p=3)
        return cohomology, cohomology.betti_numbers()

    cohomology_seconds, (cohomology, betti) = median_time(
        compute_cohomology,
        repeats,
    )
    print(f"\n{label}")
    print(f"  dimension: {model.dimension}")
    print(f"  cells: {total_cells(model):,}")
    print(f"  build median: {build_seconds:.6f}s")
    print(f"  cohomology median: {cohomology_seconds:.6f}s")
    print(f"  mod-3 Betti numbers: {betti}")
    for operation_label, degree, r, bockstein in operations:
        rank = cohomology.operation_rank(
            degree,
            r,
            bockstein=bockstein,
        )
        operation_seconds, measured_rank = median_time(
            lambda: cohomology.operation_rank(
                degree,
                r,
                bockstein=bockstein,
            ),
            repeats,
        )
        assert measured_rank == rank
        print(
            f"  {operation_label}: rank {rank}, "
            f"warm median {operation_seconds:.6f}s"
        )


def cases():
    cp3 = spaces.complex_projective_space(3)
    return {
        "suspension-1": (
            "Suspension Sigma CP3",
            lambda: cp3.suspension(),
            (("P1: H3 -> H7", 3, 1, False),),
        ),
        "suspension-2": (
            "Double suspension Sigma^2 CP3",
            lambda: cp3.suspension(2),
            (("P1: H4 -> H8", 4, 1, False),),
        ),
        "lens": (
            "Compact lens space L7(3)",
            lambda: spaces.lens_space(7, 3),
            (
                ("beta P0: H1 -> H2", 1, 0, True),
                ("P1: H2 -> H6", 2, 1, False),
            ),
        ),
        "symmetric-sphere": (
            "Symmetric cube Sym^3(S2) = CP3",
            lambda: spaces.minimal_simplicial_sphere(2).symmetric_power(3),
            (("P1: H2 -> H6", 2, 1, False),),
        ),
        "symmetric-torus": (
            "Symmetric cube Sym^3(T2)",
            lambda: spaces.symmetric_product_of_torus(3),
            (("P1: H2 -> H6", 2, 1, False),),
        ),
    }


def main() -> None:
    available = cases()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "case",
        nargs="*",
        choices=tuple(available),
        default=tuple(available),
    )
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()

    for name in args.case:
        report(*available[name], repeats=args.repeats)


if __name__ == "__main__":
    main()
