"""Benchmark symmetric cubes of closed orientable surfaces."""

from __future__ import annotations

import argparse
import sys
import time
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
if (ROOT.parent / "oddp").exists():
    sys.path.insert(0, str(ROOT.parent / "oddp"))

from fastop import spaces  # noqa: E402


def expected_betti_numbers(genus: int) -> tuple[int, ...]:
    """Return Macdonald's Betti numbers for ``Sym^3(Sigma_genus)``."""
    twice_genus = 2 * genus
    return (
        1,
        twice_genus,
        comb(twice_genus, 2) + 1,
        comb(twice_genus, 3) + twice_genus,
        comb(twice_genus, 2) + 1,
        twice_genus,
        1,
    )


def run(genus: int, *, full_cohomology: bool) -> None:
    start = time.perf_counter()
    model = spaces.symmetric_product_of_surface(genus)
    build_seconds = time.perf_counter() - start
    expected = expected_betti_numbers(genus)

    print(f"\nSym^3(Sigma_{genus})")
    print(f"  cells: {sum(model.f_vector()):,} {model.f_vector()}")
    print(f"  build: {build_seconds:.6f}s")
    print(f"  expected Betti numbers: {expected}")

    if full_cohomology:
        cohomology = model.cohomology(p=3)
        start = time.perf_counter()
        measured = tuple(
            cohomology.betti_number(degree)
            for degree in range(model.dimension + 1)
        )
        cohomology_seconds = time.perf_counter() - start
        assert measured == expected
        print(f"  full cohomology: {cohomology_seconds:.6f}s")
        print(f"  measured Betti numbers: {measured}")
    else:
        print("  full cohomology: skipped")

    cohomology = model.cohomology(p=3)
    start = time.perf_counter()
    rank = cohomology.operation_rank(2, 1)
    operation_seconds = time.perf_counter() - start
    assert rank == 1
    print(f"  P1: H2 -> H6: rank {rank} in {operation_seconds:.6f}s")
    print(
        "  measured source/target Betti numbers: "
        f"({cohomology.betti_number(2)}, {cohomology.betti_number(6)})"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("genus", nargs="*", type=int, default=(0, 1, 2))
    parser.add_argument(
        "--full-cohomology",
        action="store_true",
        help="also compute every Betti number; expensive at genus three",
    )
    args = parser.parse_args()

    for genus in args.genus:
        if genus < 0:
            parser.error("genus must be nonnegative")
        run(
            genus,
            full_cohomology=args.full_cohomology or genus <= 2,
        )


if __name__ == "__main__":
    main()
