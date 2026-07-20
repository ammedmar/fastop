"""Benchmark P1 on fifth symmetric products of surfaces."""

from __future__ import annotations

import argparse
import resource
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
if (ROOT.parent / "oddp").exists():
    sys.path.insert(0, str(ROOT.parent / "oddp"))

from fastop import spaces  # noqa: E402


def peak_memory_mib() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return usage / (1024 * 1024)
    return usage / 1024


def run(genus: int, *, max_cells: int) -> None:
    surface = spaces.minimal_simplicial_surface(genus)
    predicted = surface.symmetric_power_f_vector(5)
    total_cells = sum(predicted)

    print(f"\nSym^5(Sigma_{genus})")
    print(f"  predicted cells: {total_cells:,} {predicted}")
    if total_cells > max_cells:
        print(f"  skipped: exceeds the {max_cells:,}-cell limit")
        return

    start = time.perf_counter()
    model = surface.symmetric_power(5)
    build_seconds = time.perf_counter() - start
    assert model.f_vector() == predicted

    cohomology = model.cohomology(p=5)
    start = time.perf_counter()
    source_betti = cohomology.betti_number(2)
    source_seconds = time.perf_counter() - start
    start = time.perf_counter()
    target_betti = cohomology.betti_number(10)
    target_seconds = time.perf_counter() - start
    start = time.perf_counter()
    rank = cohomology.operation_rank(2, 1)
    operation_seconds = time.perf_counter() - start

    assert source_betti == (1 if genus == 0 else 2)
    assert target_betti == 1
    assert rank == 1
    print(f"  build: {build_seconds:.6f}s")
    print(f"  b2 = {source_betti}: {source_seconds:.6f}s")
    print(f"  b10 = {target_betti}: {target_seconds:.6f}s")
    print(f"  P1: H2 -> H10: rank {rank} in {operation_seconds:.6f}s")
    print(f"  peak memory: {peak_memory_mib():.1f} MiB")


def run_lens_space() -> None:
    print("\nL^11(5)")
    start = time.perf_counter()
    model = spaces.lens_space(11, 5)
    build_seconds = time.perf_counter() - start
    cohomology = model.cohomology(p=5)

    start = time.perf_counter()
    bockstein_rank = cohomology.operation_rank(1, 0, bockstein=True)
    bockstein_seconds = time.perf_counter() - start
    start = time.perf_counter()
    reduced_power_rank = cohomology.operation_rank(2, 1)
    reduced_power_seconds = time.perf_counter() - start

    assert bockstein_rank == 1
    assert reduced_power_rank == 1
    print(f"  cells: {sum(model.f_vector()):,} {model.f_vector()}")
    print(f"  build: {build_seconds:.6f}s")
    print(
        "  beta P0: H1 -> H2: "
        f"rank {bockstein_rank} in {bockstein_seconds:.6f}s"
    )
    print(
        "  P1: H2 -> H10: "
        f"rank {reduced_power_rank} in {reduced_power_seconds:.6f}s"
    )
    print(f"  peak memory: {peak_memory_mib():.1f} MiB")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("genus", nargs="*", type=int, default=(0, 1))
    parser.add_argument("--max-cells", type=int, default=5_000_000)
    parser.add_argument(
        "--lens",
        action="store_true",
        help="also run the larger L^11(5) face-map validation",
    )
    args = parser.parse_args()

    for genus in args.genus:
        if genus < 0:
            parser.error("genus must be nonnegative")
        run(genus, max_cells=args.max_cells)
    if args.lens:
        run_lens_space()


if __name__ == "__main__":
    main()
