"""Benchmark the first reduced power at p=7 on CP^7."""

from __future__ import annotations

import resource
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fastop import spaces  # noqa: E402


def peak_memory_mib() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return usage / (1024 * 1024) if sys.platform == "darwin" else usage / 1024


def main() -> None:
    start = time.perf_counter()
    model = spaces.symmetric_product_of_surface(0, power=7)
    build_seconds = time.perf_counter() - start
    predicted = model.f_vector()
    print("Sym^7(S^2) = CP^7")
    print(f"  predicted cells: {sum(predicted):,} {predicted}")
    print(f"  lazy model: {build_seconds:.6f}s")

    cohomology = model.cohomology(p=7)
    start = time.perf_counter()
    source_betti = cohomology.betti_number(2)
    source_seconds = time.perf_counter() - start
    start = time.perf_counter()
    target_betti = cohomology.betti_number(14)
    target_seconds = time.perf_counter() - start
    start = time.perf_counter()
    rank = cohomology.operation_rank(2, 1)
    operation_seconds = time.perf_counter() - start

    assert (source_betti, target_betti, rank) == (1, 1, 1)
    print(f"  b2 = 1: {source_seconds:.6f}s")
    print(f"  b14 = 1: {target_seconds:.6f}s")
    print(f"  P1: H2 -> H14: rank 1 in {operation_seconds:.6f}s")
    print(f"  peak memory: {peak_memory_mib():.1f} MiB")


if __name__ == "__main__":
    main()
