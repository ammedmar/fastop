"""Profile real-world Steenrod computations on symmetric surface powers.

The symmetric-power model is lazy, so timing its constructor alone is not a
useful measure of setup cost.  This benchmark follows the public workflow and
separates the cohomology work needed by the source and target degrees from the
actual cochain-level operation.
"""

from __future__ import annotations

import argparse
import gc
import resource
import statistics
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


def timed(function):
    start = time.perf_counter()
    output = function()
    return time.perf_counter() - start, output


def median_time(function, repeats: int):
    samples = []
    output = None
    for _ in range(repeats):
        start = time.perf_counter()
        output = function()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples), output


def peak_memory_mib() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return usage / (1024 * 1024)
    return usage / 1024


def percentage(part: float, total: float) -> float:
    return 0.0 if total == 0 else 100 * part / total


def run(genus: int, *, full_cohomology: bool, hot_repeats: int) -> None:
    expected = expected_betti_numbers(genus)
    setup_seconds, model = timed(
        lambda: spaces.orientable_surface(genus).symmetric_power(3)
    )
    cohomology_init_seconds, cohomology = timed(lambda: model.cohomology(p=3))

    print(f"\nSym^3(Sigma_{genus})", flush=True)
    print(f"  cells: {sum(model.f_vector()):,} {model.f_vector()}", flush=True)
    print(f"  expected Betti numbers: {expected}", flush=True)
    print(f"  lazy model setup: {setup_seconds:.6f}s", flush=True)
    print(
        f"  cohomology object initialization: {cohomology_init_seconds:.6f}s",
        flush=True,
    )

    source_seconds, source_betti = timed(lambda: cohomology.betti_number(2))
    print(f"  source H2 (b2 = {source_betti}): {source_seconds:.6f}s", flush=True)
    target_seconds, target_betti = timed(lambda: cohomology.betti_number(6))
    print(f"  target H6 (b6 = {target_betti}): {target_seconds:.6f}s", flush=True)
    operation_seconds, rank = timed(lambda: cohomology.operation_rank(2, 1))
    print(f"  first P1 evaluation: rank {rank} in {operation_seconds:.6f}s", flush=True)
    if hot_repeats:
        hot_operation_seconds, hot_rank = median_time(
            lambda: cohomology.operation_rank(2, 1),
            hot_repeats,
        )
    else:
        hot_operation_seconds = None
        hot_rank = rank

    assert source_betti == expected[2]
    assert target_betti == expected[6]
    assert rank == hot_rank == 1
    focused_seconds = (
        setup_seconds
        + cohomology_init_seconds
        + source_seconds
        + target_seconds
        + operation_seconds
    )

    if hot_operation_seconds is None:
        print("  cached P1 median: skipped")
    else:
        print(
            f"  cached P1 median ({hot_repeats} runs): "
            f"{hot_operation_seconds:.6f}s"
        )
    print(f"  focused cold total: {focused_seconds:.6f}s")
    print(
        "  focused cold split: "
        f"source H2 {percentage(source_seconds, focused_seconds):.1f}%, "
        f"target H6 {percentage(target_seconds, focused_seconds):.1f}%, "
        f"P1 evaluation {percentage(operation_seconds, focused_seconds):.1f}%"
    )

    if full_cohomology:
        del cohomology, model
        gc.collect()

        def compute_full_cohomology():
            full_model = spaces.orientable_surface(genus).symmetric_power(3)
            full_cohomology = full_model.cohomology(p=3)
            return tuple(
                full_cohomology.betti_number(degree)
                for degree in range(full_model.dimension + 1)
            )

        cohomology_seconds, measured = timed(compute_full_cohomology)
        assert measured == expected
        print(f"  full cold cohomology workflow: {cohomology_seconds:.6f}s")
        print(f"  measured Betti numbers: {measured}")
    else:
        print("  full cold cohomology workflow: skipped")
    print(f"  process peak RSS (cumulative): {peak_memory_mib():.1f} MiB")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("genus", nargs="*", type=int, default=(0, 1, 2))
    parser.add_argument(
        "--full-cohomology",
        action="store_true",
        help="also compute every Betti number; expensive at genus three",
    )
    parser.add_argument(
        "--hot-repeats",
        type=int,
        default=5,
        help="number of cached operation runs used for the median; zero skips them",
    )
    args = parser.parse_args()
    if args.hot_repeats < 0:
        parser.error("--hot-repeats must be nonnegative")

    for genus in args.genus:
        if genus < 0:
            parser.error("genus must be nonnegative")
        run(
            genus,
            full_cohomology=args.full_cohomology or genus <= 2,
            hot_repeats=args.hot_repeats,
        )


if __name__ == "__main__":
    main()
