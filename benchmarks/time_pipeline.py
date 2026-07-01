"""Time cohomology-operation pipeline stages on representative examples."""

from __future__ import annotations

import argparse
import gc
import statistics
import sys
import time
import timeit
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
if (ROOT.parent / "oddp").exists():
    sys.path.insert(0, str(ROOT.parent / "oddp"))

from fastop import spaces  # noqa: E402
from fastop._cochain_evaluation import (  # noqa: E402
    evaluate_all_targets,
    evaluate_source_focused,
    evaluate_source_mod_3,
)
from fastop._universal import universal_operation  # noqa: E402


def operation_kwargs(p: int, r: int, source_degree: int, bockstein: bool = False):
    missing = 2 * r * (p - 1) + int(bockstein)
    return {
        "p": p,
        "r": r,
        "source_degree": source_degree,
        "bockstein": bockstein,
        "target_degree": source_degree + missing,
        "missing_vertices_per_factor": missing,
        "oddp_s": -r,
        "oddp_q": -source_degree,
    }


def median_time(function, repeats: int):
    values = []
    output = None
    for _ in range(repeats):
        gc.collect()
        start = time.perf_counter()
        output = function()
        values.append(time.perf_counter() - start)
    return statistics.median(values), min(values), output


def autorange(function):
    count, total = timeit.Timer(function).autorange()
    return total / count


def percentage(part: float, total: float) -> float:
    if total == 0:
        return 0.0
    return 100 * part / total


def report_case(label, complex_, *, p: int, degree: int, r: int, bockstein: bool, repeats: int):
    print(f"\n== {label} ==")
    cohomology_median, cohomology_best, cohomology = median_time(
        lambda: complex_.cohomology(p=p),
        repeats,
    )
    element = cohomology.basis(degree)[0]
    cocycle_median, _, cochain = median_time(lambda: element.cocycle(), repeats)
    universal_median, _, universal = median_time(
        lambda: universal_operation(**operation_kwargs(p, r, degree, bockstein)),
        repeats,
    )
    signature_median, _, signatures = median_time(lambda: universal.signature_table(), repeats)

    target_faces = complex_.faces(universal.target_degree)
    target_data = cohomology._degree_data[universal.target_degree]
    print(
        "sizes "
        f"support={len(cochain)} "
        f"target_faces={len(target_faces)} "
        f"terms={len(universal.terms)} "
        f"patterns={len(signatures.coefficients)}"
    )
    print(
        f"cohomology median {cohomology_median:.8f}s "
        f"best {cohomology_best:.8f}s"
    )

    evaluators = {
        "all_targets": lambda: evaluate_all_targets(target_faces, cochain, universal),
        "source_mod_3": lambda: evaluate_source_mod_3(target_faces, cochain, signatures),
        "source_focused": lambda: evaluate_source_focused(target_faces, cochain, signatures),
    }
    for name, evaluator in evaluators.items():
        evaluation = autorange(evaluator)
        result = evaluator()
        target_vector = {
            target_data.face_to_index[simplex]: coefficient
            for simplex, coefficient in result.items()
        }
        projection = autorange(
            lambda target_vector=target_vector: cohomology.project_cocycle(
                universal.target_degree,
                target_vector,
            )
        )
        warm = cocycle_median + universal_median + signature_median + evaluation + projection
        cold = cohomology_median + warm
        print(f"{name}: warm {warm:.8f}s; cold incl cohomology {cold:.8f}s")
        print(
            f"  cohomology cold share {percentage(cohomology_median, cold):6.2f}%"
        )
        print(
            "  warm split: "
            f"cocycle {percentage(cocycle_median, warm):5.1f}% "
            f"universal {percentage(universal_median, warm):5.1f}% "
            f"signature {percentage(signature_median, warm):5.1f}% "
            f"evaluation {percentage(evaluation, warm):5.1f}% "
            f"projection {percentage(projection, warm):5.1f}%"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=9)
    args = parser.parse_args()

    report_case(
        "CP3 P1",
        spaces.complex_projective_space(3),
        p=3,
        degree=2,
        r=1,
        bockstein=False,
        repeats=args.repeats,
    )
    report_case(
        "Moore beta",
        spaces.moore_space(3),
        p=3,
        degree=1,
        r=0,
        bockstein=True,
        repeats=args.repeats,
    )


if __name__ == "__main__":
    main()
