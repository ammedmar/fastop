import sys
from pathlib import Path

import pytest

from fastop import spaces
from fastop._cochain_evaluation import cochain_operation_vector_from_universal
from fastop._oddp_bridge import (
    cochain_operation_vector_oddp,
)
from fastop._precomputed_universal import precomputed_terms
from fastop._universal import (
    native_universal_operation,
    oddp_universal_operation,
    universal_operation,
)


ODDP_ROOT = Path(__file__).resolve().parents[2] / "oddp"
if ODDP_ROOT.exists():
    sys.path.insert(0, str(ODDP_ROOT))

pytest.importorskip("oddp")


def _universal_kwargs(p, r, source_degree, bockstein=False):
    missing = 2 * r * (p - 1) + int(bockstein)
    return {
        "p": p,
        "r": r,
        "source_degree": source_degree,
        "bockstein": bockstein,
        "target_degree": source_degree + missing,
        "missing_vertices_per_factor": missing,
    }


def _bridge_kwargs(p, r, source_degree, bockstein=False):
    return {
        **_universal_kwargs(p, r, source_degree, bockstein),
        "oddp_s": -r,
        "oddp_q": -source_degree,
    }


def test_source_mod_3_operation_on_cp3_generator():
    cohomology = spaces.complex_projective_space(3).cohomology(p=3)
    x = cohomology.basis(2)[0]

    assert cohomology.betti_numbers() == {0: 1, 2: 1, 4: 1, 6: 1}
    assert x.operation(1, algorithm="source_mod_3") == cohomology.basis(6)[0]
    assert cohomology.operation_rank(2, 1, algorithm="source_mod_3") == 1


def test_source_mod_3_operation_with_bockstein_on_moore_space_generator():
    cohomology = spaces.moore_space(3).cohomology(p=3)
    x = cohomology.basis(1)[0]

    assert cohomology.betti_numbers() == {0: 1, 1: 1, 2: 1}
    assert x.operation(0, bockstein=True, algorithm="source_mod_3") == cohomology.basis(2)[0]
    assert cohomology.operation_rank(1, 0, bockstein=True, algorithm="source_mod_3") == 1


def test_reference_universal_operation_for_cp3_generator():
    universal = universal_operation(**_bridge_kwargs(p=3, r=1, source_degree=2))

    assert universal.target_degree == 6
    assert universal.terms == {((0, 1, 2), (2, 3, 4), (4, 5, 6)): 2}


def test_reference_universal_operation_for_moore_bockstein():
    universal = universal_operation(
        **_bridge_kwargs(p=3, r=0, source_degree=1, bockstein=True)
    )

    assert universal.target_degree == 2
    assert universal.terms == {((0, 2), (0, 1), (1, 2)): 2}


@pytest.mark.parametrize("p,r,source_degree", [(3, 1, 2), (5, 1, 2), (3, 2, 4)])
def test_native_universal_top_reduced_power_matches_oddp_oracle(
    p,
    r,
    source_degree,
):
    operation = _bridge_kwargs(p=p, r=r, source_degree=source_degree)

    assert native_universal_operation(**_universal_kwargs(p, r, source_degree)) == (
        oddp_universal_operation(**operation)
    )


def test_oddp_universal_oracle_covers_unimplemented_bockstein_slice():
    operation = _bridge_kwargs(p=3, r=0, source_degree=1, bockstein=True)

    assert native_universal_operation(
        **_universal_kwargs(p=3, r=0, source_degree=1, bockstein=True)
    ) is None
    assert oddp_universal_operation(**operation).terms == {
        ((0, 2), (0, 1), (1, 2)): 2
    }


@pytest.mark.parametrize(
    "p,r,source_degree,bockstein",
    [
        (3, 0, 1, True),
        (3, 0, 2, True),
        (5, 0, 1, True),
        (3, 1, 3, False),
        (3, 1, 4, False),
    ],
)
def test_precomputed_universal_terms_match_oddp_oracle(
    p,
    r,
    source_degree,
    bockstein,
):
    operation = _bridge_kwargs(
        p=p,
        r=r,
        source_degree=source_degree,
        bockstein=bockstein,
    )

    assert precomputed_terms(
        p=p,
        r=r,
        source_degree=source_degree,
        bockstein=bockstein,
    ) == oddp_universal_operation(**operation).terms


def test_all_targets_evaluator_matches_oddp_direct_on_cp3_generator():
    complex_ = spaces.complex_projective_space(3)
    cohomology = complex_.cohomology(p=3)
    operation = _bridge_kwargs(p=3, r=1, source_degree=2)
    target_data = cohomology._degree_data[operation["target_degree"]]
    cochain = cohomology.basis(2)[0].cocycle()
    universal = universal_operation(**operation)

    assert cochain_operation_vector_from_universal(
        complex_,
        cochain,
        universal,
        target_data.face_to_index,
        algorithm="all_targets",
    ) == cochain_operation_vector_oddp(
        complex_,
        cochain,
        p=operation["p"],
        bockstein=operation["bockstein"],
        oddp_s=operation["oddp_s"],
        oddp_q=operation["oddp_q"],
        target_face_to_index=target_data.face_to_index,
        algorithm="direct",
    )


def test_all_targets_evaluator_matches_oddp_direct_on_moore_bockstein():
    complex_ = spaces.moore_space(3)
    cohomology = complex_.cohomology(p=3)
    operation = _bridge_kwargs(p=3, r=0, source_degree=1, bockstein=True)
    target_data = cohomology._degree_data[operation["target_degree"]]
    cochain = cohomology.basis(1)[0].cocycle()
    universal = universal_operation(**operation)

    assert cochain_operation_vector_from_universal(
        complex_,
        cochain,
        universal,
        target_data.face_to_index,
        algorithm="all_targets",
    ) == cochain_operation_vector_oddp(
        complex_,
        cochain,
        p=operation["p"],
        bockstein=operation["bockstein"],
        oddp_s=operation["oddp_s"],
        oddp_q=operation["oddp_q"],
        target_face_to_index=target_data.face_to_index,
        algorithm="direct",
    )
