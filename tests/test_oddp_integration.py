import sys
from pathlib import Path

import pytest

from fastop import spaces
from fastop._odd_primary.indices import OperationIndex
from fastop._odd_primary.reference import universal_operation


ODDP_ROOT = Path(__file__).resolve().parents[2] / "oddp"
if ODDP_ROOT.exists():
    sys.path.insert(0, str(ODDP_ROOT))

pytest.importorskip("oddp")


def test_prime_three_operation_on_cp3_generator():
    cohomology = spaces.complex_projective_space(3).cohomology(p=3)
    x = cohomology.basis(2)[0]

    assert cohomology.betti_numbers() == {0: 1, 2: 1, 4: 1, 6: 1}
    assert x.operation(1, algorithm="prime-three") == cohomology.basis(6)[0]
    assert cohomology.operation_rank(2, 1, algorithm="prime-three") == 1


def test_prime_three_operation_with_bockstein_on_moore_space_generator():
    cohomology = spaces.moore_space(3).cohomology(p=3)
    x = cohomology.basis(1)[0]

    assert cohomology.betti_numbers() == {0: 1, 1: 1, 2: 1}
    assert x.operation(0, bockstein=True, algorithm="prime-three") == cohomology.basis(2)[0]
    assert cohomology.operation_rank(1, 0, bockstein=True, algorithm="prime-three") == 1


def test_reference_universal_operation_for_cp3_generator():
    universal = universal_operation(OperationIndex(p=3, r=1, source_degree=2))

    assert universal.target_degree == 6
    assert universal.terms == {((0, 1, 2), (2, 3, 4), (4, 5, 6)): 2}


def test_reference_universal_operation_for_moore_bockstein():
    universal = universal_operation(
        OperationIndex(p=3, r=0, source_degree=1, bockstein=True)
    )

    assert universal.target_degree == 2
    assert universal.terms == {((0, 2), (0, 1), (1, 2)): 2}
