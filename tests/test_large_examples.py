"""Memory-intensive regressions for the two showcase families."""

import pytest

from fastop import spaces


pytestmark = pytest.mark.large


def test_fifth_symmetric_power_of_torus_at_prime_five():
    model = spaces.symmetric_product_of_curve(genus=1, power=5)
    cohomology = model.cohomology(p=5)

    assert sum(model.f_vector()) == 1_797_894
    assert cohomology.betti_number(2) == 2
    assert cohomology.betti_number(10) == 1
    assert cohomology.operation_rank(2, 1) == 1


def test_eleven_dimensional_lens_space_at_prime_five():
    model = spaces.lens_space(dimension=11, order=5)
    cohomology = model.cohomology(p=5)

    assert model.f_vector() == (
        6,
        81,
        650,
        3450,
        12750,
        33625,
        63750,
        86250,
        81250,
        50625,
        18750,
        3125,
    )
    assert cohomology.operation_rank(1, 0, bockstein=True) == 1
    assert cohomology.operation_rank(2, 1) == 1
