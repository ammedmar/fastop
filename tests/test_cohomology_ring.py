"""Cup products in the common finite-cell cohomology interface."""

from __future__ import annotations

import pytest

from fastop import spaces


@pytest.mark.parametrize(
    "model",
    [
        spaces.complex_projective_space(3),
        spaces.complex_projective_space(3).as_delta_complex(),
        spaces.minimal_simplicial_sphere(2).symmetric_power(3),
    ],
)
def test_top_reduced_power_is_cup_power_on_every_input_model(model):
    cohomology = model.cohomology(p=3)
    generator = cohomology.basis(2)[0]

    assert generator**3 == generator.operation(1)
    assert generator**3
    assert cohomology.one() * generator == generator
    assert generator * cohomology.one() == generator


def test_projective_plane_degree_one_class_squares_to_fundamental_class():
    cohomology = spaces.nonorientable_surface(1).cohomology(p=2)
    generator = cohomology.basis(1)[0]

    assert generator * generator == cohomology.basis(2)[0]
    assert generator**2 == generator.sq(1)
    assert cohomology.cup_product_matrix(1, 1) == [{0: 1}]


def test_torus_product_is_graded_commutative_at_an_odd_prime():
    cohomology = spaces.orientable_surface(1).cohomology(p=3)
    left, right = cohomology.basis(1)

    assert left * right
    assert right * left == 2 * (left * right)
    assert (left * right) * cohomology.one() == left * right


def test_cup_product_is_associative_and_distributive():
    cohomology = spaces.complex_projective_space(3).cohomology(p=3)
    unit = cohomology.one()
    generator = cohomology.basis(2)[0]

    assert (generator * generator) * generator == generator * (generator * generator)
    assert generator * (generator + unit) == generator**2 + generator


def test_cup_product_rejects_classes_from_another_parent():
    first = spaces.sphere(1).cohomology(p=3)
    second = spaces.sphere(1).cohomology(p=3)

    with pytest.raises(ValueError, match="belong"):
        first.cup_product(first.one(), second.one())


def test_reduced_cohomology_has_no_unit():
    cohomology = spaces.sphere(2).cohomology(p=3, reduced=True)

    with pytest.raises(ValueError, match="no multiplicative unit"):
        cohomology.one()
