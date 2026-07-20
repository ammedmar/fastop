import pytest

from fastop import DeltaComplex, FiniteGroupAction, spaces


def test_face_map_circle_keeps_repeated_faces_and_cancels_its_boundary():
    circle = DeltaComplex([
        [()],
        [(0, 0)],
    ])

    assert circle.dimension == 1
    assert circle.f_vector() == (1, 1)
    assert circle.cohomology(p=2).betti_numbers() == {0: 1, 1: 1}
    assert circle.cohomology(p=3).betti_numbers() == {0: 1, 1: 1}


def test_one_vertex_torus_has_the_expected_odd_prime_cohomology():
    torus = DeltaComplex([
        [()],
        [(0, 0), (0, 0), (0, 0)],
        [(2, 1, 0), (0, 1, 2)],
    ])

    assert torus.f_vector() == (1, 3, 2)
    assert torus.cohomology(p=3).betti_numbers() == {0: 1, 1: 2, 2: 1}


def test_abstract_complex_conversion_preserves_cohomology_and_operations():
    rp2 = spaces.real_projective_plane()
    cp3 = spaces.complex_projective_space(3)
    delta_rp2 = rp2.as_delta_complex()
    delta_cp3 = cp3.as_delta_complex()

    assert delta_rp2.f_vector() == tuple(
        len(rp2.faces(degree)) for degree in range(rp2.dimension + 1)
    )
    assert delta_rp2.cohomology().betti_numbers() == {0: 1, 1: 1, 2: 1}
    assert delta_rp2.cohomology().operation_rank(1, 1) == 1
    assert delta_cp3.cohomology(p=3).operation_rank(2, 1) == 1


def test_face_restrictions_support_nontrivial_odd_primary_formulas():
    suspended_cp3 = (
        spaces.complex_projective_space(3)
        .suspension()
        .as_delta_complex()
        .cohomology(p=3)
    )

    assert suspended_cp3.operation_rank(3, 1) == 1


def test_face_restriction_uses_local_positions():
    triangle = spaces.simplex(2).as_delta_complex()
    top = triangle.cells(2)[0]
    edge_indices = {
        positions: triangle.restrict(2, top, positions)
        for positions in ((0, 1), (0, 2), (1, 2))
    }

    assert len(set(edge_indices.values())) == 3


def test_sage_cells_adapter_copies_dense_face_tables():
    class FakeSageDeltaComplex:
        @staticmethod
        def cells():
            return {
                -1: ((),),
                0: ((),),
                1: ((0, 0),),
            }

    circle = DeltaComplex.from_sage(FakeSageDeltaComplex())

    assert circle.face_maps == (((),), ((0, 0),))
    assert circle.cohomology(p=3).betti_numbers() == {0: 1, 1: 1}


def test_public_group_action_builds_a_cyclic_face_map_quotient():
    order = 5
    circle = DeltaComplex([
        [() for _ in range(order)],
        [((index + 1) % order, index) for index in range(order)],
    ])
    rotation = FiniteGroupAction.from_cell_maps(
        circle,
        lambda degree, cell: (cell + 1) % order,
    )

    assert rotation.order(circle) == order
    assert rotation.is_free(circle)

    quotient = circle.quotient(rotation, require_free=True)

    assert quotient.face_maps == (((),), ((0, 0),))
    assert quotient.cohomology(p=5).betti_numbers() == {0: 1, 1: 1}

    explicit_rotation = FiniteGroupAction.cyclic(rotation.generators[0])
    assert explicit_rotation.order(circle) == order


def test_public_group_action_requires_a_generator():
    with pytest.raises(ValueError, match="at least one generator"):
        FiniteGroupAction(())


def test_delta_complex_validates_face_tables():
    with pytest.raises(ValueError, match="at least one vertex"):
        DeltaComplex([])
    with pytest.raises(ValueError, match="exactly 2 faces"):
        DeltaComplex([[()], [(0,)]])
    with pytest.raises(ValueError, match="outside"):
        DeltaComplex([[()], [(0, 1)]])
    with pytest.raises(ValueError, match="identities"):
        DeltaComplex([
            [(), ()],
            [(0, 0), (0, 1), (1, 1)],
            [(0, 1, 2)],
        ])


def test_face_map_input_rejects_vertex_only_evaluators():
    cp3 = spaces.complex_projective_space(3).as_delta_complex()
    cohomology = cp3.cohomology(p=3)

    with pytest.raises(ValueError, match="abstract simplicial complex"):
        cohomology.operation_rank(2, 1, algorithm="prime-three")


def test_strict_cyclic_quotient_of_a_circle():
    circle = DeltaComplex([
        [(), (), ()],
        [(1, 0), (2, 1), (0, 2)],
    ])
    rotation = (
        (1, 2, 0),
        (1, 2, 0),
    )

    quotient = circle.quotient([rotation], require_free=True)

    assert quotient.face_maps == (((),), ((0, 0),))
    assert quotient.cohomology(p=3).betti_numbers() == {0: 1, 1: 1}


def test_quotient_validates_actions_and_optional_freeness():
    interval = DeltaComplex([
        [(), ()],
        [(1, 0)],
    ])

    with pytest.raises(ValueError, match="permutation"):
        interval.quotient([((0, 0), (0,))])
    with pytest.raises(ValueError, match="commute"):
        interval.quotient([((1, 0), (0,))])

    fixed_vertex = DeltaComplex([
        [()],
        [(0, 0), (0, 0)],
    ])
    with pytest.raises(ValueError, match="not free"):
        fixed_vertex.quotient([((0,), (1, 0))], require_free=True)


def test_compact_lens_space_detects_ground_truth_mod_three_operations():
    lens = spaces.lens_space(7, 3)
    cohomology = lens.cohomology(p=3)

    assert lens.f_vector() == (4, 22, 72, 153, 216, 198, 108, 27)
    assert cohomology.betti_numbers() == {degree: 1 for degree in range(8)}
    assert cohomology.operation_rank(1, 0, bockstein=True) == 1
    assert cohomology.operation_rank(2, 1) == 1
