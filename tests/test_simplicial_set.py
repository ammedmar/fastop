import pytest

from fastop import DeltaComplex, SimplexReference, SimplicialSet, spaces


def test_one_vertex_sphere_uses_normalized_chains():
    sphere = SimplicialSet.sphere(2)

    assert sphere.f_vector() == (1, 0, 1)
    assert sphere.face(2, 0, 0) is None
    assert sphere.cohomology(p=3).betti_numbers() == {0: 1, 2: 1}


def test_a_degenerate_face_can_have_a_nondegenerate_later_face():
    sphere = SimplicialSet.sphere(2)
    top = sphere.nondegenerate(2, 0)
    degenerate_edge = sphere.face_reference(top, 0)

    assert not degenerate_edge.is_nondegenerate
    assert sphere.face_reference(degenerate_edge, 0) == sphere.nondegenerate(0, 0)


def test_delta_complex_conversion_freely_adds_degeneracies():
    circle = DeltaComplex([
        [()],
        [(0, 0)],
    ])
    simplicial_circle = SimplicialSet.from_delta_complex(circle)

    assert simplicial_circle.f_vector() == (1, 1)
    assert simplicial_circle.cohomology(p=3).betti_numbers() == {0: 1, 1: 1}


def test_sage_adapter_retains_degenerate_faces():
    class FakeSimplex:
        def __init__(self, dimension, *, underlying=None, degeneracies=()):
            self._dimension = dimension
            self._underlying = underlying or self
            self._degeneracies = degeneracies

        def dimension(self):
            return self._dimension

        def nondegenerate(self):
            return self._underlying

        def degeneracies(self):
            return self._degeneracies

    vertex = FakeSimplex(0)
    top = FakeSimplex(2)
    degenerate_edge = FakeSimplex(
        1,
        underlying=vertex,
        degeneracies=(0,),
    )

    class FakeSageSimplicialSet:
        @staticmethod
        def cells():
            return {0: (vertex,), 1: (), 2: (top,)}

        @staticmethod
        def faces(simplex):
            assert simplex is top
            return (degenerate_edge,) * 3

    sphere = SimplicialSet.from_sage(FakeSageSimplicialSet())

    assert sphere.f_vector() == (1, 0, 1)
    assert sphere.face_maps[2][0][0].operator == (0, 0)
    assert sphere.cohomology(p=3).betti_numbers() == {0: 1, 2: 1}


def test_cartesian_cube_of_one_vertex_sphere_has_known_symmetric_quotient():
    sphere = SimplicialSet.sphere(2)
    cube = sphere.cartesian_product(sphere, sphere)
    symmetric_cube = sphere.symmetric_power(3)

    assert cube.f_vector() == (1, 0, 7, 42, 132, 180, 90)
    assert symmetric_cube.f_vector() == (1, 0, 3, 10, 25, 30, 15)
    assert symmetric_cube.cohomology(p=3).betti_numbers() == {0: 1, 2: 1, 4: 1, 6: 1}
    assert symmetric_cube.cohomology(p=3).operation_rank(2, 1) == 1


def test_simplicial_set_validates_references_and_identities():
    with pytest.raises(ValueError, match="at least one vertex"):
        SimplicialSet([])
    with pytest.raises(ValueError, match="wrong dimension"):
        SimplicialSet([
            [()],
            [(SimplexReference(0, 0, (0, 0)),) * 2],
        ])
    with pytest.raises(ValueError, match="surjective"):
        SimplicialSet([
            [()],
            [(SimplexReference(0, 0, (1,)),) * 2],
        ])


def test_symmetric_power_validates_power():
    circle = SimplicialSet.from_delta_complex(DeltaComplex([
        [()],
        [(0, 0)],
    ]))

    assert circle.symmetric_power(1) is circle
    with pytest.raises(TypeError, match="integer"):
        circle.symmetric_power(True)
    with pytest.raises(ValueError, match="positive"):
        circle.symmetric_power(0)


def test_symmetric_power_cell_counts_do_not_require_model_construction():
    sphere = spaces.orientable_surface(0)
    torus = spaces.orientable_surface(1)

    assert sphere.symmetric_power_f_vector(5) == (
        1,
        0,
        5,
        40,
        271,
        1197,
        3381,
        5985,
        6405,
        3780,
        945,
    )
    assert sum(torus.symmetric_power_f_vector(5)) == 1_797_894
    assert sum(
        spaces.orientable_surface(2).symmetric_power_f_vector(5)
    ) == 414_092_094


def test_symmetric_power_materializes_only_requested_degrees():
    model = SimplicialSet.sphere(2).symmetric_power(7)

    assert model._labels == {}
    assert sum(model.f_vector()) == 13_478_264
    assert model._labels == {}
    assert len(model.cells(2)) == 7
    assert set(model._labels) == {2}


@pytest.mark.parametrize(
    ("crosscaps", "f_vector"),
    [(1, (1, 1, 1)), (2, (1, 3, 2)), (3, (1, 6, 4))],
)
def test_compact_nonorientable_surfaces(crosscaps, f_vector):
    surface = spaces.nonorientable_surface(crosscaps)

    assert surface.f_vector() == f_vector
    assert surface.cohomology(p=2).betti_numbers() == {
        0: 1,
        1: crosscaps,
        2: 1,
    }
    expected_odd = {0: 1}
    if crosscaps > 1:
        expected_odd[1] = crosscaps - 1
    assert surface.cohomology(p=3).betti_numbers() == expected_odd


def test_symmetric_cube_of_projective_plane_is_rp6():
    model = spaces.nonorientable_surface(1).symmetric_power(3)
    mod_two = model.cohomology(p=2)
    mod_three = model.cohomology(p=3)

    assert model.f_vector() == (1, 3, 13, 35, 55, 45, 15)
    assert mod_two.betti_numbers() == {degree: 1 for degree in range(7)}
    assert mod_two.operation_rank(1, 1) == 1
    assert mod_two.operation_rank(2, 2) == 1
    assert mod_three.betti_numbers() == {0: 1}
    assert mod_three.operation_rank(2, 1) == 0


@pytest.mark.parametrize(
    ("genus", "f_vector", "betti"),
    [
        (0, (1, 0, 1), {0: 1, 2: 1}),
        (1, (1, 3, 2), {0: 1, 1: 2, 2: 1}),
        (2, (1, 9, 6), {0: 1, 1: 4, 2: 1}),
        (3, (1, 15, 10), {0: 1, 1: 6, 2: 1}),
    ],
)
def test_compact_surface_models(genus, f_vector, betti):
    surface = spaces.orientable_surface(genus)

    assert surface.f_vector() == f_vector
    assert surface.cohomology(p=3).betti_numbers() == betti


def test_orientable_surface_validates_genus():
    with pytest.raises(TypeError, match="integer"):
        spaces.orientable_surface(True)
    with pytest.raises(ValueError, match="nonnegative"):
        spaces.orientable_surface(-1)


def test_symmetric_cube_of_torus_is_a_small_six_manifold_with_nonzero_p1():
    symmetric_cube = spaces.orientable_surface(1).symmetric_power(3)
    cohomology = symmetric_cube.cohomology(p=3)

    assert symmetric_cube.f_vector() == (1, 19, 126, 380, 572, 420, 120)
    assert cohomology.betti_numbers() == {
        0: 1,
        1: 2,
        2: 2,
        3: 2,
        4: 2,
        5: 2,
        6: 1,
    }
    assert cohomology.operation_rank(2, 1) == 1


def test_symmetric_cube_of_genus_two_surface_matches_macdonald_betti_numbers():
    symmetric_cube = spaces.orientable_surface(2).symmetric_power(3)
    cohomology = symmetric_cube.cohomology(p=3)

    assert symmetric_cube.f_vector() == (
        1,
        219,
        2486,
        9180,
        15012,
        11340,
        3240,
    )
    assert cohomology.betti_numbers() == {
        0: 1,
        1: 4,
        2: 7,
        3: 8,
        4: 7,
        5: 4,
        6: 1,
    }
    assert cohomology.operation_rank(2, 1) == 1


def test_symmetric_fifth_power_of_sphere_is_cp5_at_prime_five():
    symmetric_fifth = spaces.orientable_surface(0).symmetric_power(5)
    cohomology = symmetric_fifth.cohomology(p=5)

    assert symmetric_fifth.f_vector() == (
        1,
        0,
        5,
        40,
        271,
        1197,
        3381,
        5985,
        6405,
        3780,
        945,
    )
    assert cohomology.betti_numbers() == {
        0: 1,
        2: 1,
        4: 1,
        6: 1,
        8: 1,
        10: 1,
    }
    assert cohomology.operation_rank(2, 1) == 1
