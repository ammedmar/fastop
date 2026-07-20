import pytest

from fastop import DeltaComplex, SimplexReference, SimplicialSet, spaces


def test_minimal_sphere_uses_normalized_chains():
    sphere = SimplicialSet.minimal_sphere(2)

    assert sphere.f_vector() == (1, 0, 1)
    assert sphere.face(2, 0, 0) is None
    assert sphere.cohomology(p=3).betti_numbers() == {0: 1, 2: 1}


def test_a_degenerate_face_can_have_a_nondegenerate_later_face():
    sphere = SimplicialSet.minimal_sphere(2)
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


def test_cartesian_cube_of_minimal_sphere_has_known_symmetric_quotient():
    sphere = SimplicialSet.minimal_sphere(2)
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


def test_symmetric_cube_of_torus_is_a_small_six_manifold_with_nonzero_p1():
    symmetric_cube = spaces.symmetric_product_of_torus(3)
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
