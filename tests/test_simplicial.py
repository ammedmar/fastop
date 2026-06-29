import sys
import types

import pytest

from fastop import SimplicialComplex, __version__, spaces


def test_version_is_development_version():
    assert __version__ == "0.1.0.dev0"


def test_complex_generates_all_faces():
    complex_ = SimplicialComplex.from_facets([(2, 0, 1)])
    assert complex_.dimension == 2
    assert complex_.vertices == (0, 1, 2)
    assert complex_.faces(0) == frozenset({(0,), (1,), (2,)})
    assert complex_.faces(1) == frozenset({(0, 1), (0, 2), (1, 2)})
    assert complex_.faces(2) == frozenset({(0, 1, 2)})


def test_rejects_invalid_facets():
    with pytest.raises(ValueError, match="at least one"):
        SimplicialComplex([])
    with pytest.raises(ValueError, match="repeated"):
        SimplicialComplex([(0, 0, 1)])
    with pytest.raises(ValueError, match="non-empty"):
        SimplicialComplex([()])


def test_cohomology_requires_prime_characteristic():
    complex_ = SimplicialComplex([(0, 1)])
    with pytest.raises(ValueError, match="prime"):
        complex_.cohomology(p=1)
    with pytest.raises(ValueError, match="prime"):
        complex_.cohomology(p=4)


def test_basic_mod_two_betti_numbers():
    interval = SimplicialComplex([(0, 1)])
    circle = SimplicialComplex([(0, 1), (0, 2), (1, 2)])
    sphere = spaces.sphere(2)

    assert interval.cohomology().betti_numbers() == {0: 1}
    assert interval.cohomology(reduced=True).betti_numbers() == {}
    assert circle.cohomology().betti_numbers() == {0: 1, 1: 1}
    assert sphere.cohomology().betti_numbers() == {0: 1, 2: 1}


def test_mod_p_betti_numbers_for_odd_primes():
    assert spaces.sphere(2).cohomology(p=3).betti_numbers() == {0: 1, 2: 1}
    assert spaces.complex_projective_plane().cohomology(p=5).betti_numbers() == {
        0: 1,
        2: 1,
        4: 1,
    }


def test_cohomology_object_exposes_prime():
    for p in (2, 3, 5):
        cohomology = spaces.complex_projective_plane().cohomology(p=p)
        basis_element = cohomology.basis(2)[0]

        assert cohomology.p == p
        assert basis_element.p == p
        assert cohomology.betti_numbers() == {0: 1, 2: 1, 4: 1}
        assert f"p={p}" in repr(cohomology)


def test_mod_p_class_arithmetic_uses_plain_int_coefficients():
    cohomology = spaces.sphere(2).cohomology(p=3)
    x = cohomology.basis(2)[0]

    assert (x + x)._coordinates == {2: {0: 2}}
    assert (x + x + x).is_zero()
    assert (2 * x)._coordinates == {2: {0: 2}}


def test_steenrod_square_spelling_is_only_available_mod_two():
    cohomology = spaces.complex_projective_plane().cohomology(p=3)
    x = cohomology.basis(2)[0]

    with pytest.raises(NotImplementedError, match="p=2"):
        x.sq(2)


def test_bockstein_operations_are_only_available_at_odd_primes():
    cohomology = spaces.complex_projective_plane().cohomology()
    x = cohomology.basis(2)[0]

    with pytest.raises(NotImplementedError, match="odd prime"):
        x.operation(1, bockstein=True)
    with pytest.raises(NotImplementedError, match="odd prime"):
        cohomology.operation_rank(2, 1, bockstein=True)


def test_odd_primary_operation_uses_oddp_and_projects(monkeypatch):
    cohomology = spaces.complex_projective_plane().cohomology(p=3)
    one = cohomology.basis(0)[0]
    top_cocycle = cohomology.basis(4)[0].cocycle()
    calls = []

    class FakeSteenrod:
        @staticmethod
        def cochain_operation(complex_, cochain, p, s, q, *, bockstein, algorithm):
            calls.append((complex_, cochain, p, s, q, bockstein, algorithm))
            return top_cocycle

    fake_oddp = types.ModuleType("oddp")
    fake_oddp.Steenrod = FakeSteenrod
    monkeypatch.setitem(sys.modules, "oddp", fake_oddp)

    assert one.operation(1, algorithm="direct") == cohomology.basis(4)[0]
    complex_, cochain, p, s, q, bockstein, algorithm = calls[0]
    assert complex_[4] == set(cohomology.complex.faces(4))
    assert cochain == one.cocycle()
    assert (p, s, q, bockstein, algorithm) == (3, -1, 0, False, "direct")
    assert cohomology.operation_rank(0, 1, algorithm="direct") == 1


def test_odd_primary_operation_uses_oddp_bockstein_conventions(monkeypatch):
    cohomology = spaces.complex_projective_plane().cohomology(p=5)
    one = cohomology.basis(0)[0]
    calls = []

    class FakeSteenrod:
        @staticmethod
        def cochain_operation(complex_, cochain, p, s, q, *, bockstein, algorithm):
            calls.append((p, s, q, bockstein, algorithm))
            return {}

    fake_oddp = types.ModuleType("oddp")
    fake_oddp.Steenrod = FakeSteenrod
    monkeypatch.setitem(sys.modules, "oddp", fake_oddp)

    assert one.operation(0, bockstein=True, algorithm="support").is_zero()
    assert calls == [(5, 0, 0, True, "support")]
    assert cohomology.operation_rank(0, 0, bockstein=True) == 0


def test_odd_primary_vanishing_target_degree_does_not_require_oddp(monkeypatch):
    cohomology = spaces.complex_projective_plane().cohomology(p=3)
    x = cohomology.basis(2)[0]
    monkeypatch.delitem(sys.modules, "oddp", raising=False)

    assert x.operation(1).is_zero()
    assert cohomology.operation_rank(2, 1) == 0


def test_basis_elements_have_python_style_squares():
    circle = SimplicialComplex([(0, 1), (0, 2), (1, 2)])
    cohomology = circle.cohomology()
    x = cohomology.basis(1)[0]

    assert x.sq(0) == x
    assert x.sq(1).is_zero()
    with pytest.raises(ValueError, match="nonnegative"):
        x.sq(-1)


def test_projective_plane_has_nonzero_sq1_rank():
    cohomology = spaces.real_projective_plane().cohomology()

    assert cohomology.betti_numbers() == {0: 1, 1: 1, 2: 1}
    assert cohomology.operation_rank(1, 1) == 1
    assert cohomology.basis(1)[0].sq(1) == cohomology.basis(2)[0]


def test_projective_3_space_betti_numbers_and_sq1():
    cohomology = spaces.real_projective_space(3).cohomology()

    assert cohomology.betti_numbers() == {0: 1, 1: 1, 2: 1, 3: 1}
    assert cohomology.operation_rank(1, 1) == 1
    assert cohomology.operation_rank(2, 1) == 0


def test_complex_projective_plane_betti_numbers_and_sq2():
    cohomology = spaces.complex_projective_plane().cohomology()

    assert cohomology.betti_numbers() == {0: 1, 2: 1, 4: 1}
    assert cohomology.operation_rank(2, 2) == 1
    assert cohomology.basis(2)[0].sq(2) == cohomology.basis(4)[0]


def test_complex_projective_space_catalog_includes_cp2_and_cp3():
    cp2 = spaces.complex_projective_space(2)
    cp3 = spaces.complex_projective_space(3)

    assert cp2 == spaces.complex_projective_plane()
    assert cp3.dimension == 6
    assert len(cp3.vertices) == 18
    assert len(cp3.facets) == 622


def test_moore_space_catalog_includes_mod_three_example():
    moore = spaces.moore_space(3)

    assert moore.dimension == 2
    assert len(moore.vertices) == 9
    assert len(moore.facets) == 19
    assert moore.cohomology(p=3).betti_numbers() == {0: 1, 1: 1, 2: 1}


def test_space_catalog_rejects_unknown_projective_spaces():
    with pytest.raises(NotImplementedError, match="dimensions 2 and 3"):
        spaces.real_projective_space(4)
    with pytest.raises(NotImplementedError, match="dimensions 2 and 3"):
        spaces.complex_projective_space(4)
    with pytest.raises(NotImplementedError, match="mod-3"):
        spaces.moore_space(5)
