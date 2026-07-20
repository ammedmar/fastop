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


def test_cohomology_accepts_only_supported_conventions():
    complex_ = SimplicialComplex([(0, 1)])

    assert complex_.cohomology(convention=1).convention == 1
    assert complex_.cohomology(convention=-1).convention == -1
    with pytest.raises(ValueError, match="convention"):
        complex_.cohomology(convention=0)


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


def test_odd_primary_operation_uses_universal_data_and_projects(monkeypatch):
    cohomology = spaces.moore_space(3).cohomology(p=3)
    x = cohomology.basis(1)[0]
    calls = []

    class FakeSteenrod:
        @staticmethod
        def chain_operations(p, s, q, *, bockstein):
            calls.append((p, s, q, bockstein))
            return {((0, 2), (0, 1), (1, 2)): 2}

    fake_oddp = types.ModuleType("oddp")
    fake_oddp.Steenrod = FakeSteenrod
    monkeypatch.setitem(sys.modules, "oddp", fake_oddp)

    assert x.operation(
        0,
        bockstein=True,
        algorithm="direct",
        formula_source="computed",
    ) == cohomology.basis(2)[0]
    assert calls[0] == (3, 0, -1, True)
    assert cohomology.operation_rank(1, 0, bockstein=True, algorithm="direct") == 1


def test_top_reduced_power_does_not_require_oddp(monkeypatch):
    cohomology = spaces.complex_projective_space(3).cohomology(p=3)
    x = cohomology.basis(2)[0]
    fake_oddp = types.ModuleType("oddp")
    monkeypatch.setitem(sys.modules, "oddp", fake_oddp)

    assert x.operation(1) == cohomology.basis(6)[0]


def test_odd_primary_default_auto_matches_all_targets():
    cohomology = spaces.complex_projective_space(3).cohomology(p=3)
    x = cohomology.basis(2)[0]

    assert x.operation(1) == x.operation(1, algorithm="all_targets")


def test_odd_primary_operation_uses_oddp_universal_conventions(monkeypatch):
    cohomology = spaces.moore_space(3).cohomology(p=3)
    x = cohomology.basis(1)[0]
    calls = []

    class FakeSteenrod:
        @staticmethod
        def chain_operations(p, s, q, *, bockstein):
            calls.append((p, s, q, bockstein))
            return {}

    fake_oddp = types.ModuleType("oddp")
    fake_oddp.Steenrod = FakeSteenrod
    monkeypatch.setitem(sys.modules, "oddp", fake_oddp)

    assert x.operation(
        0,
        bockstein=True,
        algorithm="source_focused",
        formula_source="computed",
    ).is_zero()
    assert calls == [(3, 0, -1, True)]


def test_odd_primary_vanishing_target_degree_does_not_require_oddp(monkeypatch):
    cohomology = spaces.complex_projective_plane().cohomology(p=3)
    x = cohomology.basis(2)[0]
    monkeypatch.delitem(sys.modules, "oddp", raising=False)

    assert x.operation(1).is_zero()
    assert cohomology.operation_rank(2, 1) == 0


def test_operation_instability_bounds_are_enforced_before_evaluation():
    mod_two = spaces.real_projective_plane().cohomology(p=2)
    odd = spaces.complex_projective_plane().cohomology(p=3)

    assert mod_two._satisfies_instability_bound(1, 1, bockstein=False)
    assert not mod_two._satisfies_instability_bound(1, 2, bockstein=False)
    assert odd._satisfies_instability_bound(2, 1, bockstein=False)
    assert not odd._satisfies_instability_bound(1, 1, bockstein=True)


def test_basis_elements_have_python_style_squares():
    circle = SimplicialComplex([(0, 1), (0, 2), (1, 2)])
    cohomology = circle.cohomology()
    x = cohomology.basis(1)[0]

    assert x.sq(0) == x
    assert x.sq(1).is_zero()
    assert x.sq(-1).is_zero()


def test_homological_convention_uses_negative_operation_indices():
    cohomology = spaces.real_projective_plane().cohomology(convention=-1)
    x = cohomology.basis(1)[0]

    assert x.sq(-1) == cohomology.basis(2)[0]
    assert x.sq(1).is_zero()


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
