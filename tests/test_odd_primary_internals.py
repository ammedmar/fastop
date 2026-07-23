import pytest
import sys
import types

from fastop import SimplicialComplex, spaces
import fastop.cohomology as cohomology_module
import fastop._cochain_evaluation as cochain_evaluation_module
from fastop._cochain_evaluation import (
    _auto_evaluation_algorithm,
    cochain_operation_vector,
    cochain_operation_vector_from_universal,
    evaluate_all_targets,
    evaluate_source_focused,
    evaluate_source_mod_3,
    evaluate_target_omissions,
)
from tests._oddp_oracle import (
    cochain_operation_vector_oddp,
)
from fastop._universal import (
    SignatureTable,
    UniversalOperation,
    native_universal_operation,
    universal_operation,
)


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


def test_cohomology_uses_public_odd_primary_degree_conventions():
    cohomology = SimplicialComplex.simplex_boundary(2).cohomology(p=5)

    assert cohomology._odd_primary_missing_vertices(2, bockstein=True) == 17
    assert cohomology._operation_target_degree(3, 2, bockstein=True) == 20


def test_cohomology_caches_universal_operations(monkeypatch):
    cohomology = SimplicialComplex.standard_simplex(2).cohomology(p=3)
    calls = []

    def fake_universal_operation(**kwargs):
        calls.append(kwargs)
        return UniversalOperation.from_terms(
            p=3,
            r=0,
            source_degree=1,
            bockstein=True,
            target_degree=2,
            missing_vertices_per_factor=1,
            terms={((0, 1), (1, 2), (0, 2)): 1},
        )

    monkeypatch.setattr(
        cohomology_module,
        "universal_operation",
        fake_universal_operation,
    )

    first = cohomology._universal_operation(
        operation_degree=0,
        source_degree=1,
        bockstein=True,
        target_degree=2,
        missing_vertices_per_factor=1,
        formula_source="auto",
    )
    second = cohomology._universal_operation(
        operation_degree=0,
        source_degree=1,
        bockstein=True,
        target_degree=2,
        missing_vertices_per_factor=1,
        formula_source="auto",
    )

    assert first is second
    assert len(calls) == 1


def test_cohomology_builds_boundary_data_only_for_requested_degree():
    cohomology = spaces.complex_projective_space(3).cohomology(p=5)

    assert cohomology._boundary_columns == {}
    assert cohomology._coboundary_columns == {}

    assert cohomology.betti_number(2) == 1
    assert set(cohomology._boundary_columns) == {2, 3}
    assert set(cohomology._coboundary_columns) == {1, 2}


def test_reference_oddp_cochain_bridge_returns_fastop_sparse_vector(monkeypatch):
    complex_ = spaces.complex_projective_space(2)
    target_faces = sorted(complex_.faces(4))
    target_face_to_index = {face: i for i, face in enumerate(target_faces)}
    target = target_faces[0]
    calls = []

    class FakeSteenrod:
        @staticmethod
        def cochain_operation(complex_by_degree, cochain, p, s, q, *, bockstein, algorithm):
            calls.append((complex_by_degree, cochain, p, s, q, bockstein, algorithm))
            return {target: 2}

    fake_oddp = types.ModuleType("oddp")
    fake_oddp.Steenrod = FakeSteenrod
    monkeypatch.setitem(sys.modules, "oddp", fake_oddp)

    vector = cochain_operation_vector_oddp(
        complex_,
        {(1,): 1},
        p=3,
        bockstein=False,
        oddp_s=-1,
        oddp_q=0,
        target_face_to_index=target_face_to_index,
        algorithm="direct",
    )

    assert vector == {target_face_to_index[target]: 2}
    complex_by_degree, cochain, p, s, q, bockstein, algorithm = calls[0]
    assert complex_by_degree[4] == set(complex_.faces(4))
    assert cochain == {(1,): 1}
    assert (p, s, q, bockstein, algorithm) == (3, -1, 0, False, "direct")


def test_computed_formula_uses_universal_data_and_native_evaluator():
    complex_ = spaces.complex_projective_space(2)
    target_faces = sorted(complex_.faces(2))
    target_face_to_index = {face: i for i, face in enumerate(target_faces)}
    target = target_faces[0]

    vector = cochain_operation_vector(
        complex_,
        {
            (target[0], target[2]): 1,
            (target[0], target[1]): 1,
            (target[1], target[2]): 1,
        },
        **_universal_kwargs(p=3, r=0, source_degree=1, bockstein=True),
        target_face_to_index=target_face_to_index,
        algorithm="source_focused",
        formula_source="computed",
    )

    assert vector == {target_face_to_index[target]: 2}


def test_universal_operation_reduces_coefficients():
    universal = UniversalOperation.from_terms(
        **_universal_kwargs(p=3, r=1, source_degree=2),
        terms={
            ((0, 1, 2), (2, 3, 4), (4, 5, 6)): 5,
            ((0,), (1,), (2,)): 3,
        },
    )

    assert universal.p == 3
    assert universal.r == 1
    assert universal.source_degree == 2
    assert universal.target_degree == 6
    assert universal.missing_vertices_per_factor == 4
    assert universal.terms == {((0, 1, 2), (2, 3, 4), (4, 5, 6)): 2}


def test_universal_operation_converts_to_omission_signature_table():
    universal = UniversalOperation.from_terms(
        **_universal_kwargs(p=3, r=1, source_degree=2),
        terms={
            ((0, 1, 2), (2, 3, 4), (4, 5, 6)): 2,
        },
    )

    assert universal.signature_table() == SignatureTable(
        p=3,
        r=1,
        source_degree=2,
        bockstein=False,
        target_degree=6,
        missing_vertices_per_factor=4,
        coefficients={
            (
                (3, 4, 5, 6),
                (0, 1, 5, 6),
                (0, 1, 2, 3),
            ): 2,
        },
    )


def test_signature_table_validates_homogeneous_tensor_factors():
    universal = UniversalOperation.from_terms(
        **_universal_kwargs(p=3, r=1, source_degree=2),
        terms={((0, 1), (2, 3, 4), (4, 5, 6)): 1},
    )

    with pytest.raises(ValueError, match="source degree"):
        universal.signature_table()


def test_computed_source_builds_native_universal_operation():
    universal = universal_operation(
        **_universal_kwargs(p=3, r=0, source_degree=1, bockstein=True),
        formula_source="computed",
    )

    assert universal.target_degree == 2
    assert universal.terms == {((0, 2), (0, 1), (1, 2)): 2}


def test_universal_operation_uses_low_dimensional_catalog_by_default(monkeypatch):
    class BrokenSteenrod:
        @staticmethod
        def chain_operations(p, s, q, *, bockstein):
            raise AssertionError("catalog lookup should not call oddp")

    fake_oddp = types.ModuleType("oddp")
    fake_oddp.Steenrod = BrokenSteenrod
    monkeypatch.setitem(sys.modules, "oddp", fake_oddp)

    universal = universal_operation(
        **_universal_kwargs(p=3, r=0, source_degree=1, bockstein=True)
    )

    assert universal.terms == {((0, 2), (0, 1), (1, 2)): 2}


def test_universal_operation_can_require_or_skip_catalog():
    operation = _universal_kwargs(p=3, r=0, source_degree=1, bockstein=True)

    assert universal_operation(**operation, formula_source="catalog").terms == {
        ((0, 2), (0, 1), (1, 2)): 2
    }
    with pytest.raises(NotImplementedError, match="catalog"):
        universal_operation(
            **_universal_kwargs(p=3, r=1, source_degree=5),
            formula_source="catalog",
        )


def test_native_universal_operation_builds_top_reduced_power_family():
    assert native_universal_operation(
        **_universal_kwargs(p=3, r=1, source_degree=2)
    ).terms == {((0, 1, 2), (2, 3, 4), (4, 5, 6)): 2}
    assert native_universal_operation(
        **_universal_kwargs(p=5, r=1, source_degree=2)
    ).terms == {
        ((0, 1, 2), (2, 3, 4), (4, 5, 6), (6, 7, 8), (8, 9, 10)): 4
    }
    assert native_universal_operation(
        **_universal_kwargs(p=3, r=2, source_degree=4)
    ).terms == {
        ((0, 1, 2, 3, 4), (4, 5, 6, 7, 8), (8, 9, 10, 11, 12)): 1
    }


def test_native_universal_operation_builds_bockstein_and_non_top_families():
    assert native_universal_operation(
        **_universal_kwargs(p=3, r=0, source_degree=1, bockstein=True)
    ).terms == {((0, 2), (0, 1), (1, 2)): 2}
    non_top = native_universal_operation(
        **_universal_kwargs(p=3, r=1, source_degree=3)
    )

    assert len(non_top.terms) == 19


def test_all_targets_evaluator_applies_tensor_terms():
    universal = UniversalOperation(
        p=3,
        r=1,
        source_degree=0,
        bockstein=False,
        target_degree=2,
        missing_vertices_per_factor=2,
        terms={((0,), (1,), (2,)): 2},
    )

    assert evaluate_all_targets(
        {(1, 2, 3), (1, 2, 4)},
        {(1,): 1, (2,): 2, (3,): 1, (4,): 0},
        universal,
    ) == {(1, 2, 3): 1}


def test_native_all_targets_matches_python_fallback(monkeypatch):
    if cochain_evaluation_module._native_evaluate_all_targets is None:
        pytest.skip("native extension is not built")

    universal = UniversalOperation(
        p=5,
        r=0,
        source_degree=2,
        bockstein=True,
        target_degree=3,
        missing_vertices_per_factor=1,
        terms={
            (
                (0, 1, 2),
                (0, 1, 3),
                (0, 2, 3),
                (1, 2, 3),
                (0, 1, 2),
            ): 4,
            (
                (0, 1, 3),
                (0, 1, 3),
                (0, 1, 2),
                (0, 2, 3),
                (0, 1, 3),
            ): 1,
        },
    )
    target_faces = {
        (0, 1, 2, 3),
        (0, 1, 2, 4),
        (0, 1, 3, 4),
    }
    cochain = {
        (0, 1, 2): 2,
        (0, 1, 3): 1,
        (0, 2, 3): 4,
        (1, 2, 3): 3,
        (0, 1, 4): 1,
        (0, 2, 4): 1,
        (1, 2, 4): 1,
    }

    native = evaluate_all_targets(target_faces, cochain, universal)
    monkeypatch.setattr(cochain_evaluation_module, "_native_evaluate_all_targets", None)

    assert native == evaluate_all_targets(target_faces, cochain, universal)


def test_source_focused_evaluator_applies_omission_signatures():
    universal = UniversalOperation(
        p=3,
        r=1,
        source_degree=0,
        bockstein=False,
        target_degree=2,
        missing_vertices_per_factor=2,
        terms={((0,), (1,), (2,)): 2},
    )

    assert evaluate_source_focused(
        {(1, 2, 3), (1, 2, 4)},
        {(1,): 1, (2,): 2, (3,): 1, (4,): 0},
        universal.signature_table(),
    ) == {(1, 2, 3): 1}


def test_source_focused_evaluator_matches_all_targets():
    universal = UniversalOperation(
        p=3,
        r=0,
        source_degree=1,
        bockstein=True,
        target_degree=2,
        missing_vertices_per_factor=1,
        terms={((0, 2), (0, 1), (1, 2)): 2},
    )
    target_faces = {(0, 1, 2), (0, 1, 3)}
    cochain = {(0, 1): 1, (0, 2): 2, (1, 2): 1, (0, 3): 1}

    assert evaluate_source_focused(
        target_faces,
        cochain,
        universal.signature_table(),
    ) == evaluate_all_targets(target_faces, cochain, universal)


def test_source_mod_3_evaluator_matches_all_targets():
    universal = UniversalOperation(
        p=3,
        r=0,
        source_degree=1,
        bockstein=True,
        target_degree=2,
        missing_vertices_per_factor=1,
        terms={((0, 2), (0, 1), (1, 2)): 2},
    )
    target_faces = {(0, 1, 2), (0, 1, 3)}
    cochain = {(0, 1): 1, (0, 2): 2, (1, 2): 1, (0, 3): 1}

    assert evaluate_source_mod_3(
        target_faces,
        cochain,
        universal.signature_table(),
    ) == evaluate_all_targets(target_faces, cochain, universal)


def test_native_source_mod_3_matches_python_fallback(monkeypatch):
    if cochain_evaluation_module._native_evaluate_source_mod_3_covered is None:
        pytest.skip("native extension is not built")

    universal = UniversalOperation(
        p=3,
        r=0,
        source_degree=1,
        bockstein=True,
        target_degree=2,
        missing_vertices_per_factor=1,
        terms={
            ((0, 2), (0, 1), (1, 2)): 2,
            ((0, 1), (0, 2), (1, 2)): 1,
        },
    )
    target_faces = {(0, 1, 2), (0, 1, 3)}
    cochain = {(0, 1): 1, (0, 2): 2, (1, 2): 1, (0, 3): 1}

    native = evaluate_source_mod_3(
        target_faces,
        cochain,
        universal.signature_table(),
    )
    monkeypatch.setattr(
        cochain_evaluation_module,
        "_native_evaluate_source_mod_3_covered",
        None,
    )

    assert native == evaluate_source_mod_3(
        target_faces,
        cochain,
        universal.signature_table(),
    )


def test_source_mod_3_evaluator_allows_repeated_source_factors():
    universal = UniversalOperation(
        p=3,
        r=0,
        source_degree=0,
        bockstein=True,
        target_degree=1,
        missing_vertices_per_factor=1,
        terms={((0,), (0,), (1,)): 1},
    )
    target_faces = {(1, 2)}
    cochain = {(1,): 2, (2,): 1}

    assert evaluate_source_mod_3(
        target_faces,
        cochain,
        universal.signature_table(),
    ) == evaluate_all_targets(target_faces, cochain, universal)


def test_source_mod_3_falls_back_for_uncovered_target_positions():
    universal = UniversalOperation(
        p=3,
        r=0,
        source_degree=0,
        bockstein=True,
        target_degree=3,
        missing_vertices_per_factor=3,
        terms={((0,), (0,), (0,)): 1},
    )
    target_faces = {(0, 1, 2, 3)}
    cochain = {(0,): 1}

    assert evaluate_source_mod_3(
        target_faces,
        cochain,
        universal.signature_table(),
    ) == evaluate_all_targets(target_faces, cochain, universal)


def test_target_omissions_evaluator_matches_all_targets():
    universal = UniversalOperation(
        p=3,
        r=0,
        source_degree=1,
        bockstein=True,
        target_degree=2,
        missing_vertices_per_factor=1,
        terms={((0, 2), (0, 1), (1, 2)): 2},
    )
    target_faces = {(0, 1, 2), (0, 1, 3)}
    cochain = {(0, 1): 1, (0, 2): 2, (1, 2): 1, (0, 3): 1}

    assert evaluate_target_omissions(
        target_faces,
        cochain,
        universal.signature_table(),
    ) == evaluate_all_targets(target_faces, cochain, universal)


def test_auto_evaluator_prefers_all_targets_for_dense_support():
    universal = UniversalOperation(
        p=3,
        r=1,
        source_degree=2,
        bockstein=False,
        target_degree=6,
        missing_vertices_per_factor=4,
        terms={((0, 1, 2), (2, 3, 4), (4, 5, 6)): 1},
    )
    target_faces = {tuple(range(7))}
    cochain = {
        face: 1
        for face in (
            (0, 1, 2),
            (0, 1, 3),
            (0, 1, 4),
            (0, 1, 5),
            (0, 1, 6),
        )
    }

    assert _auto_evaluation_algorithm(target_faces, cochain, universal) == "all_targets"


def test_auto_evaluator_keeps_all_targets_for_prime_three():
    universal = UniversalOperation(
        p=3,
        r=1,
        source_degree=2,
        bockstein=False,
        target_degree=6,
        missing_vertices_per_factor=4,
        terms={((0, 1, 2), (2, 3, 4), (4, 5, 6)): 1},
    )
    target_faces = {
        tuple(range(offset, offset + 7))
        for offset in range(20)
    }
    cochain = {(0, 1, 2): 1}

    assert _auto_evaluation_algorithm(target_faces, cochain, universal) == "all_targets"


def test_source_focused_evaluator_allows_repeated_source_factors():
    universal = UniversalOperation(
        p=3,
        r=0,
        source_degree=0,
        bockstein=True,
        target_degree=1,
        missing_vertices_per_factor=1,
        terms={((0,), (0,), (1,)): 1},
    )
    target_faces = {(1, 2)}
    cochain = {(1,): 2, (2,): 1}

    assert evaluate_source_focused(
        target_faces,
        cochain,
        universal.signature_table(),
    ) == evaluate_all_targets(target_faces, cochain, universal)
    assert evaluate_target_omissions(
        target_faces,
        cochain,
        universal.signature_table(),
    ) == evaluate_all_targets(target_faces, cochain, universal)
    assert evaluate_source_focused(
        target_faces,
        cochain,
        universal.signature_table(),
    ) == {(1, 2): 1}


def test_reference_vector_from_universal_uses_native_evaluator():
    complex_ = spaces.complex_projective_space(2)
    target_faces = sorted(complex_.faces(2))
    target_face_to_index = {face: i for i, face in enumerate(target_faces)}
    target = target_faces[0]
    universal = UniversalOperation(
        p=3,
        r=0,
        source_degree=1,
        bockstein=True,
        target_degree=2,
        missing_vertices_per_factor=1,
        terms={((0, 2), (0, 1), (1, 2)): 1},
    )
    cochain = {
        (target[0], target[2]): 1,
        (target[0], target[1]): 1,
        (target[1], target[2]): 1,
    }

    assert cochain_operation_vector_from_universal(
        complex_,
        cochain,
        universal,
        target_face_to_index,
    ) == {target_face_to_index[target]: 1}


def test_reference_vector_from_universal_can_use_target_omissions():
    complex_ = spaces.complex_projective_space(2)
    target_faces = sorted(complex_.faces(2))
    target_face_to_index = {face: i for i, face in enumerate(target_faces)}
    target = target_faces[0]
    universal = UniversalOperation(
        p=3,
        r=0,
        source_degree=1,
        bockstein=True,
        target_degree=2,
        missing_vertices_per_factor=1,
        terms={((0, 2), (0, 1), (1, 2)): 1},
    )
    cochain = {
        (target[0], target[2]): 1,
        (target[0], target[1]): 1,
        (target[1], target[2]): 1,
    }

    assert cochain_operation_vector_from_universal(
        complex_,
        cochain,
        universal,
        target_face_to_index,
        algorithm="target_omissions",
    ) == {target_face_to_index[target]: 1}


def test_source_mod_3_name_is_restricted_to_prime_three():
    complex_ = spaces.complex_projective_space(2)
    target_faces = sorted(complex_.faces(2))
    target_face_to_index = {face: i for i, face in enumerate(target_faces)}
    universal = UniversalOperation(
        p=5,
        r=0,
        source_degree=1,
        bockstein=True,
        target_degree=2,
        missing_vertices_per_factor=1,
        terms={((0, 2), (0, 2), (0, 2), (0, 1), (1, 2)): 1},
    )

    with pytest.raises(ValueError, match="p=3"):
        cochain_operation_vector_from_universal(
            complex_,
            {},
            universal,
            target_face_to_index,
            algorithm="source_mod_3",
        )


def test_reference_vector_from_universal_can_use_all_targets():
    complex_ = spaces.complex_projective_space(2)
    target_faces = sorted(complex_.faces(4))
    target_face_to_index = {face: i for i, face in enumerate(target_faces)}
    target = target_faces[0]
    universal = UniversalOperation(
        p=3,
        r=1,
        source_degree=0,
        bockstein=False,
        target_degree=4,
        missing_vertices_per_factor=4,
        terms={tuple((i,) for i in range(5)): 1},
    )
    cochain = {vertex: 1 for vertex in zip(target)}

    assert cochain_operation_vector_from_universal(
        complex_,
        cochain,
        universal,
        target_face_to_index,
        algorithm="all_targets",
    ) == {target_face_to_index[target]: 1}
