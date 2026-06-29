import sys
import types

from fastop import spaces
from fastop._odd_primary.indices import OperationIndex
from fastop._odd_primary.reference import cochain_operation_vector


def test_operation_index_translates_fastop_to_oddp_conventions():
    index = OperationIndex(p=5, r=2, source_degree=3, bockstein=True)

    assert index.oddp_s == -2
    assert index.oddp_q == -3
    assert index.missing_vertices_per_factor == 17
    assert index.target_degree == 20


def test_reference_bridge_returns_fastop_sparse_vector(monkeypatch):
    complex_ = spaces.complex_projective_plane()
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

    vector = cochain_operation_vector(
        complex_,
        {(1,): 1},
        OperationIndex(p=3, r=1, source_degree=0),
        target_face_to_index,
        algorithm="direct",
    )

    assert vector == {target_face_to_index[target]: 2}
    complex_by_degree, cochain, p, s, q, bockstein, algorithm = calls[0]
    assert complex_by_degree[4] == set(complex_.faces(4))
    assert cochain == {(1,): 1}
    assert (p, s, q, bockstein, algorithm) == (3, -1, 0, False, "direct")
