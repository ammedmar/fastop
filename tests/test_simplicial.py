import pytest

from fastop import SimplicialComplex, __version__


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
