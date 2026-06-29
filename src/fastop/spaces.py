"""Small catalog of finite simplicial complexes."""

from __future__ import annotations

from itertools import combinations

from fastop.simplicial import SimplicialComplex


def simplex(dimension: int) -> SimplicialComplex:
    """Return the standard ``dimension``-simplex."""
    _validate_dimension(dimension)
    return SimplicialComplex([tuple(range(dimension + 1))])


def sphere(dimension: int) -> SimplicialComplex:
    """Return the boundary of the standard ``(dimension + 1)``-simplex."""
    _validate_dimension(dimension)
    vertices = tuple(range(dimension + 2))
    return SimplicialComplex(combinations(vertices, dimension + 1))


def real_projective_plane() -> SimplicialComplex:
    """Return Sage's six-vertex triangulation of real projective plane."""
    return SimplicialComplex(_REAL_PROJECTIVE_PLANE_FACETS)


def real_projective_space(dimension: int) -> SimplicialComplex:
    """Return a catalog triangulation of real projective ``dimension``-space.

    This first catalog includes the small Sage triangulations of dimensions
    two and three.
    """
    if dimension == 2:
        return real_projective_plane()
    if dimension == 3:
        return SimplicialComplex(_REAL_PROJECTIVE_3_SPACE_FACETS)
    raise NotImplementedError("only real projective spaces of dimensions 2 and 3 are included")


def complex_projective_plane() -> SimplicialComplex:
    """Return Sage's nine-vertex triangulation of complex projective plane."""
    return SimplicialComplex(_COMPLEX_PROJECTIVE_PLANE_FACETS)


def _validate_dimension(dimension: int) -> None:
    if not isinstance(dimension, int) or isinstance(dimension, bool):
        raise TypeError("dimension must be an integer")
    if dimension < 0:
        raise ValueError("dimension must be nonnegative")


_REAL_PROJECTIVE_PLANE_FACETS = (
    (0, 1, 2),
    (0, 1, 5),
    (0, 2, 3),
    (0, 3, 4),
    (0, 4, 5),
    (1, 2, 4),
    (1, 3, 4),
    (1, 3, 5),
    (2, 3, 5),
    (2, 4, 5),
)

_REAL_PROJECTIVE_3_SPACE_FACETS = (
    (1, 2, 3, 7),
    (1, 2, 3, 11),
    (1, 2, 6, 9),
    (1, 2, 6, 11),
    (1, 2, 7, 9),
    (1, 3, 5, 10),
    (1, 3, 5, 11),
    (1, 3, 7, 10),
    (1, 4, 7, 9),
    (1, 4, 7, 10),
    (1, 4, 8, 9),
    (1, 4, 8, 10),
    (1, 5, 6, 8),
    (1, 5, 6, 11),
    (1, 5, 8, 10),
    (1, 6, 8, 9),
    (2, 3, 4, 8),
    (2, 3, 4, 11),
    (2, 3, 7, 8),
    (2, 4, 6, 10),
    (2, 4, 6, 11),
    (2, 4, 8, 10),
    (2, 5, 7, 8),
    (2, 5, 7, 9),
    (2, 5, 8, 10),
    (2, 5, 9, 10),
    (2, 6, 9, 10),
    (3, 4, 5, 9),
    (3, 4, 5, 11),
    (3, 4, 8, 9),
    (3, 5, 9, 10),
    (3, 6, 7, 8),
    (3, 6, 7, 10),
    (3, 6, 8, 9),
    (3, 6, 9, 10),
    (4, 5, 6, 7),
    (4, 5, 6, 11),
    (4, 5, 7, 9),
    (4, 6, 7, 10),
    (5, 6, 7, 8),
)

_COMPLEX_PROJECTIVE_PLANE_FACETS = (
    (1, 2, 3, 7, 8),
    (1, 2, 3, 7, 9),
    (1, 2, 3, 8, 9),
    (1, 2, 4, 5, 6),
    (1, 2, 4, 5, 9),
    (1, 2, 4, 6, 7),
    (1, 2, 4, 7, 9),
    (1, 2, 5, 6, 8),
    (1, 2, 5, 8, 9),
    (1, 2, 6, 7, 8),
    (1, 3, 4, 5, 6),
    (1, 3, 4, 5, 7),
    (1, 3, 4, 6, 8),
    (1, 3, 4, 7, 8),
    (1, 3, 5, 6, 9),
    (1, 3, 5, 7, 9),
    (1, 3, 6, 8, 9),
    (1, 4, 5, 7, 9),
    (1, 4, 6, 7, 8),
    (1, 5, 6, 8, 9),
    (2, 3, 4, 5, 6),
    (2, 3, 4, 5, 8),
    (2, 3, 4, 6, 9),
    (2, 3, 4, 8, 9),
    (2, 3, 5, 6, 7),
    (2, 3, 5, 7, 8),
    (2, 3, 6, 7, 9),
    (2, 4, 5, 8, 9),
    (2, 4, 6, 7, 9),
    (2, 5, 6, 7, 8),
    (3, 4, 5, 7, 8),
    (3, 4, 6, 8, 9),
    (3, 5, 6, 7, 9),
    (4, 5, 7, 8, 9),
    (4, 6, 7, 8, 9),
    (5, 6, 7, 8, 9),
)
