"""Finite Delta-complexes represented by face maps."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, TYPE_CHECKING

from fastop._group_action import (
    CellAction,
    cell_orbits,
    normalize_actions,
    validate_free_action,
    validate_permutations,
)
from fastop.cohomology import PrimeFieldCohomology

if TYPE_CHECKING:
    from fastop.simplicial import SimplicialComplex

Face = tuple[int, ...]
FaceTable = tuple[tuple[Face, ...], ...]


@dataclass(frozen=True)
class DeltaComplex:
    """A finite semi-simplicial set with cells indexed in each degree.

    The entry ``face_maps[d][sigma][i]`` is the index of the face
    :math:`d_i(sigma)` among the ``(d-1)``-cells.  This is the same dense
    representation returned by Sage's ``DeltaComplex.cells()`` method after
    its degree ``-1`` entry is removed.
    """

    face_maps: FaceTable

    def __init__(
        self,
        face_maps: Iterable[Iterable[Iterable[int]]],
        *,
        check: bool = True,
    ):
        normalized = tuple(
            tuple(tuple(face) for face in degree_faces)
            for degree_faces in face_maps
        )
        while normalized and not normalized[-1]:
            normalized = normalized[:-1]
        if not normalized or not normalized[0]:
            raise ValueError("a Delta-complex needs at least one vertex")
        self._validate_shape(normalized)
        if check:
            self._validate_identities(normalized)
        object.__setattr__(self, "face_maps", normalized)

    @classmethod
    def from_face_maps(
        cls,
        face_maps: Iterable[Iterable[Iterable[int]]],
        *,
        check: bool = True,
    ) -> "DeltaComplex":
        """Create a Delta-complex from dense face-index tables."""
        return cls(face_maps, check=check)

    @classmethod
    def from_simplicial_complex(cls, complex_: "SimplicialComplex") -> "DeltaComplex":
        """Forget global vertices and retain the induced face maps."""
        cells = {
            degree: tuple(sorted(complex_.faces(degree)))
            for degree in range(complex_.dimension + 1)
        }
        indices = {
            degree: {cell: index for index, cell in enumerate(degree_cells)}
            for degree, degree_cells in cells.items()
        }
        face_maps = [tuple(() for _ in cells[0])]
        for degree in range(1, complex_.dimension + 1):
            face_maps.append(tuple(
                tuple(
                    indices[degree - 1][cell[:i] + cell[i + 1 :]]
                    for i in range(degree + 1)
                )
                for cell in cells[degree]
            ))
        return cls(face_maps, check=False)

    @classmethod
    def from_sage(cls, complex_) -> "DeltaComplex":
        """Copy a finite Sage ``DeltaComplex`` without retaining Sage objects."""
        cells = {int(degree): values for degree, values in complex_.cells().items()}
        nonnegative_degrees = [degree for degree in cells if degree >= 0]
        if not nonnegative_degrees:
            raise ValueError("the Sage Delta-complex has no vertices")
        dimension = max(nonnegative_degrees)
        return cls(
            tuple(
                tuple(tuple(int(face) for face in cell) for cell in cells.get(degree, ()))
                for degree in range(dimension + 1)
            )
        )

    @property
    def dimension(self) -> int:
        """Return the largest cell dimension."""
        return len(self.face_maps) - 1

    @property
    def supports_vertex_algorithms(self) -> bool:
        """Return whether cells carry globally comparable vertex sets."""
        return False

    def cells(self, dimension: int | None = None):
        """Return cell indices, either in one dimension or grouped by degree."""
        if dimension is None:
            return {
                degree: range(len(degree_faces))
                for degree, degree_faces in enumerate(self.face_maps)
            }
        if dimension < 0 or dimension > self.dimension:
            return ()
        return range(len(self.face_maps[dimension]))

    def f_vector(self) -> tuple[int, ...]:
        """Return the number of cells in every nonnegative dimension."""
        return tuple(len(degree_faces) for degree_faces in self.face_maps)

    def face(self, degree: int, cell: int, index: int) -> int:
        """Return the ``index``-th face of one cell."""
        return self.face_maps[degree][cell][index]

    def restrict(self, degree: int, cell: int, positions: tuple[int, ...]) -> int:
        """Restrict a cell to the listed local vertex positions."""
        if not positions:
            raise ValueError("a restriction needs at least one retained position")
        if tuple(sorted(set(positions))) != positions:
            raise ValueError("retained positions must be strictly increasing")
        if positions[-1] > degree:
            raise ValueError("retained position is outside the cell")

        omitted = set(range(degree + 1)).difference(positions)
        current_degree = degree
        for index in sorted(omitted, reverse=True):
            cell = self.face(current_degree, cell, index)
            current_degree -= 1
        return cell

    def cohomology(
        self,
        p: int = 2,
        *,
        reduced: bool = False,
        convention: int = 1,
    ) -> PrimeFieldCohomology:
        """Return mod-``p`` cohomology with a chosen basis."""
        return PrimeFieldCohomology(self, p=p, reduced=reduced, convention=convention)

    def quotient(
        self,
        generators: Iterable[Iterable[Iterable[int]]],
        *,
        require_free: bool = False,
    ) -> "DeltaComplex":
        """Return the quotient by strict semi-simplicial automorphisms.

        Each generator is a permutation of the cells in every degree.  It
        must commute with every face map.  When ``require_free`` is true, the
        generated finite group is enumerated and every nonidentity element is
        checked for fixed cells.
        """
        actions = normalize_actions(generators)
        if not actions:
            return self
        counts = self.f_vector()
        validate_permutations(actions, counts)
        for action in actions:
            self._validate_action(action)

        if require_free:
            validate_free_action(actions, counts)

        orbit_indices, orbit_representatives = cell_orbits(actions, counts)

        quotient_faces = [
            tuple(() for _ in orbit_representatives[0])
        ]
        for degree in range(1, self.dimension + 1):
            quotient_faces.append(tuple(
                tuple(
                    orbit_indices[degree - 1][self.face(degree, representative, i)]
                    for i in range(degree + 1)
                )
                for representative in orbit_representatives[degree]
            ))
        return DeltaComplex(quotient_faces)

    def _validate_action(self, action: CellAction) -> None:
        for degree, permutation in enumerate(action):
            cell_count = len(self.face_maps[degree])
            if degree == 0:
                continue
            for cell in range(cell_count):
                for i in range(degree + 1):
                    acted_face = action[degree - 1][self.face(degree, cell, i)]
                    face_of_acted = self.face(degree, permutation[cell], i)
                    if acted_face != face_of_acted:
                        raise ValueError("cell action does not commute with face maps")

    @staticmethod
    def _validate_shape(face_maps: FaceTable) -> None:
        for degree, degree_faces in enumerate(face_maps):
            if degree and not degree_faces:
                raise ValueError("cell dimensions cannot have gaps")
            for face in degree_faces:
                expected_faces = 0 if degree == 0 else degree + 1
                if len(face) != expected_faces:
                    raise ValueError(
                        f"each {degree}-cell needs exactly {expected_faces} faces"
                    )
                if any(not isinstance(value, int) or isinstance(value, bool) for value in face):
                    raise TypeError("face indices must be integers")
                if degree and any(
                    value < 0 or value >= len(face_maps[degree - 1])
                    for value in face
                ):
                    raise ValueError("face index is outside the preceding dimension")

    @staticmethod
    def _validate_identities(face_maps: FaceTable) -> None:
        for degree in range(2, len(face_maps)):
            lower_faces = face_maps[degree - 1]
            for cell, faces in enumerate(face_maps[degree]):
                for i in range(degree):
                    for j in range(i + 1, degree + 1):
                        left = lower_faces[faces[j]][i]
                        right = lower_faces[faces[i]][j - 1]
                        if left != right:
                            raise ValueError(
                                "face maps violate the semi-simplicial identities "
                                f"at degree {degree}, cell {cell}, i={i}, j={j}"
                            )
