"""Cohomology over the field with two elements."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import TYPE_CHECKING

from fastop._f2 import CoordinateBasis, iter_bits, nullspace, rank, rref

if TYPE_CHECKING:
    from fastop.simplicial import Simplex, SimplicialComplex


@dataclass(frozen=True)
class _DegreeData:
    faces: tuple["Simplex", ...]
    face_to_index: dict["Simplex", int]
    cocycle_basis: tuple[int, ...]
    projector: CoordinateBasis


class Mod2Cohomology:
    """Mod-2 cohomology of a finite simplicial complex."""

    def __init__(self, complex_: "SimplicialComplex", *, reduced: bool = False):
        self.complex = complex_
        self.p = 2
        self.reduced = reduced
        self._faces = {
            degree: tuple(sorted(complex_.faces(degree)))
            for degree in range(complex_.dimension + 1)
        }
        self._face_to_index = {
            degree: {face: index for index, face in enumerate(faces)}
            for degree, faces in self._faces.items()
        }
        self._boundary_columns = self._build_boundary_columns()
        self._coboundary_columns = self._build_coboundary_columns()
        self._degree_data = self._build_degree_data()

    @property
    def dimension(self) -> int:
        """Return the dimension of the underlying complex."""
        return self.complex.dimension

    def basis(self, degree: int | None = None):
        """Return cohomology basis elements, optionally in one degree."""
        if degree is None:
            return tuple(
                self.element({n: 1 << i})
                for n in range(self.dimension + 1)
                for i in range(self.betti_number(n))
            )
        return tuple(self.element({degree: 1 << i}) for i in range(self.betti_number(degree)))

    def betti_number(self, degree: int) -> int:
        """Return the mod-2 Betti number in ``degree``."""
        data = self._degree_data.get(degree)
        if data is None:
            return 0
        return len(data.cocycle_basis)

    def betti_numbers(self) -> dict[int, int]:
        """Return all nonzero mod-2 Betti numbers."""
        return {
            degree: betti
            for degree in range(self.dimension + 1)
            if (betti := self.betti_number(degree))
        }

    def zero(self) -> "Mod2CohomologyElement":
        """Return the zero cohomology element."""
        return self.element({})

    def element(self, coordinates: dict[int, int]) -> "Mod2CohomologyElement":
        """Create a cohomology element from degree-index bit coordinates."""
        clean = {}
        for degree, vector in coordinates.items():
            mask = (1 << self.betti_number(degree)) - 1
            vector &= mask
            if vector:
                clean[degree] = vector
        return Mod2CohomologyElement(self, clean)

    def cocycle_vector(self, element: "Mod2CohomologyElement", degree: int) -> int:
        """Return the representative cocycle vector for one homogeneous part."""
        vector = 0
        data = self._degree_data[degree]
        coordinates = element._coordinates.get(degree, 0)
        for index in iter_bits(coordinates):
            vector ^= data.cocycle_basis[index]
        return vector

    def cocycle(self, element: "Mod2CohomologyElement", degree: int):
        """Return the representative cocycle as a simplex-to-coefficient dict."""
        data = self._degree_data[degree]
        vector = self.cocycle_vector(element, degree)
        return {data.faces[index]: 1 for index in iter_bits(vector)}

    def project_cocycle(self, degree: int, vector: int) -> "Mod2CohomologyElement":
        """Project a cocycle vector to cohomology coordinates."""
        data = self._degree_data.get(degree)
        if data is None:
            if vector:
                raise ValueError("cannot project a nonzero vector outside the complex")
            return self.zero()
        return self.element({degree: data.projector.coordinates(vector)})

    def square(self, element: "Mod2CohomologyElement", k: int) -> "Mod2CohomologyElement":
        """Apply the Steenrod square ``Sq^k``."""
        if not isinstance(k, int) or isinstance(k, bool):
            raise TypeError("k must be an integer")
        if k < 0:
            raise ValueError("k must be nonnegative")
        answer = self.zero()
        for degree in sorted(element._coordinates):
            answer += self._square_homogeneous(element, degree, k)
        return answer

    def operation_matrix(self, degree: int, k: int) -> list[int]:
        """Return the columns of ``Sq^k: H^degree -> H^(degree+k)``."""
        columns = []
        for basis_element in self.basis(degree):
            columns.append(basis_element.sq(k)._coordinates.get(degree + k, 0))
        return columns

    def operation_rank(self, degree: int, k: int) -> int:
        """Return the rank of ``Sq^k: H^degree -> H^(degree+k)``."""
        return rank(self.operation_matrix(degree, k))

    def _square_homogeneous(
        self, element: "Mod2CohomologyElement", degree: int, k: int
    ) -> "Mod2CohomologyElement":
        if k == 0:
            return self.element({degree: element._coordinates.get(degree, 0)})
        if k > degree:
            return self.zero()

        target_degree = degree + k
        target_data = self._degree_data.get(target_degree)
        if target_data is None:
            return self.zero()

        source_data = self._degree_data[degree]
        cocycle_vector = self.cocycle_vector(element, degree)
        support = [source_data.faces[index] for index in iter_bits(cocycle_vector)]
        target_support = _steenrod_square_support(
            target_degree + 1,
            support,
            set(target_data.faces),
        )
        target_vector = 0
        for simplex in target_support:
            target_vector ^= 1 << target_data.face_to_index[simplex]
        return self.project_cocycle(target_degree, target_vector)

    def _build_boundary_columns(self) -> dict[int, list[int]]:
        columns: dict[int, list[int]] = {0: [0 for _ in self._faces.get(0, ())]}
        for degree in range(1, self.dimension + 1):
            lower_index = self._face_to_index[degree - 1]
            degree_columns = []
            for simplex in self._faces[degree]:
                column = 0
                for face in _codimension_one_faces(simplex):
                    column ^= 1 << lower_index[face]
                degree_columns.append(column)
            columns[degree] = degree_columns
        return columns

    def _build_coboundary_columns(self) -> dict[int, list[int]]:
        columns = {}
        for degree in range(self.dimension + 1):
            domain_faces = self._faces.get(degree, ())
            target_faces = self._faces.get(degree + 1, ())
            degree_columns = []
            for domain_index in range(len(domain_faces)):
                column = 0
                for target_index, boundary in enumerate(
                    self._boundary_columns.get(degree + 1, ())
                ):
                    if (boundary >> domain_index) & 1:
                        column |= 1 << target_index
                degree_columns.append(column)
            columns[degree] = degree_columns
            if not target_faces:
                columns[degree] = [0 for _ in domain_faces]
        return columns

    def _build_degree_data(self) -> dict[int, _DegreeData]:
        all_data = {}
        for degree in range(self.dimension + 1):
            faces = self._faces[degree]
            coboundary = self._coboundary_columns[degree]
            cycles = nullspace(coboundary, len(self._faces.get(degree + 1, ())))
            boundary_vectors = list(rref(self._coboundary_columns.get(degree - 1, [])).values())
            if self.reduced and degree == 0 and faces:
                boundary_vectors.append((1 << len(faces)) - 1)

            projector = CoordinateBasis()
            for boundary in boundary_vectors:
                projector.add(boundary, 0)

            cocycle_basis = []
            for cycle in sorted(cycles):
                basis_coordinate = 1 << len(cocycle_basis)
                if projector.add(cycle, basis_coordinate):
                    cocycle_basis.append(cycle)

            all_data[degree] = _DegreeData(
                faces=faces,
                face_to_index=self._face_to_index[degree],
                cocycle_basis=tuple(cocycle_basis),
                projector=projector,
            )
        return all_data


class Mod2CohomologyElement:
    """An element of a mod-2 cohomology vector space."""

    def __init__(self, parent: Mod2Cohomology, coordinates: dict[int, int]):
        self.parent = parent
        self._coordinates = dict(coordinates)

    def sq(self, k: int) -> "Mod2CohomologyElement":
        """Return ``Sq^k`` applied to this cohomology element."""
        return self.parent.square(self, k)

    def cocycle(self, degree: int | None = None):
        """Return the chosen representative cocycle of a homogeneous element."""
        if degree is None:
            degree = self.degree()
        return self.parent.cocycle(self, degree)

    def degree(self) -> int:
        """Return the degree of a nonzero homogeneous element."""
        if len(self._coordinates) != 1:
            raise ValueError("element is not homogeneous")
        return next(iter(self._coordinates))

    def is_zero(self) -> bool:
        """Return whether this is the zero element."""
        return not self._coordinates

    def __add__(self, other: "Mod2CohomologyElement") -> "Mod2CohomologyElement":
        if self.parent is not other.parent:
            return NotImplemented
        coordinates = dict(self._coordinates)
        for degree, vector in other._coordinates.items():
            coordinates[degree] = coordinates.get(degree, 0) ^ vector
            if not coordinates[degree]:
                del coordinates[degree]
        return self.parent.element(coordinates)

    __sub__ = __add__

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, Mod2CohomologyElement)
            and self.parent is other.parent
            and self._coordinates == other._coordinates
        )

    def __bool__(self) -> bool:
        return bool(self._coordinates)

    def __repr__(self) -> str:
        if not self._coordinates:
            return "0"
        terms = []
        for degree in sorted(self._coordinates):
            for index in iter_bits(self._coordinates[degree]):
                terms.append(f"h^{degree},{index}")
        return " + ".join(terms)


def _codimension_one_faces(simplex: "Simplex") -> tuple["Simplex", ...]:
    return tuple(simplex[:index] + simplex[index + 1 :] for index in range(len(simplex)))


def _steenrod_square_support(
    target_length: int,
    cocycle_support: list["Simplex"],
    target_simplices: set["Simplex"],
) -> set["Simplex"]:
    answer: set["Simplex"] = set()
    for left, right in combinations(cocycle_support, 2):
        left_vertices = set(left)
        right_vertices = set(right)
        union = left_vertices | right_vertices
        if len(union) != target_length:
            continue

        simplex = tuple(sorted(union))
        if simplex not in target_simplices:
            continue

        left_only = left_vertices - right_vertices
        right_only = right_vertices - left_vertices
        symmetric_difference = sorted(left_only | right_only)
        indices = {}
        for vertex in symmetric_difference:
            indices[vertex] = (
                simplex.index(vertex) + symmetric_difference.index(vertex)
            ) % 2
        left_indices = {indices[vertex] for vertex in left_only}
        right_indices = {indices[vertex] for vertex in right_only}
        if (left_indices == {0} and right_indices == {1}) or (
            left_indices == {1} and right_indices == {0}
        ):
            if simplex in answer:
                answer.remove(simplex)
            else:
                answer.add(simplex)
    return answer
