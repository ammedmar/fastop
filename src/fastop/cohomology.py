"""Cohomology over finite prime fields."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from fastop._cochain_evaluation import (
    cochain_operation_vector_from_universal,
    evaluate_source_mod_2,
)
from fastop._linear_algebra import (
    CoordinateBasis,
    Vector,
    clean_vector,
    column_image_and_kernel_basis,
    is_prime,
    rank,
    vector_add,
    vector_scale,
)
from fastop._universal import universal_operation

if TYPE_CHECKING:
    from fastop.simplicial import Simplex, SimplicialComplex


@dataclass(frozen=True)
class _DegreeData:
    faces: tuple["Simplex", ...]
    face_to_index: dict["Simplex", int]
    cocycle_basis: tuple[Vector, ...]
    projector: CoordinateBasis


class PrimeFieldCohomology:
    """Cohomology of a finite simplicial complex over ``F_p``."""

    def __init__(
        self,
        complex_: "SimplicialComplex",
        p: int = 2,
        *,
        reduced: bool = False,
        convention: int = 1,
    ):
        if not is_prime(p):
            raise ValueError("p must be a prime")
        if convention not in {1, -1}:
            raise ValueError("convention must be 1 or -1")
        self.complex = complex_
        self.p = p
        self.reduced = reduced
        self.convention = convention
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
        self._universal_operations = {}

    @property
    def dimension(self) -> int:
        """Return the dimension of the underlying complex."""
        return self.complex.dimension

    def basis(self, degree: int | None = None):
        """Return cohomology basis elements, optionally in one degree."""
        if degree is None:
            return tuple(
                self.element({n: {i: 1}})
                for n in range(self.dimension + 1)
                for i in range(self.betti_number(n))
            )
        return tuple(self.element({degree: {i: 1}}) for i in range(self.betti_number(degree)))

    def betti_number(self, degree: int) -> int:
        """Return the mod-``p`` Betti number in ``degree``."""
        data = self._degree_data.get(degree)
        if data is None:
            return 0
        return len(data.cocycle_basis)

    def betti_numbers(self) -> dict[int, int]:
        """Return all nonzero mod-``p`` Betti numbers."""
        return {
            degree: betti
            for degree in range(self.dimension + 1)
            if (betti := self.betti_number(degree))
        }

    def __repr__(self) -> str:
        return (
            f"PrimeFieldCohomology(p={self.p}, "
            f"dimension={self.dimension}, "
            f"convention={self.convention}, "
            f"betti_numbers={self.betti_numbers()})"
        )

    def zero(self) -> "PrimeFieldCohomologyElement":
        """Return the zero cohomology element."""
        return self.element({})

    def element(self, coordinates: dict[int, Vector]) -> "PrimeFieldCohomologyElement":
        """Create a cohomology element from degree-index coordinates."""
        clean = {}
        for degree, vector in coordinates.items():
            dimension = self.betti_number(degree)
            reduced = {
                index: coefficient % self.p
                for index, coefficient in vector.items()
                if 0 <= index < dimension and coefficient % self.p
            }
            if reduced:
                clean[degree] = reduced
        return PrimeFieldCohomologyElement(self, clean)

    def cocycle_vector(
        self, element: "PrimeFieldCohomologyElement", degree: int
    ) -> Vector:
        """Return the representative cocycle vector for one homogeneous part."""
        vector: Vector = {}
        data = self._degree_data[degree]
        coordinates = element._coordinates.get(degree, {})
        for index, coefficient in coordinates.items():
            vector = vector_add(vector, data.cocycle_basis[index], self.p, coefficient)
        return vector

    def cocycle(self, element: "PrimeFieldCohomologyElement", degree: int):
        """Return the representative cocycle as a simplex-to-coefficient dict."""
        data = self._degree_data[degree]
        vector = self.cocycle_vector(element, degree)
        return {data.faces[index]: coefficient for index, coefficient in sorted(vector.items())}

    def project_cocycle(self, degree: int, vector: Vector) -> "PrimeFieldCohomologyElement":
        """Project a cocycle vector to cohomology coordinates."""
        vector = clean_vector(vector, self.p)
        data = self._degree_data.get(degree)
        if data is None:
            if vector:
                raise ValueError("cannot project a nonzero vector outside the complex")
            return self.zero()
        return self.element({degree: data.projector.coordinates(vector)})

    def square(self, element: "PrimeFieldCohomologyElement", k: int) -> "PrimeFieldCohomologyElement":
        """Apply the Steenrod square ``Sq^k``.

        Steenrod squares are only available over ``F_2``.
        """
        if self.p != 2:
            raise NotImplementedError("Steenrod squares are only implemented for p=2")
        if not isinstance(k, int) or isinstance(k, bool):
            raise TypeError("k must be an integer")
        answer = self.zero()
        for degree in sorted(element._coordinates):
            answer += self._square_homogeneous(element, degree, k)
        return answer

    def operation(
        self,
        element: "PrimeFieldCohomologyElement",
        r: int,
        *,
        bockstein: bool = False,
        algorithm: str = "auto",
    ) -> "PrimeFieldCohomologyElement":
        """Apply the Steenrod operation selected by the coefficient prime."""
        self._validate_operation_index(r)
        if not isinstance(bockstein, bool):
            raise TypeError("bockstein must be a boolean")
        if self.p == 2:
            if bockstein:
                raise NotImplementedError("Bockstein operations require an odd prime")
            return self.square(element, r)

        answer = self.zero()
        for degree in sorted(element._coordinates):
            answer += self._odd_primary_operation_homogeneous(
                element,
                degree,
                r,
                bockstein=bockstein,
                algorithm=algorithm,
            )
        return answer

    def operation_matrix(
        self,
        degree: int,
        r: int,
        *,
        bockstein: bool = False,
        algorithm: str = "auto",
    ) -> list[Vector]:
        """Return columns of the selected Steenrod operation from ``H^degree``."""
        self._validate_operation_index(r)
        if not isinstance(bockstein, bool):
            raise TypeError("bockstein must be a boolean")
        target_degree = self._operation_target_degree(degree, r, bockstein=bockstein)
        return [
            basis_element.operation(
                r,
                bockstein=bockstein,
                algorithm=algorithm,
            )._coordinates.get(target_degree, {})
            for basis_element in self.basis(degree)
        ]

    def operation_rank(
        self,
        degree: int,
        r: int,
        *,
        bockstein: bool = False,
        algorithm: str = "auto",
    ) -> int:
        """Return the rank of the selected Steenrod operation from ``H^degree``."""
        return rank(
            self.operation_matrix(
                degree,
                r,
                bockstein=bockstein,
                algorithm=algorithm,
            ),
            self.p,
        )

    def _square_homogeneous(
        self, element: "PrimeFieldCohomologyElement", degree: int, k: int
    ) -> "PrimeFieldCohomologyElement":
        k = self.convention * k
        if k < 0:
            return self.zero()
        if k == 0:
            return self.element({degree: element._coordinates.get(degree, {})})
        if k > degree:
            return self.zero()

        target_degree = degree + k
        target_data = self._degree_data.get(target_degree)
        if target_data is None:
            return self.zero()

        source_data = self._degree_data[degree]
        cocycle_vector = self.cocycle_vector(element, degree)
        support = [source_data.faces[index] for index in cocycle_vector]
        target_support = evaluate_source_mod_2(
            target_degree + 1,
            support,
            set(target_data.faces),
        )
        target_vector = {
            target_data.face_to_index[simplex]: 1
            for simplex in target_support
        }
        return self.project_cocycle(target_degree, target_vector)

    def _validate_operation_index(self, r: int) -> None:
        if not isinstance(r, int) or isinstance(r, bool):
            raise TypeError("r must be an integer")

    def _operation_target_degree(self, degree: int, r: int, *, bockstein: bool) -> int:
        if self.p == 2:
            if bockstein:
                raise NotImplementedError("Bockstein operations require an odd prime")
            return degree + self.convention * r
        operation_degree = self.convention * r
        return degree + self._odd_primary_missing_vertices(
            operation_degree,
            bockstein=bockstein,
        )

    def _odd_primary_missing_vertices(self, r: int, *, bockstein: bool) -> int:
        return 2 * r * (self.p - 1) + int(bockstein)

    def _oddp_operation_index(self, r: int) -> int:
        return -r

    def _oddp_source_degree(self, degree: int) -> int:
        return -degree

    def _universal_operation(
        self,
        *,
        operation_degree: int,
        source_degree: int,
        bockstein: bool,
        target_degree: int,
        missing_vertices_per_factor: int,
    ):
        key = (
            operation_degree,
            source_degree,
            bockstein,
            target_degree,
            missing_vertices_per_factor,
        )
        cached = self._universal_operations.get(key)
        if cached is None:
            cached = universal_operation(
                p=self.p,
                r=operation_degree,
                source_degree=source_degree,
                bockstein=bockstein,
                target_degree=target_degree,
                missing_vertices_per_factor=missing_vertices_per_factor,
                oddp_s=self._oddp_operation_index(operation_degree),
                oddp_q=self._oddp_source_degree(source_degree),
            )
            self._universal_operations[key] = cached
        return cached

    def _odd_primary_operation_homogeneous(
        self,
        element: "PrimeFieldCohomologyElement",
        degree: int,
        r: int,
        *,
        bockstein: bool,
        algorithm: str,
    ) -> "PrimeFieldCohomologyElement":
        operation_degree = self.convention * r
        if operation_degree < 0:
            return self.zero()

        missing_vertices_per_factor = self._odd_primary_missing_vertices(
            operation_degree,
            bockstein=bockstein,
        )
        target_degree = degree + missing_vertices_per_factor
        target_data = self._degree_data.get(target_degree)
        if target_data is None:
            return self.zero()
        if operation_degree == 0 and not bockstein:
            return self.element({degree: element._coordinates.get(degree, {})})

        universal = self._universal_operation(
            operation_degree=operation_degree,
            source_degree=degree,
            bockstein=bockstein,
            target_degree=target_degree,
            missing_vertices_per_factor=missing_vertices_per_factor,
        )
        target_vector = cochain_operation_vector_from_universal(
            self.complex,
            self.cocycle(element, degree),
            universal,
            target_face_to_index=target_data.face_to_index,
            algorithm=algorithm,
        )
        return self.project_cocycle(target_degree, target_vector)

    def _build_boundary_columns(self) -> dict[int, list[Vector]]:
        columns: dict[int, list[Vector]] = {0: [{} for _ in self._faces.get(0, ())]}
        for degree in range(1, self.dimension + 1):
            lower_index = self._face_to_index[degree - 1]
            signs = tuple(((-1) ** index) % self.p for index in range(degree + 1))
            degree_columns = []
            for simplex in self._faces[degree]:
                column = {
                    lower_index[simplex[:index] + simplex[index + 1 :]]: signs[index]
                    for index in range(degree + 1)
                    if signs[index]
                }
                degree_columns.append(column)
            columns[degree] = degree_columns
        return columns

    def _build_coboundary_columns(self) -> dict[int, list[Vector]]:
        columns = {}
        for degree in range(self.dimension + 1):
            domain_faces = self._faces.get(degree, ())
            degree_columns = [{} for _ in domain_faces]
            for target_index, boundary in enumerate(
                self._boundary_columns.get(degree + 1, ())
            ):
                for domain_index, coefficient in boundary.items():
                    degree_columns[domain_index][target_index] = coefficient % self.p
            columns[degree] = degree_columns
        return columns

    def _build_degree_data(self) -> dict[int, _DegreeData]:
        image_basis_by_degree = {}
        cycles_by_degree = {}
        for degree in range(self.dimension + 1):
            image_basis_by_degree[degree], cycles_by_degree[degree] = (
                column_image_and_kernel_basis(
                    self._coboundary_columns[degree],
                    self.p,
                )
            )

        all_data = {}
        for degree in range(self.dimension + 1):
            faces = self._faces[degree]
            cycles = cycles_by_degree[degree]
            boundary_vectors = image_basis_by_degree.get(degree - 1, [])
            if self.reduced and degree == 0 and faces:
                boundary_vectors = [
                    *boundary_vectors,
                    {index: 1 for index in range(len(faces))},
                ]

            projector = CoordinateBasis(self.p)
            for boundary in boundary_vectors:
                projector.add_vector(boundary)

            cocycle_basis = []
            for cycle in sorted(cycles, key=_vector_sort_key):
                basis_coordinate = {len(cocycle_basis): 1}
                if projector.add(cycle, basis_coordinate):
                    cocycle_basis.append(cycle)

            all_data[degree] = _DegreeData(
                faces=faces,
                face_to_index=self._face_to_index[degree],
                cocycle_basis=tuple(cocycle_basis),
                projector=projector,
            )
        return all_data


class PrimeFieldCohomologyElement:
    """An element of a cohomology vector space over ``F_p``."""

    def __init__(self, parent: PrimeFieldCohomology, coordinates: dict[int, Vector]):
        self.parent = parent
        self._coordinates = {degree: dict(vector) for degree, vector in coordinates.items()}

    @property
    def p(self) -> int:
        """Return the characteristic."""
        return self.parent.p

    def sq(self, k: int) -> "PrimeFieldCohomologyElement":
        """Return ``Sq^k`` applied to this cohomology element."""
        return self.parent.square(self, k)

    def operation(
        self,
        r: int,
        *,
        bockstein: bool = False,
        algorithm: str = "auto",
    ) -> "PrimeFieldCohomologyElement":
        """Return the Steenrod operation selected by the coefficient prime."""
        return self.parent.operation(
            self,
            r,
            bockstein=bockstein,
            algorithm=algorithm,
        )

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

    def __add__(self, other: "PrimeFieldCohomologyElement") -> "PrimeFieldCohomologyElement":
        if self.parent is not other.parent:
            return NotImplemented
        coordinates = {degree: dict(vector) for degree, vector in self._coordinates.items()}
        for degree, vector in other._coordinates.items():
            coordinates[degree] = vector_add(coordinates.get(degree, {}), vector, self.p)
            if not coordinates[degree]:
                del coordinates[degree]
        return self.parent.element(coordinates)

    def __sub__(self, other: "PrimeFieldCohomologyElement") -> "PrimeFieldCohomologyElement":
        if self.parent is not other.parent:
            return NotImplemented
        coordinates = {degree: dict(vector) for degree, vector in self._coordinates.items()}
        for degree, vector in other._coordinates.items():
            coordinates[degree] = vector_add(coordinates.get(degree, {}), vector, self.p, -1)
            if not coordinates[degree]:
                del coordinates[degree]
        return self.parent.element(coordinates)

    def __rmul__(self, scalar: int) -> "PrimeFieldCohomologyElement":
        return self.parent.element({
            degree: vector_scale(vector, scalar, self.p)
            for degree, vector in self._coordinates.items()
        })

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, PrimeFieldCohomologyElement)
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
            for index, coefficient in sorted(self._coordinates[degree].items()):
                basis = f"h^{degree},{index}"
                if coefficient == 1:
                    terms.append(basis)
                else:
                    terms.append(f"{coefficient}*{basis}")
        return " + ".join(terms)


Mod2Cohomology = PrimeFieldCohomology
Mod2CohomologyElement = PrimeFieldCohomologyElement


def _codimension_one_faces(simplex: "Simplex", p: int) -> tuple[tuple[int, "Simplex"], ...]:
    return tuple(
        (((-1) ** index) % p, simplex[:index] + simplex[index + 1 :])
        for index in range(len(simplex))
    )


def _vector_sort_key(vector: Vector) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(vector.items()))
