"""Cohomology over finite prime fields."""

from __future__ import annotations

from dataclasses import dataclass
from fastop._cochain_evaluation import (
    cochain_operation_vector_from_universal,
    evaluate_all_targets_mod_2,
    evaluate_source_mod_2,
)
from fastop._linear_algebra import (
    CoordinateBasis,
    Vector,
    add_to_basis,
    clean_vector,
    column_image_and_kernel_basis,
    coordinate_basis_from_vectors,
    is_prime,
    vector_add,
    vector_scale,
)
from fastop._universal import universal_operation

@dataclass(frozen=True)
class _DegreeData:
    faces: tuple["Simplex", ...]
    face_to_index: dict["Simplex", int]
    cocycle_basis: tuple[Vector, ...]
    projector: CoordinateBasis


class _DegreeDataCache(dict[int, _DegreeData]):
    """Lazy degree-data cache for a cohomology object."""

    def __init__(self, parent: "PrimeFieldCohomology"):
        super().__init__()
        self._parent = parent

    def __missing__(self, degree: int) -> _DegreeData:
        data = self._parent._build_degree_data_for_degree(degree)
        if data is None:
            raise KeyError(degree)
        self[degree] = data
        return data

    def get(self, degree: int, default=None):
        if isinstance(degree, int) and 0 <= degree <= self._parent.dimension:
            return self[degree]
        return default


class PrimeFieldCohomology:
    """Cohomology of a supported finite simplicial model over ``F_p``."""

    def __init__(
        self,
        complex_,
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
        self._faces = {}
        self._face_to_index = {}
        self._boundary_columns: dict[int, list[Vector]] = {}
        self._coboundary_columns: dict[int, list[Vector]] = {}
        self._image_basis_by_degree: dict[int, list[Vector]] = {}
        self._cycles_by_degree: dict[int, list[Vector]] = {}
        self._degree_data = _DegreeDataCache(self)
        self._universal_operations = {}

    def _ensure_faces(self, degree: int) -> None:
        if degree in self._faces or degree < 0 or degree > self.dimension:
            return
        cells = self.complex.cells(degree)
        faces = (
            cells
            if isinstance(cells, range) and cells == range(len(cells))
            else tuple(sorted(cells))
        )
        self._faces[degree] = faces
        self._face_to_index[degree] = (
            faces
            if isinstance(faces, range)
            else {face: index for index, face in enumerate(faces)}
        )

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

    def one(self) -> "PrimeFieldCohomologyElement":
        """Return the multiplicative unit in unreduced cohomology."""
        if self.reduced:
            raise ValueError("reduced cohomology has no multiplicative unit")
        self._ensure_faces(0)
        return self.project_cocycle(
            0,
            {index: 1 for index in range(len(self._faces[0]))},
        )

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

    def cup_product(
        self,
        left: "PrimeFieldCohomologyElement",
        right: "PrimeFieldCohomologyElement",
    ) -> "PrimeFieldCohomologyElement":
        """Return the Alexander--Whitney cup product ``left * right``."""
        if left.parent is not self or right.parent is not self:
            raise ValueError("cup-product factors must belong to this cohomology object")
        answer = self.zero()
        for left_degree in sorted(left._coordinates):
            for right_degree in sorted(right._coordinates):
                answer += self._cup_product_homogeneous(
                    left,
                    left_degree,
                    right,
                    right_degree,
                )
        return answer

    def cup_product_matrix(self, left_degree: int, right_degree: int) -> list[Vector]:
        """Return cup-product columns for lexicographically ordered basis pairs.

        The column for ``(i, j)`` records the product of basis class ``i`` in
        ``left_degree`` with basis class ``j`` in ``right_degree``.
        """
        target_degree = left_degree + right_degree
        return [
            (
                left * right
            )._coordinates.get(target_degree, {})
            for left in self.basis(left_degree)
            for right in self.basis(right_degree)
        ]

    def _cup_product_homogeneous(
        self,
        left: "PrimeFieldCohomologyElement",
        left_degree: int,
        right: "PrimeFieldCohomologyElement",
        right_degree: int,
    ) -> "PrimeFieldCohomologyElement":
        target_degree = left_degree + right_degree
        target_data = self._degree_data.get(target_degree)
        if target_data is None:
            return self.zero()

        left_cochain = self.cocycle(left, left_degree)
        right_cochain = self.cocycle(right, right_degree)
        front_positions = tuple(range(left_degree + 1))
        back_positions = tuple(range(left_degree, target_degree + 1))
        product_vector = {}
        for target_index, target in enumerate(target_data.faces):
            front = self.complex.restrict(
                target_degree,
                target,
                front_positions,
            )
            if front is None:
                continue
            back = self.complex.restrict(
                target_degree,
                target,
                back_positions,
            )
            if back is None:
                continue
            coefficient = (
                left_cochain.get(front, 0)
                * right_cochain.get(back, 0)
            ) % self.p
            if coefficient:
                product_vector[target_index] = coefficient
        return self.project_cocycle(target_degree, product_vector)

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
        formula_source: str = "auto",
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
                formula_source=formula_source,
            )
        return answer

    def operation_matrix(
        self,
        degree: int,
        r: int,
        *,
        bockstein: bool = False,
        algorithm: str = "auto",
        formula_source: str = "auto",
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
                formula_source=formula_source,
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
        formula_source: str = "auto",
    ) -> int:
        """Return the rank of the selected Steenrod operation from ``H^degree``."""
        self._validate_operation_index(r)
        if not isinstance(bockstein, bool):
            raise TypeError("bockstein must be a boolean")
        if self.p == 2 and bockstein:
            raise NotImplementedError("Bockstein operations require an odd prime")
        if not self._satisfies_instability_bound(degree, r, bockstein=bockstein):
            return 0

        source_rank = self.betti_number(degree)
        target_degree = self._operation_target_degree(degree, r, bockstein=bockstein)
        target_rank = self.betti_number(target_degree)
        if source_rank == 0 or target_rank == 0:
            return 0

        image_basis: dict[int, Vector] = {}
        for index in range(source_rank):
            basis_element = self.element({degree: {index: 1}})
            image = basis_element.operation(
                r,
                bockstein=bockstein,
                algorithm=algorithm,
                formula_source=formula_source,
            )._coordinates.get(target_degree, {})
            if add_to_basis(image_basis, image, self.p) and len(image_basis) == target_rank:
                return target_rank
        return len(image_basis)

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
        if self.complex._supports_vertex_algorithms:
            target_support = evaluate_source_mod_2(
                target_degree + 1,
                support,
                set(target_data.faces),
            )
        else:
            cochain = {
                source_data.faces[index]: coefficient
                for index, coefficient in cocycle_vector.items()
            }
            target_support = evaluate_all_targets_mod_2(
                self.complex,
                target_data.faces,
                cochain,
                source_degree=degree,
                target_degree=target_degree,
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

    def _satisfies_instability_bound(self, degree: int, r: int, *, bockstein: bool) -> bool:
        operation_degree = self.convention * r
        if operation_degree < 0:
            return False
        if self.p == 2:
            return not bockstein and operation_degree <= degree
        return 2 * operation_degree + int(bockstein) <= degree

    def _universal_operation(
        self,
        *,
        operation_degree: int,
        source_degree: int,
        bockstein: bool,
        target_degree: int,
        missing_vertices_per_factor: int,
        formula_source: str,
    ):
        key = (
            operation_degree,
            source_degree,
            bockstein,
            target_degree,
            missing_vertices_per_factor,
            formula_source,
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
                formula_source=formula_source,
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
        formula_source: str,
    ) -> "PrimeFieldCohomologyElement":
        operation_degree = self.convention * r
        if operation_degree < 0:
            return self.zero()
        if 2 * operation_degree + int(bockstein) > degree:
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
            formula_source=formula_source,
        )
        target_vector = cochain_operation_vector_from_universal(
            self.complex,
            self.cocycle(element, degree),
            universal,
            target_face_to_index=target_data.face_to_index,
            algorithm=algorithm,
        )
        if operation_degree % 2:
            target_vector = vector_scale(target_vector, -1, self.p)
        return self.project_cocycle(target_degree, target_vector)

    def _boundary_columns_for_degree(self, degree: int) -> list[Vector]:
        cached = self._boundary_columns.get(degree)
        if cached is not None:
            return cached
        if degree < 0 or degree > self.dimension:
            return []
        self._ensure_faces(degree)
        if degree == 0:
            columns = [{} for _ in self._faces.get(0, ())]
        else:
            self._ensure_faces(degree - 1)
            lower_index = self._face_to_index[degree - 1]
            signs = tuple(((-1) ** index) % self.p for index in range(degree + 1))
            columns = []
            for simplex in self._faces[degree]:
                column = {}
                for index, sign in enumerate(signs):
                    face = self.complex.face(degree, simplex, index)
                    if face is None:
                        continue
                    face_index = lower_index[face]
                    coefficient = (column.get(face_index, 0) + sign) % self.p
                    if coefficient:
                        column[face_index] = coefficient
                    elif face_index in column:
                        del column[face_index]
                columns.append(column)
        self._boundary_columns[degree] = columns
        return columns

    def _coboundary_columns_for_degree(self, degree: int) -> list[Vector]:
        cached = self._coboundary_columns.get(degree)
        if cached is not None:
            return cached
        if degree < 0 or degree > self.dimension:
            return []
        self._ensure_faces(degree)
        columns = [{} for _ in self._faces.get(degree, ())]
        for target_index, boundary in enumerate(
            self._boundary_columns_for_degree(degree + 1)
        ):
            for domain_index, coefficient in boundary.items():
                columns[domain_index][target_index] = coefficient % self.p
        self._coboundary_columns[degree] = columns
        return columns

    def _image_and_cycles(self, degree: int) -> tuple[list[Vector], list[Vector]]:
        cached_image = self._image_basis_by_degree.get(degree)
        cached_cycles = self._cycles_by_degree.get(degree)
        if cached_image is not None and cached_cycles is not None:
            return cached_image, cached_cycles
        image_basis, cycles = column_image_and_kernel_basis(
            self._coboundary_columns_for_degree(degree),
            self.p,
        )
        self._image_basis_by_degree[degree] = image_basis
        self._cycles_by_degree[degree] = cycles
        return image_basis, cycles

    def _build_degree_data_for_degree(self, degree: int) -> _DegreeData | None:
        if degree < 0 or degree > self.dimension:
            return None
        self._ensure_faces(degree)
        _, cycles = self._image_and_cycles(degree)
        if degree == 0:
            boundary_vectors = []
        else:
            boundary_vectors, _ = self._image_and_cycles(degree - 1)
        if self.reduced and degree == 0 and self._faces[degree]:
            boundary_vectors = [
                *boundary_vectors,
                {index: 1 for index in range(len(self._faces[degree]))},
            ]

        cocycle_basis, projector = coordinate_basis_from_vectors(
            self.p,
            boundary_vectors,
            sorted(cycles, key=_vector_sort_key),
        )

        return _DegreeData(
            faces=self._faces[degree],
            face_to_index=self._face_to_index[degree],
            cocycle_basis=tuple(cocycle_basis),
            projector=projector,
        )

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
        formula_source: str = "auto",
    ) -> "PrimeFieldCohomologyElement":
        """Return the Steenrod operation selected by the coefficient prime."""
        return self.parent.operation(
            self,
            r,
            bockstein=bockstein,
            algorithm=algorithm,
            formula_source=formula_source,
        )

    def cup(
        self,
        other: "PrimeFieldCohomologyElement",
    ) -> "PrimeFieldCohomologyElement":
        """Return the cup product with ``other``."""
        return self.parent.cup_product(self, other)

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
        if not isinstance(scalar, int) or isinstance(scalar, bool):
            return NotImplemented
        return self.parent.element({
            degree: vector_scale(vector, scalar, self.p)
            for degree, vector in self._coordinates.items()
        })

    def __mul__(self, other):
        if isinstance(other, int) and not isinstance(other, bool):
            return other * self
        if isinstance(other, PrimeFieldCohomologyElement):
            if self.parent is not other.parent:
                return NotImplemented
            return self.parent.cup_product(self, other)
        return NotImplemented

    def __pow__(self, exponent: int) -> "PrimeFieldCohomologyElement":
        if not isinstance(exponent, int) or isinstance(exponent, bool):
            raise TypeError("cohomology-class exponent must be an integer")
        if exponent < 0:
            raise ValueError("cohomology-class exponent must be nonnegative")
        answer = self.parent.one()
        factor = self
        while exponent:
            if exponent & 1:
                answer = answer * factor
            exponent >>= 1
            if exponent:
                factor = factor * factor
        return answer

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


def _vector_sort_key(vector: Vector) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(vector.items()))
