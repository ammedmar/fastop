"""Finite groups acting by permutations on graded simplicial cells."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

CellAction = tuple[tuple[int, ...], ...]
CellMap = Callable[[int, int], int]

__all__ = ["FiniteGroupAction"]


@dataclass(frozen=True)
class FiniteGroupAction:
    """A finite cell action specified by a set of permutation generators.

    Each generator contains one permutation for every cell degree.  The
    quotient validates that these permutations commute with all face maps.
    """

    generators: tuple[CellAction, ...]

    def __post_init__(self) -> None:
        normalized = _normalize_raw_actions(self.generators)
        if not normalized:
            raise ValueError("a finite group action needs at least one generator")
        object.__setattr__(self, "generators", normalized)

    @classmethod
    def cyclic(cls, generator: Iterable[Iterable[int]]) -> "FiniteGroupAction":
        """Create a cyclic action from one graded cell permutation."""
        return cls((_normalize_raw_action(generator),))

    @classmethod
    def from_cell_maps(
        cls,
        model,
        *generators: CellMap,
    ) -> "FiniteGroupAction":
        """Create generators from callables ``(degree, cell) -> image``."""
        if not generators:
            raise ValueError("a finite group action needs at least one generator")
        counts = model.f_vector()
        return cls(
            tuple(
                tuple(
                    tuple(generator(degree, cell) for cell in range(cell_count))
                    for degree, cell_count in enumerate(counts)
                )
                for generator in generators
            )
        )

    def order(self, model) -> int:
        """Return the order of the generated permutation group."""
        counts = model.f_vector()
        validate_permutations(self.generators, counts)
        return len(generated_actions(self.generators, counts))

    def is_free(self, model) -> bool:
        """Return whether every nonidentity element is fixed-cell free."""
        counts = model.f_vector()
        validate_permutations(self.generators, counts)
        try:
            validate_free_action(self.generators, counts)
        except ValueError:
            return False
        return True


def normalize_actions(
    generators: FiniteGroupAction | Iterable[Iterable[Iterable[int]]],
) -> tuple[CellAction, ...]:
    if isinstance(generators, FiniteGroupAction):
        return generators.generators
    return _normalize_raw_actions(generators)


def _normalize_raw_actions(
    generators: Iterable[Iterable[Iterable[int]]],
) -> tuple[CellAction, ...]:
    return tuple(_normalize_raw_action(generator) for generator in generators)


def _normalize_raw_action(
    generator: Iterable[Iterable[int]],
) -> CellAction:
    return tuple(tuple(permutation) for permutation in generator)


def validate_permutations(
    actions: tuple[CellAction, ...],
    counts: tuple[int, ...],
) -> None:
    for action in actions:
        if len(action) != len(counts):
            raise ValueError("a cell action needs one permutation per dimension")
        for permutation, cell_count in zip(action, counts):
            if (
                len(permutation) != cell_count
                or set(permutation) != set(range(cell_count))
            ):
                raise ValueError("each cell action must be a permutation")


def generated_actions(
    generators: tuple[CellAction, ...],
    counts: tuple[int, ...],
) -> set[CellAction]:
    identity = tuple(tuple(range(count)) for count in counts)
    actions = {identity}
    pending = [identity]
    while pending:
        current = pending.pop()
        for generator in generators:
            composite = tuple(
                tuple(
                    generator[degree][current[degree][cell]]
                    for cell in range(counts[degree])
                )
                for degree in range(len(counts))
            )
            if composite not in actions:
                actions.add(composite)
                pending.append(composite)
    return actions


def validate_free_action(
    generators: tuple[CellAction, ...],
    counts: tuple[int, ...],
) -> None:
    for action in generated_actions(generators, counts):
        if is_identity_action(action):
            continue
        if any(
            image == cell
            for permutation in action
            for cell, image in enumerate(permutation)
        ):
            raise ValueError("the generated action is not free on cells")


def cell_orbits(
    generators: tuple[CellAction, ...],
    counts: tuple[int, ...],
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    orbit_indices = []
    orbit_representatives = []
    for degree, cell_count in enumerate(counts):
        parent = list(range(cell_count))

        def find(cell: int) -> int:
            while parent[cell] != cell:
                parent[cell] = parent[parent[cell]]
                cell = parent[cell]
            return cell

        def union(left: int, right: int) -> None:
            left_root = find(left)
            right_root = find(right)
            if left_root == right_root:
                return
            if left_root < right_root:
                parent[right_root] = left_root
            else:
                parent[left_root] = right_root

        for action in generators:
            for cell, image in enumerate(action[degree]):
                union(cell, image)

        roots = tuple(find(cell) for cell in range(cell_count))
        representatives = tuple(sorted(set(roots)))
        root_to_orbit = {
            root: orbit for orbit, root in enumerate(representatives)
        }
        orbit_indices.append(tuple(root_to_orbit[root] for root in roots))
        orbit_representatives.append(representatives)
    return tuple(orbit_indices), tuple(orbit_representatives)


def is_identity_action(action: CellAction) -> bool:
    return all(
        all(cell == image for cell, image in enumerate(permutation))
        for permutation in action
    )
