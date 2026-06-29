"""Native evaluators for odd-primary universal operation data."""

from __future__ import annotations

from typing import TYPE_CHECKING

from fastop._odd_primary.universal import UniversalOperation

if TYPE_CHECKING:
    from fastop.simplicial import Simplex


def evaluate_all_targets(
    target_faces: set["Simplex"] | frozenset["Simplex"],
    cochain: dict["Simplex", int],
    universal: UniversalOperation,
) -> dict["Simplex", int]:
    """Evaluate every target simplex against every universal tensor term."""
    answer: dict["Simplex", int] = {}
    for target in target_faces:
        coefficient = 0
        for tensor, tensor_coefficient in universal.terms.items():
            product = tensor_coefficient
            for factor in tensor:
                source = tuple(target[index] for index in factor)
                source_coefficient = cochain.get(source)
                if source_coefficient is None:
                    break
                product *= source_coefficient
            else:
                coefficient += product
        coefficient %= universal.p
        if coefficient:
            answer[target] = coefficient
    return answer
