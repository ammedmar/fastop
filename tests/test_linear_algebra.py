import pytest

import fastop._linear_algebra as linear_algebra_module
from fastop._linear_algebra import (
    CoordinateBasis,
    column_image_and_kernel_basis,
    column_kernel_basis,
    nullspace,
    rank,
)


def _matrix_product(columns, vector, p):
    answer = {}
    for column_index, coefficient in vector.items():
        for row_index, value in columns[column_index].items():
            total = (answer.get(row_index, 0) + coefficient * value) % p
            if total:
                answer[row_index] = total
            else:
                answer.pop(row_index, None)
    return answer


def test_column_kernel_basis_matches_nullspace_rank_and_vectors():
    columns = [
        {0: 1, 1: 1},
        {1: 2, 2: 1},
        {0: 1, 1: 3, 2: 1},
        {3: 1},
        {3: 2},
    ]
    p = 5

    kernel = column_kernel_basis(columns, p)

    assert len(kernel) == len(nullspace(columns, 4, p))
    assert rank(kernel, p) == len(kernel)
    assert all(_matrix_product(columns, vector, p) == {} for vector in kernel)


def test_native_column_reduction_matches_python_fallback(monkeypatch):
    if linear_algebra_module._native_column_image_and_kernel_basis is None:
        pytest.skip("native extension is not built")

    columns = [
        {0: 1, 1: 1},
        {1: 2, 2: 1},
        {0: 1, 1: 3, 2: 1},
        {3: 1},
        {3: 2},
    ]

    native = column_image_and_kernel_basis(columns, 5)
    monkeypatch.setattr(
        linear_algebra_module,
        "_native_column_image_and_kernel_basis",
        None,
    )

    assert native == column_image_and_kernel_basis(columns, 5)


def test_coordinate_basis_uses_echelon_rows_for_coordinates():
    basis = CoordinateBasis(5)

    assert basis.add({2: 1, 0: 1}, {})
    assert basis.add({1: 1}, {0: 1})

    assert basis.coordinates({2: 3, 1: 4, 0: 3}) == {0: 4}
