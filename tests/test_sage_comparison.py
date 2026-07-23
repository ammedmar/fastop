import json
import os
import shutil
import subprocess
import tempfile
import textwrap

import pytest

from fastop import SimplicialComplex, spaces


EXAMPLES = [
    ("circle", [(0, 1), (0, 2), (1, 2)], 1, 1),
    ("sphere", list(SimplicialComplex.simplex_boundary(2).facets), 2, 1),
    ("real_projective_plane", list(spaces.real_projective_space(2).facets), 1, 1),
    ("real_projective_3_space", list(spaces.real_projective_space(3).facets), 1, 1),
    ("complex_projective_plane", list(spaces.complex_projective_space(2).facets), 2, 2),
]

COHOMOLOGY_EXAMPLES = [
    ("sphere", list(SimplicialComplex.simplex_boundary(2).facets)),
    ("real_projective_plane", list(spaces.real_projective_space(2).facets)),
    ("real_projective_3_space", list(spaces.real_projective_space(3).facets)),
    ("complex_projective_plane", list(spaces.complex_projective_space(2).facets)),
]


def test_operation_ranks_match_sage_when_available():
    sage = shutil.which("sage")
    if sage is None:
        pytest.skip("Sage is not installed")

    expected = _sage_operation_ranks(sage)
    actual = {
        name: SimplicialComplex(facets).cohomology().operation_rank(degree, square)
        for name, facets, degree, square in EXAMPLES
    }
    assert actual == expected


def test_mod_p_betti_numbers_match_sage_when_available():
    sage = shutil.which("sage")
    if sage is None:
        pytest.skip("Sage is not installed")

    primes = (2, 3, 5)
    expected = _sage_betti_numbers(sage, primes)
    actual = {
        f"{name}_p{p}": SimplicialComplex(facets).cohomology(p=p).betti_numbers()
        for name, facets in COHOMOLOGY_EXAMPLES
        for p in primes
    }
    assert actual == expected


def _sage_operation_ranks(sage):
    payload = [
        {"name": name, "facets": facets, "degree": degree, "square": square}
        for name, facets, degree, square in EXAMPLES
    ]
    script = textwrap.dedent(
        f"""
        import json
        from sage.all import GF, Matrix, SimplicialComplex

        payload = {payload!r}
        ranks = {{}}
        for item in payload:
            K = SimplicialComplex(item["facets"])
            K.set_immutable()
            H = K.cohomology_ring(GF(2))
            degree = item["degree"]
            square = item["square"]
            target_degree = degree + square
            target = list(H.basis(target_degree))
            target_indices = [basis.support()[0] for basis in target]
            rows = []
            for basis in H.basis(degree):
                image = basis.Sq(square).monomial_coefficients()
                rows.append([int(image.get(index, 0)) for index in target_indices])
            ranks[item["name"]] = Matrix(GF(2), rows).rank()
        print(json.dumps(ranks))
        """
    )
    with tempfile.TemporaryDirectory(prefix="fastop-sage-") as dot_sage:
        env = os.environ.copy()
        env["DOT_SAGE"] = dot_sage
        result = subprocess.run(
            [sage, "-c", script],
            check=False,
            capture_output=True,
            text=True,
            env=env,
        )
    if result.returncode:
        pytest.skip(f"Sage is present but not usable: {result.stderr.strip()}")
    return json.loads(result.stdout)


def _sage_betti_numbers(sage, primes):
    payload = [
        {"name": name, "facets": facets}
        for name, facets in COHOMOLOGY_EXAMPLES
    ]
    script = textwrap.dedent(
        f"""
        import json
        from sage.all import GF, SimplicialComplex

        payload = {payload!r}
        primes = {primes!r}
        betti = {{}}
        for item in payload:
            K = SimplicialComplex(item["facets"])
            K.set_immutable()
            for p in primes:
                H = K.cohomology_ring(GF(p))
                betti[f"{{item['name']}}_p{{p}}"] = {{
                    degree: len(H.basis(degree))
                    for degree in range(K.dimension() + 1)
                    if len(H.basis(degree))
                }}
        print(json.dumps(betti))
        """
    )
    with tempfile.TemporaryDirectory(prefix="fastop-sage-") as dot_sage:
        env = os.environ.copy()
        env["DOT_SAGE"] = dot_sage
        result = subprocess.run(
            [sage, "-c", script],
            check=False,
            capture_output=True,
            text=True,
            env=env,
        )
    if result.returncode:
        pytest.skip(f"Sage is present but not usable: {result.stderr.strip()}")
    return {
        key: {int(degree): value for degree, value in betti.items()}
        for key, betti in json.loads(result.stdout).items()
    }
