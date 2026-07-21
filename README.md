# fastop

[![CI](https://github.com/ammedmar/fastop/actions/workflows/ci.yml/badge.svg)](https://github.com/ammedmar/fastop/actions/workflows/ci.yml)

`fastop` computes cohomology and Steenrod operations on finite simplicial
models. It provides a common prime-field cohomology interface for simplicial
complexes, Delta-complexes, and simplicial sets.

The package includes:

- cohomology over $\mathbf F_p$, with cocycle representatives and cup products;
- Steenrod squares and odd-primary operations $P^r$ and $\beta P^r$;
- products, suspensions, finite group actions, quotients, and symmetric powers;
- compact face-map models suited to identifications and group-action quotients.

`fastop` has no required runtime dependencies. An optional C extension
accelerates core computations when a compiler is available.

## Quick start

Install from a checkout:

```bash
python -m pip install .
```

All supported input models use the same cohomology API:

```python
from fastop import spaces

H = spaces.complex_projective_space(3).cohomology(p=3)
u = H.basis(2)[0]

assert u.operation(1) == u**3  # P^1(u) = u^3
```

## Documentation and examples

The [Sphinx documentation](https://github.com/ammedmar/fastop/blob/main/docs/source/index.rst)
describes the input models, cohomology interface, group actions, quotients,
and Sage interoperability.

The [notebooks](https://github.com/ammedmar/fastop/tree/main/notebooks) provide
hands-on introductions and develop the larger computations on symmetric
products of surfaces and lens spaces.

## Development

```bash
python -m pip install -e '.[dev]'
python -m pytest
```

Build the HTML documentation with `cd docs && make html`. Memory-intensive
regressions are available separately through `python -m pytest -m large`.

`fastop` is distributed under the
[MIT license](https://github.com/ammedmar/fastop/blob/main/LICENSE).
