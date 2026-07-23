# fastop

[![CI](https://github.com/ammedmar/fastop/actions/workflows/ci.yml/badge.svg)](https://github.com/ammedmar/fastop/actions/workflows/ci.yml)

`fastop` computes Steenrod operations on prime-field cohomology for finite
simplicial models.

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

Begin with the
[mathematical and package overview](https://github.com/ammedmar/fastop/blob/main/docs/source/overview.rst),
then choose an [input model](https://github.com/ammedmar/fastop/blob/main/docs/source/input-models.rst)
and see the [operation guide](https://github.com/ammedmar/fastop/blob/main/docs/source/operations.rst).
The [notebooks](https://github.com/ammedmar/fastop/tree/main/notebooks) contain
the longer examples on symmetric products of surfaces and lens spaces.

## Development

```bash
python -m pip install -e '.[dev]'
python -m pytest
```

Build the HTML documentation with `cd docs && make html`. Memory-intensive
regressions are available separately through `python -m pytest -m large`.

`fastop` is distributed under the
[MIT license](https://github.com/ammedmar/fastop/blob/main/LICENSE).

Authors: Federico Cantero-Morán and Anibal M. Medina-Mardones.
