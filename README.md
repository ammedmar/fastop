# fastop

`fastop` is a new package for computing cohomology operations on spaces presented simplicially.

The goal is to provide user-facing methods that start from a simplicial or
semi-simplicial object, compute cohomology over `F_p`, choose cohomology
classes, and evaluate Steenrod operations on those classes.

## Examples

```python
from fastop import spaces

H = spaces.real_projective_plane().cohomology(p=2)
x = H.basis(1)[0]
assert x.operation(1) == H.basis(2)[0]

CP2 = spaces.complex_projective_plane().cohomology(p=2)
u = CP2.basis(2)[0]
assert u.operation(2) == CP2.basis(4)[0]
```

Odd-primary operations are exposed on the same cohomology classes. Universal
formulas are served from the built-in catalog when available and otherwise
computed by the package's native formula engine. The installed package is
self-contained; the neighboring `oddp` project is used only as a development
oracle in optional parity tests. The resulting cocycle is projected back to
cohomology.

```python
H = spaces.complex_projective_space(3).cohomology(p=3)
u = H.basis(2)[0]
assert u.operation(1, algorithm="prime-three") == H.basis(6)[0]

M = spaces.moore_space(3).cohomology(p=3)
a = M.basis(1)[0]
assert a.operation(0, bockstein=True, algorithm="prime-three") == M.basis(2)[0]
```

Simplicial suspension can move a known operation into a range where it is not
a cup power:

```python
suspended = spaces.complex_projective_space(3).suspension(2)
H = suspended.cohomology(p=3)
assert H.operation_rank(4, 1) == 1
```

The completed implementation program for face-map input, compact quotients,
and the six-manifold search is recorded in [docs/roadmap.md](docs/roadmap.md).

Finite Delta-complexes use dense face maps compatible with Sage's
`DeltaComplex.cells()` output.  This permits loops, repeated faces, and other
compact semi-simplicial models:

```python
from fastop import DeltaComplex

circle = DeltaComplex([
    [()],       # one vertex
    [(0, 0)],   # one edge with both faces equal to that vertex
])
assert circle.cohomology(p=3).betti_numbers() == {0: 1, 1: 1}
```

Strict cell actions can be quotiented without subdivision. Actions may be
specified by graded permutations or by cell maps:

```python
from fastop import DeltaComplex, FiniteGroupAction

pentagon = DeltaComplex([
    [() for _ in range(5)],
    [((i + 1) % 5, i) for i in range(5)],
])
rotation = FiniteGroupAction.from_cell_maps(
    pentagon,
    lambda degree, cell: (cell + 1) % 5,
)
circle = pentagon.quotient(rotation, require_free=True)
assert circle.f_vector() == (1, 1)
```

The catalog lens space uses the same interface for the diagonal cyclic
quotient of a join of circles:

```python
L = spaces.lens_space(7, 3)
H = L.cohomology(p=3)
assert H.operation_rank(1, 0, bockstein=True) == 1
assert H.operation_rank(2, 1) == 1
```

Finite simplicial sets retain degeneracies until a complete face restriction
has been evaluated. This makes very small symmetric-product models possible:

```python
CP3 = spaces.minimal_simplicial_sphere(2).symmetric_power(3)
assert CP3.f_vector() == (1, 0, 3, 10, 25, 30, 15)
assert CP3.cohomology(p=3).operation_rank(2, 1) == 1

X = spaces.symmetric_product_of_surface(genus=1, power=3)
assert X.f_vector() == (1, 19, 126, 380, 572, 420, 120)
assert X.cohomology(p=3).operation_rank(2, 1) == 1

genus_two = spaces.symmetric_product_of_surface(genus=2)
assert sum(genus_two.f_vector()) == 41_478
assert genus_two.cohomology(p=3).operation_rank(2, 1) == 1
```

At the prime 5, the same construction reaches the fifth symmetric power.
The sphere is the projective-space ground truth, while the torus provides a
1.8-million-cell positive-genus computation:

```python
CP5 = spaces.symmetric_product_of_surface(genus=0, power=5)
assert CP5.cohomology(p=5).operation_rank(2, 1) == 1

torus_fifth = spaces.symmetric_product_of_surface(genus=1, power=5)
assert sum(torus_fifth.f_vector()) == 1_797_894
assert torus_fifth.cohomology(p=5).operation_rank(2, 1) == 1
```

Run `python benchmarks/showcase.py` for reproducible construction,
cohomology, and operation timings on all showcase models. A reference run and
the role of each example are recorded in [docs/showcase.md](docs/showcase.md).
The genus expansion and its extended benchmark are described in
[docs/symmetric-surfaces.md](docs/symmetric-surfaces.md). Prime-five results
and their performance envelope are recorded in
[docs/prime-five.md](docs/prime-five.md).

## Notebooks

- [Steenrod operations on symmetric products of surfaces](notebooks/symmetric_products_of_surfaces.ipynb)
  develops the genus family at primes 2, 3, 5, and 7.
- [Lens spaces](notebooks/lens_spaces.ipynb) separately develops the cyclic
  join quotients, their cell counts, Bocksteins, and reduced powers.

## Development

Install the package from a checkout with:

```bash
python -m pip install .
```

The C accelerator is optional; the package retains a pure-Python fallback if
it cannot be compiled. The supported public models, constructors, quotient
interface, and operation coverage are summarized in
[docs/package-api.md](docs/package-api.md).

The package-facing documentation is generated with Sphinx from the public
docstrings and narrative sources:

```bash
python -m pip install -e '.[docs]'
cd docs
make html
```

Open `docs/build/html/index.html` after the build.

Install the development dependencies and run the routine suite with:

```bash
python -m pip install -e '.[dev]'
python -m pytest
```

The two memory-intensive prime-five showcase regressions are opt-in:

```bash
python -m pytest -m large
```

GitHub Actions runs the routine suite on Python 3.10 through 3.14, builds the
Sphinx site with warnings treated as errors, and checks both distribution
artifacts. The large regressions remain an explicit local or release-machine
run.
