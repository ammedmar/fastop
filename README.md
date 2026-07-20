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
formulas are served from the built-in catalog when available, with the
neighboring `oddp` package used as a development fallback for formulas not
yet cataloged. The resulting cocycle is projected back to cohomology.

```python
H = spaces.complex_projective_space(3).cohomology(p=3)
u = H.basis(2)[0]
assert u.operation(1, algorithm="prime-three") == H.basis(6)[0]

M = spaces.moore_space(3).cohomology(p=3)
a = M.basis(1)[0]
assert a.operation(0, bockstein=True, algorithm="prime-three") == M.basis(2)[0]
```

The catalog also contains the matching complex $M_7$, a naturally
occurring odd-primary example, and simplicial suspension can move a known
operation into a range where it is not a cup power:

```python
M7 = spaces.matching_complex(7).cohomology(p=3)
assert M7.operation_rank(1, 0, bockstein=True) == 1

suspended = spaces.complex_projective_space(3).suspension(2)
H = suspended.cohomology(p=3)
assert H.operation_rank(4, 1) == 1
```

The implementation roadmap for face-map input, compact quotients, and the
six-manifold search is recorded in [docs/roadmap.md](docs/roadmap.md).

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

Strict cell actions can be quotiented without subdivision.  For example, the
catalog lens space is the diagonal cyclic quotient of a join of circles:

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

Run `python benchmarks/showcase.py` for reproducible construction,
cohomology, and operation timings on all showcase models. A reference run and
the role of each example are recorded in [docs/showcase.md](docs/showcase.md).
The genus expansion and its extended benchmark are described in
[docs/symmetric-surfaces.md](docs/symmetric-surfaces.md).

## Development

```bash
.venv/bin/python -m pytest
```
