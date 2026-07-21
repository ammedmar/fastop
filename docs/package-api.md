# Package API

`fastop` accepts three finite combinatorial models and gives each one the same
cohomology interface.

| Model | Use it for | Essential input |
| --- | --- | --- |
| `SimplicialComplex` | ordinary triangulations | maximal simplices |
| `DeltaComplex` | quotients, loops, and repeated faces | dense face-index tables |
| `SimplicialSet` | degeneracies, products, and symmetric powers | nondegenerate simplices with face references |

All three provide `cohomology(p=...)`. The result provides `basis`,
`betti_number`, `betti_numbers`, `operation_matrix`, and `operation_rank`.
Cohomology classes provide `operation`.

## Constructors for the showcase families

```python
from fastop import spaces

surface = spaces.orientable_surface(genus=2)
symmetry_quotient = spaces.symmetric_product_of_surface(genus=2, power=3)
lens = spaces.lens_space(dimension=11, order=5)
projective_plane = spaces.nonorientable_surface(crosscaps=1)
sym3_projective_plane = spaces.symmetric_product_of_nonorientable_surface(
    crosscaps=1,
    power=3,
)
```

`symmetric_product_of_surface` uses a compact model of the closed orientable
surface of the requested genus. `minimal_simplicial_surface` exposes that
base model directly. Before constructing a large symmetric
power, use `surface.symmetric_power_f_vector(power)` to inspect its cell count.

`nonorientable_surface(crosscaps)` uses the polygon presentation
`a1 a1 ... ah ah`; one crosscap is the minimal simplicial-set model of
`RP2`, and two crosscaps give the Klein bottle. Its symmetric-product
constructor uses the same degree-lazy representation as the oriented family.

The two example families deliberately remain separate. Symmetric products are
formed as permutation quotients of products of a normalized simplicial-set
surface. Lens spaces are diagonal cyclic quotients of joins of polygonal
circles and are returned as Delta-complexes.

## Face-map Delta-complexes

For a `d`-cell `sigma`, the value `face_maps[d][sigma][i]` is the index of
`d_i(sigma)` among the cells of degree `d - 1`. Vertices have empty face
tuples. This dense format preserves distinct cells with identical boundaries,
so it can express group quotients without subdivision.

```python
from fastop import DeltaComplex

circle = DeltaComplex([
    [()],
    [(0, 0)],
])
```

`DeltaComplex.from_sage(value)` copies the tables exposed by a finite Sage
Delta-complex. The resulting object has no runtime dependency on Sage.

## Finite group actions and quotients

An action generator is a permutation in every cell degree. The convenient
cell-map constructor computes those permutations from `(degree, cell)`
callables:

```python
from fastop import DeltaComplex, FiniteGroupAction

order = 5
polygon = DeltaComplex([
    [() for _ in range(order)],
    [((i + 1) % order, i) for i in range(order)],
])
rotation = FiniteGroupAction.from_cell_maps(
    polygon,
    lambda degree, cell: (cell + 1) % order,
)

assert rotation.order(polygon) == 5
assert rotation.is_free(polygon)
quotient = polygon.quotient(rotation, require_free=True)
```

For already-computed graded permutations, use
`FiniteGroupAction.cyclic(generator)` or pass multiple generators to
`FiniteGroupAction((generator_1, generator_2))`. A quotient always validates
that every generator is a permutation and commutes with every face map.
`require_free=True` additionally enumerates the generated finite group and
rejects any nonidentity element that fixes a cell.

The same action object is accepted by `SimplicialSet.quotient`.

## Operation coverage

At `p=2`, `operation(k)` computes the Steenrod square `Sq^k`. At an odd prime,
`operation(r)` computes the reduced power `P^r`, while
`operation(r, bockstein=True)` computes `beta P^r`. Instability and a missing
source or target cohomology group are detected before cochain evaluation.

The installed package is self-contained. It uses cached formulas for a few
frequent low-dimensional cases and computes general odd-primary universal
formulas with its native builder. The sibling `oddp` project is retained only
as a development oracle for optional parity tests; it is not an installation
dependency of `fastop`.

The optional C extension accelerates sparse linear algebra and cochain
evaluation. If it is unavailable, the same public API uses pure Python.

## Test tiers

`python -m pytest` runs the routine correctness suite. The 1.8-million-cell
Sym⁵(T²) and high-memory L¹¹(5) validations are retained as regressions but are
excluded by default. Run them explicitly with `python -m pytest -m large` on a
machine with sufficient memory.
