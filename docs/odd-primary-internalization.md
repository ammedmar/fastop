# Odd-primary Steenrod operation internalization plan

This memo records how to internalize the useful parts of `oddp` into `fastop`.
It is not an implementation plan that preserves `oddp`'s module structure.
The goal is to keep the mathematics and tests, while replacing the current
verbose object model with smaller kernels that match `fastop`'s cohomology
classes and sparse vectors.

## Terminology

The public language in `fastop` should be mathematical:

- `operation(r)` over `F_2` means `Sq^r`.
- `operation(r)` over odd `F_p` means the reduced power `P^r`.
- `operation(r, bockstein=True)` over odd `F_p` means the Bockstein followed
  by `P^r`.

The names `support` and `prime-three` should not become `fastop` terminology.
They are names of two evaluation strategies in `oddp`, not mathematical
operations.

Better internal names:

- `all_targets`: evaluate the universal tensor formula on every target simplex.
  This corresponds to `oddp`'s `direct` algorithm.
- `support_enumeration`: build candidate target simplices from the support of
  the input cochain. This corresponds to `oddp`'s `support` algorithm.
- `triple_support_enumeration`: a specialization of support enumeration when
  `p = 3`, where candidate targets come from triples of source simplices. This
  corresponds to `oddp`'s `prime-three` algorithm.

The first one is simplest, but expensive. The second one is the general sparse
strategy. The third one is a p=3 optimization and should eventually disappear
behind the same internal sparse evaluator.

## What `oddp` currently does

The public operation entry point is:

```python
Steenrod.cochain_operation(complex, cochain, p, s, q, bockstein=False)
```

Its conventions differ from the `fastop` API:

- `q` is the negative cohomological degree.
- `s` is the negative reduced-power index.

For a `fastop` class in degree `d`, the bridge is:

```python
q = -d
s = -r
target_degree = d + 2*r*(p - 1) + int(bockstein)
```

The core pipeline in `oddp` is:

```text
operation parameters
  -> minimal resolution generator
  -> phi_dual
  -> psi_aux or psi_dual
  -> abc or abc_dual
  -> universal tensor/signature data
  -> simplicial evaluation on the cochain
```

The last step is the part that touches the input complex. Everything before it
is universal data depending only on `(p, r, source_degree, bockstein)`.

## What to keep

Keep the formulas and tests, not the class hierarchy.

Useful pieces:

- index conversion in `powers._operation_parameters`;
- `phi_dual`, especially its odd-prime coefficient convention;
- `psi_aux`, which gives the compact "signature" data used by sparse
  evaluation;
- `abc` and `abc_dual`, which translate between the tensor shape and the
  Milnor-resolution shape;
- the CP3 and Moore-space regression examples.

Pieces to avoid copying directly:

- the public `Steenrod` wrapper;
- the broad `ParametricCounter` arithmetic model;
- the resolution object hierarchy as public or semi-public API;
- the algorithm names `support` and `prime-three`.

`fastop` already has sparse vectors over `F_p`; native code should use plain
dicts with explicit modular arithmetic.

## Proposed internal shape

A native implementation should live behind the existing `operation(...)` API.
The likely private modules are:

```text
fastop/_odd_primary/
  indices.py
  universal.py
  evaluate.py
  reference.py
```

`indices.py` should own the convention conversion:

```python
OperationIndex(
    p: int,
    r: int,
    source_degree: int,
    bockstein: bool,
)
```

with computed fields:

```python
oddp_s = -r
oddp_q = -source_degree
target_degree = source_degree + 2*r*(p - 1) + int(bockstein)
missing_vertices_per_factor = 2*r*(p - 1) + int(bockstein)
```

`universal.py` should produce cached universal data:

```python
UniversalOperation(
    p: int,
    r: int,
    source_degree: int,
    bockstein: bool,
    target_degree: int,
    signatures: dict[tuple[tuple[int, ...], ...], int],
)
```

The cache key should be:

```python
(p, r, source_degree, bockstein)
```

`evaluate.py` should contain the sparse simplicial evaluator. It takes:

```python
faces_by_degree
cochain_support
universal_operation
```

and returns a sparse cochain in the target degree. It should not know about
cohomology projection.

`reference.py` can temporarily wrap `oddp` for comparison tests while native
kernels are being written.

## Evaluation strategies

### All-target evaluation

For every target simplex in the target degree, try every tensor summand in the
universal operation. This mirrors `oddp`'s `direct` algorithm.

This is useful as a simple reference implementation but should not be the main
path for large complexes.

### Sparse support enumeration

Choose `p` source simplices from the input cochain support. Their union is a
possible target simplex. If that union is actually present in the complex and
the omitted vertex positions have the required sizes, evaluate the universal
signature on that configuration.

This mirrors `oddp`'s `support` algorithm. In `fastop`, call it sparse support
enumeration or simply the sparse evaluator.

### Specialized p=3 sparse enumeration

When `p = 3`, the source tuple has exactly three factors. `oddp` uses a
candidate list to reduce triple enumeration. This is a good optimization idea,
but the name `prime-three` should not survive.

In `fastop`, this should be an internal branch inside the sparse evaluator:

```python
if p == 3:
    use_triple_candidate_enumeration(...)
else:
    use_general_support_enumeration(...)
```

The caller should still ask for `operation(r)`, not for an algorithm named
after the prime.

## Test oracle

Until native kernels replace the bridge, `oddp` should remain the oracle. The
initial native tests should compare:

- raw cochains from native all-target evaluation against `oddp` direct output;
- raw cochains from native sparse evaluation against `oddp` support output;
- projected cohomology classes through the current `operation(...)` API.

Regression spaces:

- CP2 over `F_3`: `P^1` on the degree-2 class is zero by dimension.
- CP3 over `F_3`: `P^1: H^2 -> H^6` has rank 1.
- Moore space over `F_3`: Bockstein `P^0: H^1 -> H^2` has rank 1.
- A future `p = 5` example before claiming confidence beyond p=3.

## Suggested next implementation slice

Do not port all of `oddp`.

First native slice:

1. Add `_odd_primary/indices.py`.
2. Add `_odd_primary/reference.py` as the current `oddp` bridge.
3. Move the existing cohomology-level bridge out of `cohomology.py` into that
   reference module.
4. Keep behavior unchanged.

Second native slice:

1. Implement universal signature extraction for a very small range of cases.
2. Compare it against `oddp`'s universal tensor chains.
3. Use CP3 and Moore only as end-to-end tests.

Third native slice:

1. Implement all-target evaluation.
2. Compare raw output with `oddp` direct output.
3. Then implement sparse support enumeration and compare raw output with the
   all-target evaluator.

Only after these slices should we consider C++ kernels.
