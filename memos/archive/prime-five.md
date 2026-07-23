# Archived: prime-five computations

The prime-five phase keeps the degree-two source class but moves the target
from degree 6 to degree 10:

`P¹: H²(X; F₅) → H¹⁰(X; F₅)`.

For a degree-two class, the unstable identity is `P¹(x) = x⁵`. The package's
formula is native and self-contained: it uses the single fivefold tensor for
the top reduced power.

## Computed examples

| Model | Cells | Relevant cohomology | Result |
| --- | ---: | --- | --- |
| Sym⁵(S²) = CP⁵ | 22,010 | one class in degrees 0, 2, 4, 6, 8, 10 | rank-one P¹ |
| Sym⁵(T²) | 1,797,894 | `b₂ = 2`, `b₁₀ = 1` | rank-one P¹ |
| L¹¹(5) | 354,312 | one mod-5 class in every degree 0–11 | rank-one βP⁰ and P¹ |

The lens space is important as an independent validation: it enters through
dense face maps and a cyclic quotient rather than through the symmetric-set
constructor.

## Reference run

On the development arm64 macOS host with Python 3.14.3:

| Model | Build | Selected cohomology | Operation | Peak memory |
| --- | ---: | ---: | ---: | ---: |
| Sym⁵(S²) | 0.025 s | 0.202 s | 0.081 s | 29 MiB |
| Sym⁵(T²) | 1.333 s | 7.666 s | 8.837 s | 780 MiB |
| L¹¹(5) | 14.228 s | included in operations | 0.018 s βP⁰; 5.003 s P¹ | 1,509 MiB |

The symmetric-product operation time is measured after constructing the
source and target cohomology projectors. The complete cold focused pipeline
for Sym⁵(T²) is roughly 18 seconds.

Reproduce the standard run with:

```bash
python benchmarks/prime_five.py
```

The larger face-map validation is opt-in because of its memory use:

```bash
python benchmarks/prime_five.py --lens
```

## Scaling strategy

Three changes make the torus case feasible:

1. Degeneracy-capacity pruning avoids enumerating tuples that cannot become
   nondegenerate.
2. Symmetric cells are stored as compact unordered labels; their faces are
   computed on demand.
3. Cohomology builds boundary and coboundary data only for requested degrees.

The preflight method predicts the next model before allocating it:

```python
surface = spaces.orientable_surface(2)
assert sum(surface.symmetric_power_f_vector(5)) == 414_092_094
```

That cell count rules out a direct genus-two construction for now. The next
prime-five search should therefore inspect other compact ten-dimensional
families, screen their cohomology algebra first, and reserve cochain-level P¹
evaluation for models with plausible nonzero targets.
