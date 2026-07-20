# Odd-primary showcase

The showcase combines familiar validation cases with compact models that
exercise the new input formats. Every entry below is constructed and tested
through the public `fastop` API.

| Model | Nondegenerate cells | Mod-3 Betti numbers | Detected operation |
| --- | ---: | --- | --- |
| Matching complex M₇ | 231 | `(1, 1, 21)` | rank-one βP⁰: H¹ → H² |
| ΣCP³ | 26,654 | `(1, 0, 0, 1, 0, 1, 0, 1)` | rank-one P¹: H³ → H⁷ |
| Σ²CP³ | 79,964 | `(1, 0, 0, 0, 1, 0, 1, 0, 1)` | rank-one P¹: H⁴ → H⁸ |
| L⁷(3) | 800 | `(1, 1, 1, 1, 1, 1, 1, 1)` | rank-one βP⁰ and P¹ |
| Sym³(S²) = CP³ | 84 | `(1, 0, 1, 0, 1, 0, 1)` | rank-one P¹: H² → H⁶ |
| Sym³(T²) | 1,638 | `(1, 2, 2, 2, 2, 2, 1)` | rank-one P¹: H² → H⁶ |

The first four cases provide ground truth and a suspension stress test.
Sym³(S²) validates the simplicial-set product and quotient pipeline against
CP³. Sym³(T²) is the more interesting example in the wild: a naturally
occurring six-manifold, built as a compact quotient, in which the computation
finds a nontrivial odd-primary operation.

## Reference timings

The following are medians of five runs on the development arm64 macOS host
with Python 3.14.3. Times are reference measurements rather than performance
guarantees.

| Model | Build | Cohomology | Warm operation |
| --- | ---: | ---: | ---: |
| M₇ | 0.000191 s | 0.000679 s | 0.000070 s |
| ΣCP³ | 0.019261 s | 0.135493 s | 0.001237 s |
| Σ²CP³ | 0.091220 s | 0.467672 s | 0.015326 s |
| L⁷(3) | 0.011986 s | 0.003781 s | 0.000100–0.000161 s |
| Sym³(S²) | 0.032990 s | 0.000212 s | 0.000129 s |
| Sym³(T²) | 0.565474 s | 0.016900 s | 0.000285 s |

Reproduce the run from the repository root with:

```bash
python benchmarks/showcase.py --repeats 5
```
