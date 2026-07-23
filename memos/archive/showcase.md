# Archived: odd-primary showcase

The showcase combines familiar validation cases with compact models that
exercise the new input formats. Every entry below is constructed and tested
through the public `fastop` API.

| Model | Nondegenerate cells | Mod-3 Betti numbers | Detected operation |
| --- | ---: | --- | --- |
| ΣCP³ | 26,654 | `(1, 0, 0, 1, 0, 1, 0, 1)` | rank-one P¹: H³ → H⁷ |
| Σ²CP³ | 79,964 | `(1, 0, 0, 0, 1, 0, 1, 0, 1)` | rank-one P¹: H⁴ → H⁸ |
| L⁷(3) | 800 | `(1, 1, 1, 1, 1, 1, 1, 1)` | rank-one βP⁰ and P¹ |
| Sym³(S²) = CP³ | 84 | `(1, 0, 1, 0, 1, 0, 1)` | rank-one P¹: H² → H⁶ |
| Sym³(T²) | 1,638 | `(1, 2, 2, 2, 2, 2, 1)` | rank-one P¹: H² → H⁶ |
| Sym³(Σ₂) | 41,478 | `(1, 4, 7, 8, 7, 4, 1)` | rank-one P¹: H² → H⁶ |
| Sym⁵(S²) = CP⁵ | 22,010 | one class in even degrees 0–10 | rank-one mod-5 P¹ |
| Sym⁵(T²) | 1,797,894 | measured `b₂=2`, `b₁₀=1` | rank-one mod-5 P¹ |

The first three cases provide ground truth and a suspension stress test.
Sym³(S²) validates the simplicial-set product and quotient pipeline against
CP³. The surface sequence Sym³(Σ_g) then holds the dimension and operation
fixed while increasing both topology and computational complexity. Genus two
is covered by the regression suite; genus three is an explicit extended run.

## Reference timings

The following are medians of five runs on the development arm64 macOS host
with Python 3.14.3. Times are reference measurements rather than performance
guarantees.

| Model | Build | Cohomology | Warm operation |
| --- | ---: | ---: | ---: |
| ΣCP³ | 0.021379 s | 0.147362 s | 0.001332 s |
| Σ²CP³ | 0.098188 s | 0.507693 s | 0.016326 s |
| L⁷(3) | 0.013023 s | 0.004228 s | 0.000111–0.000172 s |
| Sym³(S²) | 0.008316 s | 0.000245 s | 0.000142 s |
| Sym³(T²) | 0.129659 s | 0.016254 s | 0.000746 s |

Reproduce the run from the repository root with:

```bash
python benchmarks/showcase.py --repeats 5
```

Direct construction of unordered symmetric simplices reduced the median
Sym³(T²) build from 0.565 seconds in the first implementation to 0.130
seconds. Run the expanding genus family separately with:

```bash
python benchmarks/symmetric_surfaces.py
python benchmarks/symmetric_surfaces.py 3
```

The prime-five runs use degree-lazy cohomology and compact symmetric cell
labels. Reproduce them with:

```bash
python benchmarks/prime_five.py
python benchmarks/prime_five.py --lens
```
