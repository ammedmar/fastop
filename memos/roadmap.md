# Roadmap

`fastop` now has a compact public interface for finite simplicial models,
prime-field cohomology, cup products, and Steenrod operations. Completed
investigations and benchmark reports live in `memos/archive`; the notebooks
are the canonical worked examples.

## Before the first public release

- Publish the rendered Sphinx documentation.
- Check the mathematical conventions and citations against the accompanying
  papers.
- Decide whether the advanced `algorithm`, `formula_source`, and `convention`
  controls should remain public or move to a diagnostic interface.
- Validate the TestPyPI candidate on the oldest and newest supported Python
  versions.

## Structural cleanup

- Move embedded triangulation tables out of `spaces.py` into private catalog
  data.
- Separate the lazy symmetric-power implementation from the core
  `SimplicialSet` source file while retaining one conceptual simplicial-set
  interface.
- Add a small internal protocol describing the common finite-model operations
  used by cohomology.

## Mathematical development

- Continue the prime-five search with compact ten-dimensional models whose
  cell counts and cohomology rings can be screened before cochain evaluation.
- Treat the prime-seven projective-space computation as an opt-in scaling
  benchmark rather than a routine regression.
- Add further Bockstein families only when their coefficient sequences and
  conventions can be represented explicitly in the public API.
