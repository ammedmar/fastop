Mathematical conventions and scope
==================================

Coefficients and grading
------------------------

``fastop`` computes cohomology over the prime field :math:`\mathbf F_p`.
The public API uses cohomological grading and nonnegative operation indices.
For a class :math:`x\in H^d(X;\mathbf F_p)`:

* at :math:`p=2`, ``x.operation(k)`` is
  :math:`Sq^k(x)\in H^{d+k}(X;\mathbf F_2)`;
* at an odd prime, ``x.operation(r)`` is
  :math:`P^r(x)\in H^{d+2r(p-1)}(X;\mathbf F_p)`;
* ``x.operation(r, bockstein=True)`` is
  :math:`\beta P^r(x)\in H^{d+2r(p-1)+1}(X;\mathbf F_p)`.

In particular, ``operation(0, bockstein=True)`` is the mod-``p`` Bockstein.
The standard instability conditions are checked before cochain evaluation.

Cochains and representatives
-----------------------------

The chain basis consists of oriented simplices for ``SimplicialComplex`` and
indexed nondegenerate cells for ``DeltaComplex`` and ``SimplicialSet``.
Boundary signs use the usual alternating convention.  Simplicial-set chains
are normalized, so degenerate faces contribute zero.

A cohomology class is stored by coordinates in a deterministic basis.  Its
``cocycle()`` is one chosen representative of that class.  Neither the basis
nor the representative is a canonical geometric choice, and both may change
if the model, cell ordering, or reduction algorithm changes.  Equalities and
operation ranks are invariant; individual printed basis coordinates need not
be.

Products
--------

The product ``x * y`` is induced by the Alexander--Whitney diagonal and then
projected into the chosen cohomology basis.  ``x ** n`` uses this product, with
``x ** 0`` equal to the unit in unreduced cohomology.  Reduced cohomology has
no multiplicative unit.

Supported models and limitations
--------------------------------

The package accepts finite abstract simplicial complexes, finite
Delta-complexes, and finite normalized simplicial sets.  It does not currently
provide integral cohomology, arbitrary coefficient rings, infinite
simplicial objects, or implicit coercion from general Sage cell complexes.

The optional native extension changes performance but not the public API or
mathematical result.  Large symmetric powers and quotient models can still be
limited by the number of cells in the source and target degrees; preflight
their ``f_vector`` before construction.

References and citation
-----------------------

The odd-primary universal formulas are adapted from the ``oddp`` project and
the accompanying work of Federico Cantero-Moran and Anibal M.
Medina-Mardones, *Steenrod operations and Tate resolutions*.  The mod-two
evaluation is based on the cup-``i`` approach developed for fast computation
of Steenrod squares.  See the repository's ``CITATION.cff`` and
``THIRD_PARTY_NOTICES.md`` for citation and attribution information.
