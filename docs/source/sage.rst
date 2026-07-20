Working with Sage
=================

``fastop`` is compatible with Sage's *input ideas*, but it is not a drop-in
replacement for Sage's complete topology API. It is a focused finite-model
package: construct or import a model, compute mod-``p`` cohomology, and
evaluate operations.

What is shared
--------------

Both projects use familiar notions such as facets, cells, dimension, faces,
and cohomology. A ``fastop`` object has ``dimension``, ``cells`` and
``cohomology`` methods, so simple computational workflows feel similar.

What differs
------------

Sage is a large interactive algebra system with many constructors, homology
backends, categories, and conversions. ``fastop`` does not reproduce that
ecosystem or accept arbitrary Sage objects by implicit coercion. Its objects
are immutable, compact, and intentionally narrow: enough combinatorial data
to evaluate the supported cohomology operations efficiently.

Import an abstract Sage simplicial complex
------------------------------------------

For an ordinary Sage complex, extract its facets and pass them to the
``fastop`` constructor. The only requirement is that each facet can be
iterated as integer vertex labels.

.. code-block:: python

   from fastop import SimplicialComplex

   # sage_complex is a Sage SimplicialComplex.
   facets = [tuple(facet) for facet in sage_complex.facets()]
   model = SimplicialComplex(facets)

   H = model.cohomology(p=3)
   print(H.betti_numbers())

Import Sage face-map models
---------------------------

``DeltaComplex.from_sage`` copies the dense table returned by a finite Sage
Delta-complex's ``cells()`` method. ``SimplicialSet.from_sage`` copies finite
Sage simplicial sets, including degenerate faces. Neither result keeps a Sage
object alive, so Sage is not a runtime dependency after conversion.

.. code-block:: python

   from fastop import DeltaComplex, SimplicialSet

   delta_model = DeltaComplex.from_sage(sage_delta_complex)
   simplicial_set_model = SimplicialSet.from_sage(sage_simplicial_set)

   assert delta_model.cohomology(p=3).betti_numbers()

The adapters expect the corresponding finite Sage object. If the object is a
general cell complex, first choose an explicit finite simplicial,
Delta-complex, or simplicial-set representation rather than relying on an
implicit conversion.
