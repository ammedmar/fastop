Finite quotients
================

Use ``FiniteGroupAction`` when your model is naturally a quotient of a finite
action. An action generator is one permutation of the cells in every degree.
The cell-map constructor is usually clearer than writing those permutations by
hand.

.. code-block:: python

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
   circle = polygon.quotient(rotation, require_free=True)

``quotient`` verifies that each generated map is a permutation and commutes
with every face map. ``require_free=True`` additionally enumerates the finite
group and rejects an action with a fixed cell.

If you already have the graded permutations, construct a cyclic action
directly:

.. code-block:: python

   rotation = FiniteGroupAction.cyclic(
       ((1, 2, 3, 4, 0), (1, 2, 3, 4, 0))
   )

The first tuple acts on vertices; the second acts on edges. For several
generators, pass a tuple of such graded permutations to ``FiniteGroupAction``.
Both ``DeltaComplex`` and ``SimplicialSet`` accept the action object in their
``quotient`` method.
