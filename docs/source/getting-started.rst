Getting started
===============

Install from a checkout:

.. code-block:: console

   python -m pip install .

The optional native extension accelerates selected computations. If it cannot
be compiled, the pure-Python implementation retains the same public API.

For development, including generated documentation:

.. code-block:: console

   python -m pip install -e '.[dev]'
   cd docs
   make html

The HTML site is written to ``docs/build/html``. It is generated output and is
not committed.

The standard test suite is lightweight enough for routine development:

.. code-block:: console

   python -m pytest

The two high-memory prime-five validations are deliberately opt-in:

.. code-block:: console

   python -m pytest -m large
