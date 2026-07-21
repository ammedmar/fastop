Releasing
=========

The first planned public version is ``0.1.0a1``. Releases use
``.github/workflows/release.yml`` and PyPI Trusted Publishing, so no long-lived
package-index token belongs in GitHub secrets.

One-time setup
--------------

Before the first release:

#. Create or identify the GitHub repository and add it as the local ``origin``.
#. Create GitHub environments named ``testpypi`` and ``pypi``. Require manual
   approval for the ``pypi`` environment.
#. Register pending Trusted Publishers on PyPI and TestPyPI with the GitHub
   owner, repository, workflow filename ``release.yml``, and the corresponding
   environment name.
#. Protect tags matching ``v*`` so only maintainers can create release tags.

The project name ``fastop`` was unclaimed on both package indexes when the
alpha release infrastructure was prepared. Index availability must still be
checked immediately before registering the pending publishers.

TestPyPI candidate
------------------

Set ``fastop.__version__`` to a non-development version and run the ``Release``
workflow manually. A manual run builds and tests binary wheels for CPython
3.10 through 3.14 on Linux, Windows, Intel macOS, and Apple Silicon macOS,
builds the source distribution, validates the artifacts, and publishes them
only to TestPyPI.

Install the candidate in a clean environment using both the wheel path and the
source path:

.. code-block:: console

   python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ fastop==0.1.0a1
   python -c "from fastop import spaces; assert spaces.lens_space(7, 3).cohomology(p=3).operation_rank(2, 1) == 1"

   python -m pip install --no-binary fastop --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ fastop==0.1.0a1

TestPyPI cannot replace an uploaded file. Bump the prerelease version before
uploading another candidate.

PyPI alpha
----------

After the TestPyPI candidate passes:

#. Confirm CI is green on the exact commit.
#. Create and publish a GitHub prerelease tagged ``v0.1.0a1``.
#. Approve the protected ``pypi`` environment deployment.
#. Install ``fastop==0.1.0a1`` from PyPI in clean Python 3.10 and 3.14
   environments and run the public smoke example.

The release event refuses to publish when the Git tag does not exactly match
``v`` followed by the package version. The publishing jobs receive only
``id-token: write`` and contain only artifact download and publishing steps.
