"""Execute the lightweight introductory notebooks as Python workflows."""

from __future__ import annotations

import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "notebook_name",
    [
        "data_structures_and_cohomology.ipynb",
        "products_quotients_and_symmetric_powers.ipynb",
    ],
)
def test_introductory_notebook(notebook_name):
    notebook_path = ROOT / "notebooks" / notebook_name
    notebook = json.loads(notebook_path.read_text())
    namespace = {"__name__": "__main__"}

    for cell_number, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        exec(
            compile(
                source,
                f"{notebook_path.name}:cell-{cell_number}",
                "exec",
            ),
            namespace,
        )
