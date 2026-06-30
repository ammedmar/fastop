from fastop import spaces
from fastop._cochain_evaluation import evaluate_source_mod_2


def test_source_mod_2_evaluator_detects_basic_square_target():
    assert evaluate_source_mod_2(
        3,
        [(0, 1), (1, 2)],
        {(0, 1, 2)},
    ) == {(0, 1, 2)}


def test_source_mod_2_evaluator_ignores_non_targets():
    assert evaluate_source_mod_2(
        3,
        [(0, 1), (1, 2)],
        set(),
    ) == set()


def test_source_mod_2_evaluator_matches_cp2_square_projection():
    cohomology = spaces.complex_projective_plane().cohomology()
    source = cohomology.basis(2)[0]
    target_data = cohomology._degree_data[4]
    source_data = cohomology._degree_data[2]
    cocycle_vector = cohomology.cocycle_vector(source, 2)
    support = [source_data.faces[index] for index in cocycle_vector]

    target_support = evaluate_source_mod_2(5, support, set(target_data.faces))
    target_vector = {
        target_data.face_to_index[simplex]: 1
        for simplex in target_support
    }

    assert cohomology.project_cocycle(4, target_vector) == cohomology.basis(4)[0]
