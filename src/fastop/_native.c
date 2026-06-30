#define PY_SSIZE_T_CLEAN
#include <Python.h>


static int
add_mod(PyObject *answer, PyObject *key, long value, long p)
{
    PyObject *old_value = PyDict_GetItemWithError(answer, key);
    if (old_value == NULL && PyErr_Occurred()) {
        return -1;
    }

    long total = value % p;
    if (total < 0) {
        total += p;
    }
    if (old_value != NULL) {
        long old = PyLong_AsLong(old_value);
        if (old == -1 && PyErr_Occurred()) {
            return -1;
        }
        total = (old + total) % p;
    }

    if (total == 0) {
        if (old_value != NULL && PyDict_DelItem(answer, key) < 0) {
            return -1;
        }
        return 0;
    }

    PyObject *new_value = PyLong_FromLong(total);
    if (new_value == NULL) {
        return -1;
    }
    int status = PyDict_SetItem(answer, key, new_value);
    Py_DECREF(new_value);
    return status;
}


static PyObject *
selected_face(PyObject *target, PyObject *factor)
{
    Py_ssize_t length = PyTuple_GET_SIZE(factor);
    PyObject *source = PyTuple_New(length);
    if (source == NULL) {
        return NULL;
    }

    for (Py_ssize_t i = 0; i < length; i++) {
        PyObject *index_object = PyTuple_GET_ITEM(factor, i);
        Py_ssize_t index = PyLong_AsSsize_t(index_object);
        if (index == -1 && PyErr_Occurred()) {
            Py_DECREF(source);
            return NULL;
        }
        PyObject *vertex = PyTuple_GetItem(target, index);
        if (vertex == NULL) {
            Py_DECREF(source);
            return NULL;
        }
        Py_INCREF(vertex);
        PyTuple_SET_ITEM(source, i, vertex);
    }

    return source;
}


static int
tuple_contains_ssize(PyObject *tuple, Py_ssize_t value)
{
    Py_ssize_t length = PyTuple_GET_SIZE(tuple);
    for (Py_ssize_t i = 0; i < length; i++) {
        Py_ssize_t item = PyLong_AsSsize_t(PyTuple_GET_ITEM(tuple, i));
        if (item == -1 && PyErr_Occurred()) {
            return -1;
        }
        if (item == value) {
            return 1;
        }
    }
    return 0;
}


static PyObject *
selected_positions_from_omissions(PyObject *omissions, Py_ssize_t target_length)
{
    PyObject *selected = PyTuple_New(target_length - PyTuple_GET_SIZE(omissions));
    if (selected == NULL) {
        return NULL;
    }

    Py_ssize_t selected_index = 0;
    for (Py_ssize_t position = 0; position < target_length; position++) {
        int omitted = tuple_contains_ssize(omissions, position);
        if (omitted < 0) {
            Py_DECREF(selected);
            return NULL;
        }
        if (!omitted) {
            PyObject *position_object = PyLong_FromSsize_t(position);
            if (position_object == NULL) {
                Py_DECREF(selected);
                return NULL;
            }
            PyTuple_SET_ITEM(selected, selected_index, position_object);
            selected_index++;
        }
    }
    return selected;
}


static PyObject *
selected_positions_for_pattern(PyObject *pattern, Py_ssize_t target_length)
{
    if (!PyTuple_Check(pattern) || PyTuple_GET_SIZE(pattern) != 3) {
        PyErr_SetString(PyExc_TypeError, "p=3 patterns must be triples");
        return NULL;
    }

    PyObject *selected = PyTuple_New(3);
    if (selected == NULL) {
        return NULL;
    }
    for (Py_ssize_t i = 0; i < 3; i++) {
        PyObject *omissions = PyTuple_GET_ITEM(pattern, i);
        if (!PyTuple_Check(omissions)) {
            Py_DECREF(selected);
            PyErr_SetString(PyExc_TypeError, "omissions must be tuples");
            return NULL;
        }
        PyObject *factor_positions = selected_positions_from_omissions(
            omissions,
            target_length
        );
        if (factor_positions == NULL) {
            Py_DECREF(selected);
            return NULL;
        }
        PyTuple_SET_ITEM(selected, i, factor_positions);
    }
    return selected;
}


static int
pair_union_size(PyObject *positions_a, PyObject *positions_b)
{
    Py_ssize_t size = PyTuple_GET_SIZE(positions_a);
    Py_ssize_t length_b = PyTuple_GET_SIZE(positions_b);
    for (Py_ssize_t i = 0; i < length_b; i++) {
        Py_ssize_t position = PyLong_AsSsize_t(PyTuple_GET_ITEM(positions_b, i));
        if (position == -1 && PyErr_Occurred()) {
            return -1;
        }
        int found = tuple_contains_ssize(positions_a, position);
        if (found < 0) {
            return -1;
        }
        if (!found) {
            size++;
        }
    }
    return (int)size;
}


static int
best_anchor_pair(PyObject *selected, Py_ssize_t *anchor_a, Py_ssize_t *anchor_b, Py_ssize_t *remaining)
{
    static const Py_ssize_t pairs[3][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};
    int best_size = -1;
    for (Py_ssize_t i = 0; i < 3; i++) {
        PyObject *positions_a = PyTuple_GET_ITEM(selected, pairs[i][0]);
        PyObject *positions_b = PyTuple_GET_ITEM(selected, pairs[i][1]);
        int size = pair_union_size(positions_a, positions_b);
        if (size < 0) {
            return -1;
        }
        if (size > best_size) {
            best_size = size;
            *anchor_a = pairs[i][0];
            *anchor_b = pairs[i][1];
            *remaining = pairs[i][2];
        }
    }
    return 0;
}


static PyObject *
remaining_fixed_data(PyObject *positions_c, PyObject *positions_a, PyObject *positions_b)
{
    PyObject *fixed = PyTuple_New(2);
    if (fixed == NULL) {
        return NULL;
    }
    PyObject *fixed_positions = PyList_New(0);
    PyObject *fixed_indices = PyList_New(0);
    if (fixed_positions == NULL || fixed_indices == NULL) {
        Py_XDECREF(fixed_positions);
        Py_XDECREF(fixed_indices);
        Py_DECREF(fixed);
        return NULL;
    }

    Py_ssize_t length = PyTuple_GET_SIZE(positions_c);
    for (Py_ssize_t i = 0; i < length; i++) {
        PyObject *position_object = PyTuple_GET_ITEM(positions_c, i);
        Py_ssize_t position = PyLong_AsSsize_t(position_object);
        if (position == -1 && PyErr_Occurred()) {
            Py_DECREF(fixed_positions);
            Py_DECREF(fixed_indices);
            Py_DECREF(fixed);
            return NULL;
        }
        int in_a = tuple_contains_ssize(positions_a, position);
        if (in_a < 0) {
            Py_DECREF(fixed_positions);
            Py_DECREF(fixed_indices);
            Py_DECREF(fixed);
            return NULL;
        }
        int in_b = tuple_contains_ssize(positions_b, position);
        if (in_b < 0) {
            Py_DECREF(fixed_positions);
            Py_DECREF(fixed_indices);
            Py_DECREF(fixed);
            return NULL;
        }
        if (in_a || in_b) {
            PyObject *index_object = PyLong_FromSsize_t(i);
            if (index_object == NULL) {
                Py_DECREF(fixed_positions);
                Py_DECREF(fixed_indices);
                Py_DECREF(fixed);
                return NULL;
            }
            int status = PyList_Append(fixed_positions, position_object);
            if (status == 0) {
                status = PyList_Append(fixed_indices, index_object);
            }
            Py_DECREF(index_object);
            if (status < 0) {
                Py_DECREF(fixed_positions);
                Py_DECREF(fixed_indices);
                Py_DECREF(fixed);
                return NULL;
            }
        }
    }

    PyObject *fixed_positions_tuple = PyList_AsTuple(fixed_positions);
    PyObject *fixed_indices_tuple = PyList_AsTuple(fixed_indices);
    Py_DECREF(fixed_positions);
    Py_DECREF(fixed_indices);
    if (fixed_positions_tuple == NULL || fixed_indices_tuple == NULL) {
        Py_XDECREF(fixed_positions_tuple);
        Py_XDECREF(fixed_indices_tuple);
        Py_DECREF(fixed);
        return NULL;
    }
    PyTuple_SET_ITEM(fixed, 0, fixed_positions_tuple);
    PyTuple_SET_ITEM(fixed, 1, fixed_indices_tuple);
    return fixed;
}


static PyObject *
source_key_from_indices(PyObject *source, PyObject *indices)
{
    Py_ssize_t length = PyTuple_GET_SIZE(indices);
    PyObject *key = PyTuple_New(length);
    if (key == NULL) {
        return NULL;
    }
    for (Py_ssize_t i = 0; i < length; i++) {
        Py_ssize_t index = PyLong_AsSsize_t(PyTuple_GET_ITEM(indices, i));
        if (index == -1 && PyErr_Occurred()) {
            Py_DECREF(key);
            return NULL;
        }
        PyObject *vertex = PyTuple_GET_ITEM(source, index);
        Py_INCREF(vertex);
        PyTuple_SET_ITEM(key, i, vertex);
    }
    return key;
}


static PyObject *
build_source_index(PyObject *support, PyObject *indices)
{
    PyObject *index = PyDict_New();
    if (index == NULL) {
        return NULL;
    }

    Py_ssize_t support_length = PySequence_Fast_GET_SIZE(support);
    PyObject **support_items = PySequence_Fast_ITEMS(support);
    for (Py_ssize_t i = 0; i < support_length; i++) {
        PyObject *pair = support_items[i];
        PyObject *source = PyTuple_GET_ITEM(pair, 0);
        PyObject *key = source_key_from_indices(source, indices);
        if (key == NULL) {
            Py_DECREF(index);
            return NULL;
        }

        PyObject *bucket = PyDict_GetItemWithError(index, key);
        if (bucket == NULL && PyErr_Occurred()) {
            Py_DECREF(key);
            Py_DECREF(index);
            return NULL;
        }
        if (bucket == NULL) {
            bucket = PyList_New(0);
            if (bucket == NULL) {
                Py_DECREF(key);
                Py_DECREF(index);
                return NULL;
            }
            if (PyDict_SetItem(index, key, bucket) < 0) {
                Py_DECREF(bucket);
                Py_DECREF(key);
                Py_DECREF(index);
                return NULL;
            }
            Py_DECREF(bucket);
        }
        if (PyList_Append(bucket, pair) < 0) {
            Py_DECREF(key);
            Py_DECREF(index);
            return NULL;
        }
        Py_DECREF(key);
    }

    return index;
}


static int
assigned_positions_are_increasing(PyObject *target)
{
    PyObject *last = NULL;
    Py_ssize_t length = PyList_GET_SIZE(target);
    for (Py_ssize_t i = 0; i < length; i++) {
        PyObject *vertex = PyList_GET_ITEM(target, i);
        if (vertex == Py_None) {
            continue;
        }
        if (last != NULL) {
            int increasing = PyObject_RichCompareBool(last, vertex, Py_LT);
            if (increasing <= 0) {
                return increasing;
            }
        }
        last = vertex;
    }
    return 1;
}


static PyObject *
partial_target_from_pair(
    PyObject *source_a,
    PyObject *source_b,
    PyObject *positions_a,
    PyObject *positions_b,
    Py_ssize_t target_length
)
{
    PyObject *target = PyList_New(target_length);
    if (target == NULL) {
        return NULL;
    }
    for (Py_ssize_t i = 0; i < target_length; i++) {
        Py_INCREF(Py_None);
        PyList_SET_ITEM(target, i, Py_None);
    }

    Py_ssize_t length_a = PyTuple_GET_SIZE(positions_a);
    for (Py_ssize_t i = 0; i < length_a; i++) {
        Py_ssize_t position = PyLong_AsSsize_t(PyTuple_GET_ITEM(positions_a, i));
        if (position == -1 && PyErr_Occurred()) {
            Py_DECREF(target);
            return NULL;
        }
        PyObject *vertex = PyTuple_GET_ITEM(source_a, i);
        PyObject *old_vertex = PyList_GET_ITEM(target, position);
        Py_DECREF(old_vertex);
        Py_INCREF(vertex);
        PyList_SET_ITEM(target, position, vertex);
    }

    Py_ssize_t length_b = PyTuple_GET_SIZE(positions_b);
    for (Py_ssize_t i = 0; i < length_b; i++) {
        Py_ssize_t position = PyLong_AsSsize_t(PyTuple_GET_ITEM(positions_b, i));
        if (position == -1 && PyErr_Occurred()) {
            Py_DECREF(target);
            return NULL;
        }
        PyObject *vertex = PyTuple_GET_ITEM(source_b, i);
        PyObject *old_vertex = PyList_GET_ITEM(target, position);
        if (old_vertex != Py_None) {
            int equal = PyObject_RichCompareBool(old_vertex, vertex, Py_EQ);
            if (equal < 0) {
                Py_DECREF(target);
                return NULL;
            }
            if (!equal) {
                Py_DECREF(target);
                Py_RETURN_NONE;
            }
        }
        Py_INCREF(vertex);
        Py_DECREF(old_vertex);
        PyList_SET_ITEM(target, position, vertex);
    }

    int increasing = assigned_positions_are_increasing(target);
    if (increasing < 0) {
        Py_DECREF(target);
        return NULL;
    }
    if (!increasing) {
        Py_DECREF(target);
        Py_RETURN_NONE;
    }
    return target;
}


static PyObject *
key_from_partial_target(PyObject *partial_target, PyObject *positions)
{
    Py_ssize_t length = PyTuple_GET_SIZE(positions);
    PyObject *key = PyTuple_New(length);
    if (key == NULL) {
        return NULL;
    }
    for (Py_ssize_t i = 0; i < length; i++) {
        Py_ssize_t position = PyLong_AsSsize_t(PyTuple_GET_ITEM(positions, i));
        if (position == -1 && PyErr_Occurred()) {
            Py_DECREF(key);
            return NULL;
        }
        PyObject *vertex = PyList_GET_ITEM(partial_target, position);
        Py_INCREF(vertex);
        PyTuple_SET_ITEM(key, i, vertex);
    }
    return key;
}


static PyObject *
target_from_remaining_source(PyObject *partial_target, PyObject *source, PyObject *positions)
{
    PyObject *target = PySequence_List(partial_target);
    if (target == NULL) {
        return NULL;
    }

    Py_ssize_t length = PyTuple_GET_SIZE(positions);
    for (Py_ssize_t i = 0; i < length; i++) {
        Py_ssize_t position = PyLong_AsSsize_t(PyTuple_GET_ITEM(positions, i));
        if (position == -1 && PyErr_Occurred()) {
            Py_DECREF(target);
            return NULL;
        }
        PyObject *vertex = PyTuple_GET_ITEM(source, i);
        PyObject *old_vertex = PyList_GET_ITEM(target, position);
        if (old_vertex != Py_None) {
            int equal = PyObject_RichCompareBool(old_vertex, vertex, Py_EQ);
            if (equal < 0) {
                Py_DECREF(target);
                return NULL;
            }
            if (!equal) {
                Py_DECREF(target);
                Py_RETURN_NONE;
            }
        }
        Py_INCREF(vertex);
        Py_DECREF(old_vertex);
        PyList_SET_ITEM(target, position, vertex);
    }

    Py_ssize_t target_length = PyList_GET_SIZE(target);
    for (Py_ssize_t i = 0; i < target_length; i++) {
        if (PyList_GET_ITEM(target, i) == Py_None) {
            Py_DECREF(target);
            Py_RETURN_NONE;
        }
    }
    for (Py_ssize_t i = 0; i < target_length - 1; i++) {
        int increasing = PyObject_RichCompareBool(
            PyList_GET_ITEM(target, i),
            PyList_GET_ITEM(target, i + 1),
            Py_LT
        );
        if (increasing < 0) {
            Py_DECREF(target);
            return NULL;
        }
        if (!increasing) {
            Py_DECREF(target);
            Py_RETURN_NONE;
        }
    }

    PyObject *target_tuple = PyList_AsTuple(target);
    Py_DECREF(target);
    return target_tuple;
}


static PyObject *
evaluate_all_targets(PyObject *self, PyObject *args)
{
    PyObject *target_faces;
    PyObject *cochain;
    PyObject *terms;
    long p;

    if (!PyArg_ParseTuple(args, "OOlO", &target_faces, &cochain, &p, &terms)) {
        return NULL;
    }
    if (!PyDict_Check(cochain)) {
        PyErr_SetString(PyExc_TypeError, "cochain must be a dictionary");
        return NULL;
    }
    if (!PyDict_Check(terms)) {
        PyErr_SetString(PyExc_TypeError, "terms must be a dictionary");
        return NULL;
    }
    if (p <= 1) {
        PyErr_SetString(PyExc_ValueError, "p must be greater than one");
        return NULL;
    }

    PyObject *answer = PyDict_New();
    if (answer == NULL) {
        return NULL;
    }

    PyObject *target_iter = PyObject_GetIter(target_faces);
    if (target_iter == NULL) {
        Py_DECREF(answer);
        return NULL;
    }

    PyObject *target;
    while ((target = PyIter_Next(target_iter)) != NULL) {
        if (!PyTuple_Check(target)) {
            Py_DECREF(target);
            Py_DECREF(target_iter);
            Py_DECREF(answer);
            PyErr_SetString(PyExc_TypeError, "target faces must be tuples");
            return NULL;
        }

        long coefficient = 0;
        PyObject *factor_cache = PyDict_New();
        if (factor_cache == NULL) {
            Py_DECREF(target);
            Py_DECREF(target_iter);
            Py_DECREF(answer);
            return NULL;
        }

        PyObject *tensor;
        PyObject *tensor_coefficient_object;
        Py_ssize_t position = 0;
        while (PyDict_Next(terms, &position, &tensor, &tensor_coefficient_object)) {
            if (!PyTuple_Check(tensor)) {
                Py_DECREF(factor_cache);
                Py_DECREF(target);
                Py_DECREF(target_iter);
                Py_DECREF(answer);
                PyErr_SetString(PyExc_TypeError, "term tensors must be tuples");
                return NULL;
            }

            long term_value = PyLong_AsLong(tensor_coefficient_object);
            if (term_value == -1 && PyErr_Occurred()) {
                Py_DECREF(factor_cache);
                Py_DECREF(target);
                Py_DECREF(target_iter);
                Py_DECREF(answer);
                return NULL;
            }

            int complete = 1;
            Py_ssize_t factor_count = PyTuple_GET_SIZE(tensor);
            for (Py_ssize_t i = 0; i < factor_count; i++) {
                PyObject *factor = PyTuple_GET_ITEM(tensor, i);
                if (!PyTuple_Check(factor)) {
                    Py_DECREF(factor_cache);
                    Py_DECREF(target);
                    Py_DECREF(target_iter);
                    Py_DECREF(answer);
                    PyErr_SetString(PyExc_TypeError, "tensor factors must be tuples");
                    return NULL;
                }

                PyObject *source_coefficient = PyDict_GetItemWithError(factor_cache, factor);
                if (source_coefficient == NULL && PyErr_Occurred()) {
                    Py_DECREF(factor_cache);
                    Py_DECREF(target);
                    Py_DECREF(target_iter);
                    Py_DECREF(answer);
                    return NULL;
                }
                if (source_coefficient == NULL) {
                    PyObject *source = selected_face(target, factor);
                    if (source == NULL) {
                        Py_DECREF(factor_cache);
                        Py_DECREF(target);
                        Py_DECREF(target_iter);
                        Py_DECREF(answer);
                        return NULL;
                    }
                    source_coefficient = PyDict_GetItemWithError(cochain, source);
                    if (source_coefficient == NULL && PyErr_Occurred()) {
                        Py_DECREF(source);
                        Py_DECREF(factor_cache);
                        Py_DECREF(target);
                        Py_DECREF(target_iter);
                        Py_DECREF(answer);
                        return NULL;
                    }
                    if (source_coefficient == NULL) {
                        source_coefficient = Py_None;
                    }
                    if (PyDict_SetItem(factor_cache, factor, source_coefficient) < 0) {
                        Py_DECREF(source);
                        Py_DECREF(factor_cache);
                        Py_DECREF(target);
                        Py_DECREF(target_iter);
                        Py_DECREF(answer);
                        return NULL;
                    }
                    Py_DECREF(source);
                }

                if (source_coefficient == Py_None) {
                    complete = 0;
                    break;
                }

                long source_value = PyLong_AsLong(source_coefficient);
                if (source_value == -1 && PyErr_Occurred()) {
                    Py_DECREF(factor_cache);
                    Py_DECREF(target);
                    Py_DECREF(target_iter);
                    Py_DECREF(answer);
                    return NULL;
                }
                term_value *= source_value;
            }

            if (complete) {
                coefficient += term_value;
            }
        }

        coefficient %= p;
        if (coefficient < 0) {
            coefficient += p;
        }
        if (coefficient != 0 && add_mod(answer, target, coefficient, p) < 0) {
            Py_DECREF(factor_cache);
            Py_DECREF(target);
            Py_DECREF(target_iter);
            Py_DECREF(answer);
            return NULL;
        }

        Py_DECREF(factor_cache);
        Py_DECREF(target);
    }

    Py_DECREF(target_iter);
    if (PyErr_Occurred()) {
        Py_DECREF(answer);
        return NULL;
    }
    return answer;
}


static PyObject *
evaluate_source_mod_3_covered(PyObject *self, PyObject *args)
{
    PyObject *target_faces;
    PyObject *support_object;
    PyObject *coefficients;
    Py_ssize_t target_degree;

    if (!PyArg_ParseTuple(
            args,
            "OOOn",
            &target_faces,
            &support_object,
            &coefficients,
            &target_degree
        )) {
        return NULL;
    }
    if (!PyDict_Check(coefficients)) {
        PyErr_SetString(PyExc_TypeError, "coefficients must be a dictionary");
        return NULL;
    }

    PyObject *support = PySequence_Fast(support_object, "support must be iterable");
    if (support == NULL) {
        return NULL;
    }

    Py_ssize_t support_length = PySequence_Fast_GET_SIZE(support);
    PyObject **support_items = PySequence_Fast_ITEMS(support);
    for (Py_ssize_t i = 0; i < support_length; i++) {
        PyObject *pair = support_items[i];
        if (!PyTuple_Check(pair) || PyTuple_GET_SIZE(pair) != 2) {
            Py_DECREF(support);
            PyErr_SetString(PyExc_TypeError, "support entries must be pairs");
            return NULL;
        }
        if (!PyTuple_Check(PyTuple_GET_ITEM(pair, 0))) {
            Py_DECREF(support);
            PyErr_SetString(PyExc_TypeError, "support faces must be tuples");
            return NULL;
        }
    }

    PyObject *answer = PyDict_New();
    if (answer == NULL) {
        Py_DECREF(support);
        return NULL;
    }

    Py_ssize_t target_length = target_degree + 1;
    PyObject *pattern;
    PyObject *pattern_coefficient_object;
    Py_ssize_t pattern_position = 0;
    while (PyDict_Next(coefficients, &pattern_position, &pattern, &pattern_coefficient_object)) {
        long pattern_coefficient = PyLong_AsLong(pattern_coefficient_object);
        if (pattern_coefficient == -1 && PyErr_Occurred()) {
            Py_DECREF(answer);
            Py_DECREF(support);
            return NULL;
        }

        PyObject *selected = selected_positions_for_pattern(pattern, target_length);
        if (selected == NULL) {
            Py_DECREF(answer);
            Py_DECREF(support);
            return NULL;
        }

        Py_ssize_t anchor_a = 0;
        Py_ssize_t anchor_b = 1;
        Py_ssize_t remaining = 2;
        if (best_anchor_pair(selected, &anchor_a, &anchor_b, &remaining) < 0) {
            Py_DECREF(selected);
            Py_DECREF(answer);
            Py_DECREF(support);
            return NULL;
        }

        PyObject *positions_a = PyTuple_GET_ITEM(selected, anchor_a);
        PyObject *positions_b = PyTuple_GET_ITEM(selected, anchor_b);
        PyObject *positions_c = PyTuple_GET_ITEM(selected, remaining);
        PyObject *fixed = remaining_fixed_data(positions_c, positions_a, positions_b);
        if (fixed == NULL) {
            Py_DECREF(selected);
            Py_DECREF(answer);
            Py_DECREF(support);
            return NULL;
        }
        PyObject *fixed_positions_c = PyTuple_GET_ITEM(fixed, 0);
        PyObject *fixed_indices_c = PyTuple_GET_ITEM(fixed, 1);
        PyObject *source_c_by_key = build_source_index(support, fixed_indices_c);
        if (source_c_by_key == NULL) {
            Py_DECREF(fixed);
            Py_DECREF(selected);
            Py_DECREF(answer);
            Py_DECREF(support);
            return NULL;
        }

        for (Py_ssize_t i = 0; i < support_length; i++) {
            PyObject *pair_a = support_items[i];
            PyObject *source_a = PyTuple_GET_ITEM(pair_a, 0);
            long coefficient_a = PyLong_AsLong(PyTuple_GET_ITEM(pair_a, 1));
            if (coefficient_a == -1 && PyErr_Occurred()) {
                Py_DECREF(source_c_by_key);
                Py_DECREF(fixed);
                Py_DECREF(selected);
                Py_DECREF(answer);
                Py_DECREF(support);
                return NULL;
            }

            for (Py_ssize_t j = 0; j < support_length; j++) {
                PyObject *pair_b = support_items[j];
                PyObject *source_b = PyTuple_GET_ITEM(pair_b, 0);
                long coefficient_b = PyLong_AsLong(PyTuple_GET_ITEM(pair_b, 1));
                if (coefficient_b == -1 && PyErr_Occurred()) {
                    Py_DECREF(source_c_by_key);
                    Py_DECREF(fixed);
                    Py_DECREF(selected);
                    Py_DECREF(answer);
                    Py_DECREF(support);
                    return NULL;
                }

                PyObject *partial_target = partial_target_from_pair(
                    source_a,
                    source_b,
                    positions_a,
                    positions_b,
                    target_length
                );
                if (partial_target == NULL) {
                    Py_DECREF(source_c_by_key);
                    Py_DECREF(fixed);
                    Py_DECREF(selected);
                    Py_DECREF(answer);
                    Py_DECREF(support);
                    return NULL;
                }
                if (partial_target == Py_None) {
                    Py_DECREF(partial_target);
                    continue;
                }

                PyObject *key = key_from_partial_target(partial_target, fixed_positions_c);
                if (key == NULL) {
                    Py_DECREF(partial_target);
                    Py_DECREF(source_c_by_key);
                    Py_DECREF(fixed);
                    Py_DECREF(selected);
                    Py_DECREF(answer);
                    Py_DECREF(support);
                    return NULL;
                }
                PyObject *source_c_bucket = PyDict_GetItemWithError(source_c_by_key, key);
                Py_DECREF(key);
                if (source_c_bucket == NULL && PyErr_Occurred()) {
                    Py_DECREF(partial_target);
                    Py_DECREF(source_c_by_key);
                    Py_DECREF(fixed);
                    Py_DECREF(selected);
                    Py_DECREF(answer);
                    Py_DECREF(support);
                    return NULL;
                }
                if (source_c_bucket == NULL) {
                    Py_DECREF(partial_target);
                    continue;
                }

                Py_ssize_t bucket_length = PyList_GET_SIZE(source_c_bucket);
                for (Py_ssize_t k = 0; k < bucket_length; k++) {
                    PyObject *pair_c = PyList_GET_ITEM(source_c_bucket, k);
                    PyObject *source_c = PyTuple_GET_ITEM(pair_c, 0);
                    long coefficient_c = PyLong_AsLong(PyTuple_GET_ITEM(pair_c, 1));
                    if (coefficient_c == -1 && PyErr_Occurred()) {
                        Py_DECREF(partial_target);
                        Py_DECREF(source_c_by_key);
                        Py_DECREF(fixed);
                        Py_DECREF(selected);
                        Py_DECREF(answer);
                        Py_DECREF(support);
                        return NULL;
                    }

                    PyObject *target = target_from_remaining_source(
                        partial_target,
                        source_c,
                        positions_c
                    );
                    if (target == NULL) {
                        Py_DECREF(partial_target);
                        Py_DECREF(source_c_by_key);
                        Py_DECREF(fixed);
                        Py_DECREF(selected);
                        Py_DECREF(answer);
                        Py_DECREF(support);
                        return NULL;
                    }
                    if (target == Py_None) {
                        Py_DECREF(target);
                        continue;
                    }

                    int contains = PySet_Contains(target_faces, target);
                    if (contains < 0) {
                        Py_DECREF(target);
                        Py_DECREF(partial_target);
                        Py_DECREF(source_c_by_key);
                        Py_DECREF(fixed);
                        Py_DECREF(selected);
                        Py_DECREF(answer);
                        Py_DECREF(support);
                        return NULL;
                    }
                    if (contains) {
                        long term_value = (
                            pattern_coefficient * coefficient_a * coefficient_b * coefficient_c
                        ) % 3;
                        if (term_value < 0) {
                            term_value += 3;
                        }
                        if (term_value && add_mod(answer, target, term_value, 3) < 0) {
                            Py_DECREF(target);
                            Py_DECREF(partial_target);
                            Py_DECREF(source_c_by_key);
                            Py_DECREF(fixed);
                            Py_DECREF(selected);
                            Py_DECREF(answer);
                            Py_DECREF(support);
                            return NULL;
                        }
                    }
                    Py_DECREF(target);
                }

                Py_DECREF(partial_target);
            }
        }

        Py_DECREF(source_c_by_key);
        Py_DECREF(fixed);
        Py_DECREF(selected);
    }

    Py_DECREF(support);
    return answer;
}


static PyMethodDef NativeMethods[] = {
    {
        "evaluate_all_targets",
        evaluate_all_targets,
        METH_VARARGS,
        "Evaluate an odd-primary universal tensor formula on all target faces.",
    },
    {
        "evaluate_source_mod_3_covered",
        evaluate_source_mod_3_covered,
        METH_VARARGS,
        "Evaluate covered p=3 omission patterns from source support.",
    },
    {NULL, NULL, 0, NULL},
};


static struct PyModuleDef native_module = {
    PyModuleDef_HEAD_INIT,
    "_native",
    "Native kernels for fastop.",
    -1,
    NativeMethods,
};


PyMODINIT_FUNC
PyInit__native(void)
{
    return PyModule_Create(&native_module);
}
