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


static long
positive_mod(long value, long p)
{
    value %= p;
    if (value < 0) {
        value += p;
    }
    return value;
}


static long
mod_inverse(long value, long p)
{
    long t = 0;
    long new_t = 1;
    long r = p;
    long new_r = positive_mod(value, p);

    while (new_r != 0) {
        long quotient = r / new_r;
        long old_t = t;
        t = new_t;
        new_t = old_t - quotient * new_t;
        long old_r = r;
        r = new_r;
        new_r = old_r - quotient * new_r;
    }

    if (r > 1) {
        PyErr_SetString(PyExc_ValueError, "value is not invertible modulo p");
        return -1;
    }
    return positive_mod(t, p);
}


static PyObject *
clean_vector(PyObject *vector, long p)
{
    if (!PyDict_Check(vector)) {
        PyErr_SetString(PyExc_TypeError, "vectors must be dictionaries");
        return NULL;
    }

    PyObject *clean = PyDict_New();
    if (clean == NULL) {
        return NULL;
    }

    PyObject *key;
    PyObject *value_object;
    Py_ssize_t position = 0;
    while (PyDict_Next(vector, &position, &key, &value_object)) {
        long value = PyLong_AsLong(value_object);
        if (value == -1 && PyErr_Occurred()) {
            Py_DECREF(clean);
            return NULL;
        }
        value = positive_mod(value, p);
        if (value) {
            PyObject *new_value = PyLong_FromLong(value);
            if (new_value == NULL) {
                Py_DECREF(clean);
                return NULL;
            }
            int status = PyDict_SetItem(clean, key, new_value);
            Py_DECREF(new_value);
            if (status < 0) {
                Py_DECREF(clean);
                return NULL;
            }
        }
    }

    return clean;
}


static int
vector_add_inplace(PyObject *left, PyObject *right, long p, long scale)
{
    scale = positive_mod(scale, p);
    if (!scale) {
        return 0;
    }

    PyObject *key;
    PyObject *coefficient_object;
    Py_ssize_t position = 0;
    while (PyDict_Next(right, &position, &key, &coefficient_object)) {
        long coefficient = PyLong_AsLong(coefficient_object);
        if (coefficient == -1 && PyErr_Occurred()) {
            return -1;
        }
        if (add_mod(left, key, scale * coefficient, p) < 0) {
            return -1;
        }
    }

    return 0;
}


static int
vector_scale_inplace(PyObject *vector, long p, long scale)
{
    scale = positive_mod(scale, p);
    PyObject *items = PyMapping_Items(vector);
    if (items == NULL) {
        return -1;
    }
    PyObject *items_fast = PySequence_Fast(items, "vector items must be iterable");
    Py_DECREF(items);
    if (items_fast == NULL) {
        return -1;
    }

    Py_ssize_t length = PySequence_Fast_GET_SIZE(items_fast);
    PyObject **entries = PySequence_Fast_ITEMS(items_fast);
    for (Py_ssize_t i = 0; i < length; i++) {
        PyObject *entry = entries[i];
        PyObject *key = PyTuple_GET_ITEM(entry, 0);
        PyObject *value_object = PyTuple_GET_ITEM(entry, 1);
        long value = PyLong_AsLong(value_object);
        if (value == -1 && PyErr_Occurred()) {
            Py_DECREF(items_fast);
            return -1;
        }
        value = positive_mod(value * scale, p);
        if (value) {
            PyObject *new_value = PyLong_FromLong(value);
            if (new_value == NULL) {
                Py_DECREF(items_fast);
                return -1;
            }
            int status = PyDict_SetItem(vector, key, new_value);
            Py_DECREF(new_value);
            if (status < 0) {
                Py_DECREF(items_fast);
                return -1;
            }
        } else if (PyDict_DelItem(vector, key) < 0) {
            Py_DECREF(items_fast);
            return -1;
        }
    }

    Py_DECREF(items_fast);
    return 0;
}


static int
leading_index(PyObject *vector, Py_ssize_t *answer)
{
    PyObject *key;
    PyObject *value;
    Py_ssize_t position = 0;
    int seen = 0;
    Py_ssize_t largest = 0;

    while (PyDict_Next(vector, &position, &key, &value)) {
        Py_ssize_t index = PyLong_AsSsize_t(key);
        if (index == -1 && PyErr_Occurred()) {
            return -1;
        }
        if (!seen || index > largest) {
            largest = index;
            seen = 1;
        }
    }

    if (!seen) {
        PyErr_SetString(PyExc_ValueError, "empty vector has no leading index");
        return -1;
    }
    *answer = largest;
    return 0;
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
column_image_and_kernel_basis(PyObject *self, PyObject *args)
{
    PyObject *columns_object;
    long p;

    if (!PyArg_ParseTuple(args, "Ol", &columns_object, &p)) {
        return NULL;
    }
    if (p <= 1) {
        PyErr_SetString(PyExc_ValueError, "p must be greater than one");
        return NULL;
    }

    PyObject *columns = PySequence_Fast(columns_object, "columns must be iterable");
    if (columns == NULL) {
        return NULL;
    }

    PyObject *reduced_columns = PyList_New(0);
    PyObject *transforms = PyList_New(0);
    PyObject *pivot_to_column = PyDict_New();
    PyObject *cycles = PyList_New(0);
    if (
        reduced_columns == NULL
        || transforms == NULL
        || pivot_to_column == NULL
        || cycles == NULL
    ) {
        Py_XDECREF(reduced_columns);
        Py_XDECREF(transforms);
        Py_XDECREF(pivot_to_column);
        Py_XDECREF(cycles);
        Py_DECREF(columns);
        return NULL;
    }

    Py_ssize_t column_count = PySequence_Fast_GET_SIZE(columns);
    PyObject **column_items = PySequence_Fast_ITEMS(columns);
    for (Py_ssize_t column_index = 0; column_index < column_count; column_index++) {
        PyObject *reduced = clean_vector(column_items[column_index], p);
        if (reduced == NULL) {
            Py_DECREF(reduced_columns);
            Py_DECREF(transforms);
            Py_DECREF(pivot_to_column);
            Py_DECREF(cycles);
            Py_DECREF(columns);
            return NULL;
        }

        PyObject *transform = PyDict_New();
        PyObject *column_index_object = PyLong_FromSsize_t(column_index);
        PyObject *one = PyLong_FromLong(1);
        if (transform == NULL || column_index_object == NULL || one == NULL) {
            Py_XDECREF(transform);
            Py_XDECREF(column_index_object);
            Py_XDECREF(one);
            Py_DECREF(reduced);
            Py_DECREF(reduced_columns);
            Py_DECREF(transforms);
            Py_DECREF(pivot_to_column);
            Py_DECREF(cycles);
            Py_DECREF(columns);
            return NULL;
        }
        if (PyDict_SetItem(transform, column_index_object, one) < 0) {
            Py_DECREF(transform);
            Py_DECREF(column_index_object);
            Py_DECREF(one);
            Py_DECREF(reduced);
            Py_DECREF(reduced_columns);
            Py_DECREF(transforms);
            Py_DECREF(pivot_to_column);
            Py_DECREF(cycles);
            Py_DECREF(columns);
            return NULL;
        }
        Py_DECREF(column_index_object);
        Py_DECREF(one);

        while (PyDict_Size(reduced)) {
            Py_ssize_t pivot;
            if (leading_index(reduced, &pivot) < 0) {
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            PyObject *pivot_object = PyLong_FromSsize_t(pivot);
            if (pivot_object == NULL) {
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            PyObject *pivot_column_index_object = PyDict_GetItemWithError(
                pivot_to_column,
                pivot_object
            );
            if (pivot_column_index_object == NULL && PyErr_Occurred()) {
                Py_DECREF(pivot_object);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            if (pivot_column_index_object == NULL) {
                Py_DECREF(pivot_object);
                break;
            }

            Py_ssize_t pivot_column_index = PyLong_AsSsize_t(pivot_column_index_object);
            if (pivot_column_index == -1 && PyErr_Occurred()) {
                Py_DECREF(pivot_object);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            PyObject *pivot_column = PyList_GET_ITEM(reduced_columns, pivot_column_index);
            PyObject *pivot_transform = PyList_GET_ITEM(transforms, pivot_column_index);
            PyObject *reduced_coefficient_object = PyDict_GetItemWithError(
                reduced,
                pivot_object
            );
            PyObject *pivot_coefficient_object = PyDict_GetItemWithError(
                pivot_column,
                pivot_object
            );
            if (
                (reduced_coefficient_object == NULL && PyErr_Occurred())
                || (pivot_coefficient_object == NULL && PyErr_Occurred())
                || reduced_coefficient_object == NULL
                || pivot_coefficient_object == NULL
            ) {
                if (!PyErr_Occurred()) {
                    PyErr_SetString(PyExc_RuntimeError, "missing pivot coefficient");
                }
                Py_DECREF(pivot_object);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }

            long reduced_coefficient = PyLong_AsLong(reduced_coefficient_object);
            long pivot_coefficient = PyLong_AsLong(pivot_coefficient_object);
            if (
                (reduced_coefficient == -1 && PyErr_Occurred())
                || (pivot_coefficient == -1 && PyErr_Occurred())
            ) {
                Py_DECREF(pivot_object);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            long inverse = mod_inverse(pivot_coefficient, p);
            if (inverse == -1 && PyErr_Occurred()) {
                Py_DECREF(pivot_object);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            long coefficient = positive_mod(reduced_coefficient * inverse, p);
            if (
                vector_add_inplace(reduced, pivot_column, p, -coefficient) < 0
                || vector_add_inplace(transform, pivot_transform, p, -coefficient) < 0
            ) {
                Py_DECREF(pivot_object);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            Py_DECREF(pivot_object);
        }

        if (PyDict_Size(reduced)) {
            Py_ssize_t pivot;
            if (leading_index(reduced, &pivot) < 0) {
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            PyObject *pivot_object = PyLong_FromSsize_t(pivot);
            PyObject *column_position = PyLong_FromSsize_t(PyList_GET_SIZE(reduced_columns));
            if (pivot_object == NULL || column_position == NULL) {
                Py_XDECREF(pivot_object);
                Py_XDECREF(column_position);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }

            PyObject *pivot_coefficient_object = PyDict_GetItemWithError(
                reduced,
                pivot_object
            );
            if (pivot_coefficient_object == NULL) {
                if (!PyErr_Occurred()) {
                    PyErr_SetString(PyExc_RuntimeError, "missing pivot coefficient");
                }
                Py_DECREF(pivot_object);
                Py_DECREF(column_position);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            long pivot_coefficient = PyLong_AsLong(pivot_coefficient_object);
            if (pivot_coefficient == -1 && PyErr_Occurred()) {
                Py_DECREF(pivot_object);
                Py_DECREF(column_position);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            long inverse = mod_inverse(pivot_coefficient, p);
            if (inverse == -1 && PyErr_Occurred()) {
                Py_DECREF(pivot_object);
                Py_DECREF(column_position);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            if (
                vector_scale_inplace(reduced, p, inverse) < 0
                || vector_scale_inplace(transform, p, inverse) < 0
            ) {
                Py_DECREF(pivot_object);
                Py_DECREF(column_position);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            if (
                PyDict_SetItem(pivot_to_column, pivot_object, column_position) < 0
                || PyList_Append(reduced_columns, reduced) < 0
                || PyList_Append(transforms, transform) < 0
            ) {
                Py_DECREF(pivot_object);
                Py_DECREF(column_position);
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            Py_DECREF(pivot_object);
            Py_DECREF(column_position);
            Py_DECREF(transform);
            Py_DECREF(reduced);
        } else {
            if (PyList_Append(cycles, transform) < 0) {
                Py_DECREF(transform);
                Py_DECREF(reduced);
                Py_DECREF(reduced_columns);
                Py_DECREF(transforms);
                Py_DECREF(pivot_to_column);
                Py_DECREF(cycles);
                Py_DECREF(columns);
                return NULL;
            }
            Py_DECREF(transform);
            Py_DECREF(reduced);
        }
    }

    PyObject *answer = PyTuple_New(2);
    if (answer == NULL) {
        Py_DECREF(reduced_columns);
        Py_DECREF(transforms);
        Py_DECREF(pivot_to_column);
        Py_DECREF(cycles);
        Py_DECREF(columns);
        return NULL;
    }
    PyTuple_SET_ITEM(answer, 0, reduced_columns);
    PyTuple_SET_ITEM(answer, 1, cycles);
    Py_DECREF(transforms);
    Py_DECREF(pivot_to_column);
    Py_DECREF(columns);
    return answer;
}


static int
max_reducible_pivot(PyObject *vector, PyObject *rows, Py_ssize_t *pivot)
{
    PyObject *key;
    PyObject *value;
    Py_ssize_t position = 0;
    int seen = 0;
    Py_ssize_t largest = 0;

    while (PyDict_Next(vector, &position, &key, &value)) {
        PyObject *row_entry = PyDict_GetItemWithError(rows, key);
        if (row_entry == NULL && PyErr_Occurred()) {
            return -1;
        }
        if (row_entry == NULL) {
            continue;
        }
        Py_ssize_t index = PyLong_AsSsize_t(key);
        if (index == -1 && PyErr_Occurred()) {
            return -1;
        }
        if (!seen || index > largest) {
            largest = index;
            seen = 1;
        }
    }

    if (!seen) {
        return 0;
    }
    *pivot = largest;
    return 1;
}


static int
reduce_with_rows(PyObject *vector, PyObject *coordinate, PyObject *rows, long p)
{
    while (1) {
        Py_ssize_t pivot;
        int found = max_reducible_pivot(vector, rows, &pivot);
        if (found < 0) {
            return -1;
        }
        if (!found) {
            return 0;
        }

        PyObject *pivot_object = PyLong_FromSsize_t(pivot);
        if (pivot_object == NULL) {
            return -1;
        }
        PyObject *coefficient_object = PyDict_GetItemWithError(vector, pivot_object);
        PyObject *row_entry = PyDict_GetItemWithError(rows, pivot_object);
        Py_DECREF(pivot_object);
        if (
            (coefficient_object == NULL && PyErr_Occurred())
            || (row_entry == NULL && PyErr_Occurred())
            || coefficient_object == NULL
            || row_entry == NULL
        ) {
            if (!PyErr_Occurred()) {
                PyErr_SetString(PyExc_RuntimeError, "missing reducer pivot");
            }
            return -1;
        }

        long coefficient = PyLong_AsLong(coefficient_object);
        if (coefficient == -1 && PyErr_Occurred()) {
            return -1;
        }
        PyObject *row = PyTuple_GET_ITEM(row_entry, 0);
        PyObject *row_coordinate = PyTuple_GET_ITEM(row_entry, 1);
        if (vector_add_inplace(vector, row, p, -coefficient) < 0) {
            return -1;
        }
        if (
            coordinate != NULL
            && vector_add_inplace(coordinate, row_coordinate, p, -coefficient) < 0
        ) {
            return -1;
        }
    }
}


static int
add_coordinate_row(PyObject *rows, PyObject *vector_object, PyObject *coordinate_object, long p)
{
    PyObject *vector = clean_vector(vector_object, p);
    PyObject *coordinate;
    if (coordinate_object == NULL) {
        coordinate = PyDict_New();
    } else {
        coordinate = clean_vector(coordinate_object, p);
    }
    if (vector == NULL || coordinate == NULL) {
        Py_XDECREF(vector);
        Py_XDECREF(coordinate);
        return -1;
    }

    if (reduce_with_rows(vector, coordinate, rows, p) < 0) {
        Py_DECREF(vector);
        Py_DECREF(coordinate);
        return -1;
    }
    if (PyDict_Size(vector) == 0) {
        Py_DECREF(vector);
        Py_DECREF(coordinate);
        return 0;
    }

    Py_ssize_t pivot;
    if (leading_index(vector, &pivot) < 0) {
        Py_DECREF(vector);
        Py_DECREF(coordinate);
        return -1;
    }
    PyObject *pivot_object = PyLong_FromSsize_t(pivot);
    if (pivot_object == NULL) {
        Py_DECREF(vector);
        Py_DECREF(coordinate);
        return -1;
    }
    PyObject *pivot_coefficient_object = PyDict_GetItemWithError(vector, pivot_object);
    if (pivot_coefficient_object == NULL) {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "missing coordinate pivot coefficient");
        }
        Py_DECREF(pivot_object);
        Py_DECREF(vector);
        Py_DECREF(coordinate);
        return -1;
    }
    long pivot_coefficient = PyLong_AsLong(pivot_coefficient_object);
    if (pivot_coefficient == -1 && PyErr_Occurred()) {
        Py_DECREF(pivot_object);
        Py_DECREF(vector);
        Py_DECREF(coordinate);
        return -1;
    }
    long inverse = mod_inverse(pivot_coefficient, p);
    if (inverse == -1 && PyErr_Occurred()) {
        Py_DECREF(pivot_object);
        Py_DECREF(vector);
        Py_DECREF(coordinate);
        return -1;
    }
    if (
        vector_scale_inplace(vector, p, inverse) < 0
        || vector_scale_inplace(coordinate, p, inverse) < 0
    ) {
        Py_DECREF(pivot_object);
        Py_DECREF(vector);
        Py_DECREF(coordinate);
        return -1;
    }

    PyObject *row_entry = PyTuple_New(2);
    if (row_entry == NULL) {
        Py_DECREF(pivot_object);
        Py_DECREF(vector);
        Py_DECREF(coordinate);
        return -1;
    }
    PyTuple_SET_ITEM(row_entry, 0, vector);
    PyTuple_SET_ITEM(row_entry, 1, coordinate);
    int status = PyDict_SetItem(rows, pivot_object, row_entry);
    Py_DECREF(pivot_object);
    Py_DECREF(row_entry);
    if (status < 0) {
        return -1;
    }
    return 1;
}


static PyObject *
coordinate_basis_from_vectors(PyObject *self, PyObject *args)
{
    PyObject *boundary_vectors_object;
    PyObject *cycles_object;
    long p;

    if (!PyArg_ParseTuple(args, "OOl", &boundary_vectors_object, &cycles_object, &p)) {
        return NULL;
    }
    PyObject *boundary_vectors = PySequence_Fast(
        boundary_vectors_object,
        "boundary vectors must be iterable"
    );
    PyObject *cycles = PySequence_Fast(cycles_object, "cycles must be iterable");
    if (boundary_vectors == NULL || cycles == NULL) {
        Py_XDECREF(boundary_vectors);
        Py_XDECREF(cycles);
        return NULL;
    }

    PyObject *rows = PyDict_New();
    PyObject *cocycle_basis = PyList_New(0);
    if (rows == NULL || cocycle_basis == NULL) {
        Py_XDECREF(rows);
        Py_XDECREF(cocycle_basis);
        Py_DECREF(boundary_vectors);
        Py_DECREF(cycles);
        return NULL;
    }

    Py_ssize_t boundary_count = PySequence_Fast_GET_SIZE(boundary_vectors);
    PyObject **boundary_items = PySequence_Fast_ITEMS(boundary_vectors);
    for (Py_ssize_t i = 0; i < boundary_count; i++) {
        if (add_coordinate_row(rows, boundary_items[i], NULL, p) < 0) {
            Py_DECREF(rows);
            Py_DECREF(cocycle_basis);
            Py_DECREF(boundary_vectors);
            Py_DECREF(cycles);
            return NULL;
        }
    }

    Py_ssize_t cycle_count = PySequence_Fast_GET_SIZE(cycles);
    PyObject **cycle_items = PySequence_Fast_ITEMS(cycles);
    for (Py_ssize_t i = 0; i < cycle_count; i++) {
        PyObject *coordinate = PyDict_New();
        PyObject *coordinate_index = PyLong_FromSsize_t(PyList_GET_SIZE(cocycle_basis));
        PyObject *one = PyLong_FromLong(1);
        if (coordinate == NULL || coordinate_index == NULL || one == NULL) {
            Py_XDECREF(coordinate);
            Py_XDECREF(coordinate_index);
            Py_XDECREF(one);
            Py_DECREF(rows);
            Py_DECREF(cocycle_basis);
            Py_DECREF(boundary_vectors);
            Py_DECREF(cycles);
            return NULL;
        }
        if (PyDict_SetItem(coordinate, coordinate_index, one) < 0) {
            Py_DECREF(coordinate);
            Py_DECREF(coordinate_index);
            Py_DECREF(one);
            Py_DECREF(rows);
            Py_DECREF(cocycle_basis);
            Py_DECREF(boundary_vectors);
            Py_DECREF(cycles);
            return NULL;
        }
        Py_DECREF(coordinate_index);
        Py_DECREF(one);

        int added = add_coordinate_row(rows, cycle_items[i], coordinate, p);
        Py_DECREF(coordinate);
        if (added < 0) {
            Py_DECREF(rows);
            Py_DECREF(cocycle_basis);
            Py_DECREF(boundary_vectors);
            Py_DECREF(cycles);
            return NULL;
        }
        if (added) {
            if (PyList_Append(cocycle_basis, cycle_items[i]) < 0) {
                Py_DECREF(rows);
                Py_DECREF(cocycle_basis);
                Py_DECREF(boundary_vectors);
                Py_DECREF(cycles);
                return NULL;
            }
        }
    }

    PyObject *answer = PyTuple_New(2);
    if (answer == NULL) {
        Py_DECREF(rows);
        Py_DECREF(cocycle_basis);
        Py_DECREF(boundary_vectors);
        Py_DECREF(cycles);
        return NULL;
    }
    PyTuple_SET_ITEM(answer, 0, cocycle_basis);
    PyTuple_SET_ITEM(answer, 1, rows);
    Py_DECREF(boundary_vectors);
    Py_DECREF(cycles);
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
    {
        "column_image_and_kernel_basis",
        column_image_and_kernel_basis,
        METH_VARARGS,
        "Return image and kernel bases using TDA-style column reduction.",
    },
    {
        "coordinate_basis_from_vectors",
        coordinate_basis_from_vectors,
        METH_VARARGS,
        "Build cocycle representatives and coordinate pivot rows.",
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
