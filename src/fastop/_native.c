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


static PyMethodDef NativeMethods[] = {
    {
        "evaluate_all_targets",
        evaluate_all_targets,
        METH_VARARGS,
        "Evaluate an odd-primary universal tensor formula on all target faces.",
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
