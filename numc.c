#include "numc.h"
#include <structmember.h>

PyTypeObject Matrix61cType;

/* Helper functions for initalization of matrices and vectors */

/*
 * Return a tuple given rows and cols
 */
PyObject *get_shape(int rows, int cols) {
  if (rows == 1 || cols == 1) {
    return PyTuple_Pack(1, PyLong_FromLong(rows * cols));
  } else {
    return PyTuple_Pack(2, PyLong_FromLong(rows), PyLong_FromLong(cols));
  }
}
/*
 * Matrix(rows, cols, low, high). Fill a matrix random double values
 */
int init_rand(PyObject *self, int rows, int cols, unsigned int seed, double low,
              double high) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed != 0) {
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate matrix.");
        return alloc_failed;
    }
    rand_matrix(new_mat, seed, low, high);
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(rows, cols, val). Fill a matrix of dimension rows * cols with val
 */
int init_fill(PyObject *self, int rows, int cols, double val) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed != 0) {
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate matrix.");
        return alloc_failed;
    }
    else {
        fill_matrix(new_mat, val);
        ((Matrix61c *)self)->mat = new_mat;
        ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    }
    return 0;
}

/*
 * Matrix(rows, cols, 1d_list). Fill a matrix with dimension rows * cols with 1d_list values
 */
int init_1d(PyObject *self, int rows, int cols, PyObject *lst) {
    if (rows * cols != PyList_Size(lst)) {
        PyErr_SetString(PyExc_ValueError, "Incorrect number of elements in list");
        return -1;
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed != 0) {
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate matrix.");
        return alloc_failed;
    }
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j, PyFloat_AsDouble(PyList_GetItem(lst, count)));
            count++;
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(2d_list). Fill a matrix with dimension len(2d_list) * len(2d_list[0])
 */
int init_2d(PyObject *self, PyObject *lst) {
    int rows = PyList_Size(lst);
    if (rows == 0) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot initialize numc.Matrix with an empty list");
        return -1;
    }
    int cols;
    if (!PyList_Check(PyList_GetItem(lst, 0))) {
        PyErr_SetString(PyExc_ValueError, "List values not valid");
        return -1;
    } else {
        cols = PyList_Size(PyList_GetItem(lst, 0));
    }
    for (int i = 0; i < rows; i++) {
        if (!PyList_Check(PyList_GetItem(lst, i)) ||
                PyList_Size(PyList_GetItem(lst, i)) != cols) {
            PyErr_SetString(PyExc_ValueError, "List values not valid");
            return -1;
        }
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed != 0) {
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate matrix.");
        return alloc_failed;
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j,
                PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(lst, i), j)));
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * This deallocation function is called when reference count is 0
 */
void Matrix61c_dealloc(Matrix61c *self) {
    deallocate_matrix(self->mat);
    Py_TYPE(self)->tp_free(self);
}

/* For immutable types all initializations should take place in tp_new */
PyObject *Matrix61c_new(PyTypeObject *type, PyObject *args,
                        PyObject *kwds) {
    /* size of allocated memory is tp_basicsize + nitems*tp_itemsize*/
    Matrix61c *self = (Matrix61c *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

/*
 * This matrix61c type is mutable, so needs init function. Return 0 on success otherwise -1
 */
int Matrix61c_init(PyObject *self, PyObject *args, PyObject *kwds) {
    /* Generate random matrices */
    if (kwds != NULL) {
        PyObject *rand = PyDict_GetItemString(kwds, "rand");
        if (!rand) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (!PyBool_Check(rand)) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (rand != Py_True) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        PyObject *low = PyDict_GetItemString(kwds, "low");
        PyObject *high = PyDict_GetItemString(kwds, "high");
        PyObject *seed = PyDict_GetItemString(kwds, "seed");
        double double_low = 0;
        double double_high = 1;
        unsigned int unsigned_seed = 0;

        if (low) {
            if (PyFloat_Check(low)) {
                double_low = PyFloat_AsDouble(low);
            } else if (PyLong_Check(low)) {
                double_low = PyLong_AsLong(low);
            }
        }

        if (high) {
            if (PyFloat_Check(high)) {
                double_high = PyFloat_AsDouble(high);
            } else if (PyLong_Check(high)) {
                double_high = PyLong_AsLong(high);
            }
        }

        if (double_low >= double_high) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        // Set seed if argument exists
        if (seed) {
            if (PyLong_Check(seed)) {
                unsigned_seed = PyLong_AsUnsignedLong(seed);
            }
        }

        PyObject *rows = NULL;
        PyObject *cols = NULL;
        if (PyArg_UnpackTuple(args, "args", 2, 2, &rows, &cols)) {
            if (rows && cols && PyLong_Check(rows) && PyLong_Check(cols)) {
                return init_rand(self, PyLong_AsLong(rows), PyLong_AsLong(cols), unsigned_seed, double_low,
                                 double_high);
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    }
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 3, &arg1, &arg2, &arg3)) {
        /* arguments are (rows, cols, val) */
        if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && (PyLong_Check(arg3)
                || PyFloat_Check(arg3))) {
            if (PyLong_Check(arg3)) {
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyLong_AsLong(arg3));
            } else
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyFloat_AsDouble(arg3));
        } else if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && PyList_Check(arg3)) {
            /* Matrix(rows, cols, 1D list) */
            return init_1d(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), arg3);
        } else if (arg1 && PyList_Check(arg1) && arg2 == NULL && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_2d(self, arg1);
        } else if (arg1 && arg2 && PyLong_Check(arg1) && PyLong_Check(arg2) && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), 0);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return -1;
    }
}

/*
 * List of lists representations for matrices
 */
PyObject *Matrix61c_to_list(Matrix61c *self) {
    int rows = self->mat->rows;
    int cols = self->mat->cols;
    PyObject *py_lst = NULL;
    if (self->mat->is_1d) {  // If 1D matrix, print as a single list
        py_lst = PyList_New(rows * cols);
        int count = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(py_lst, count, PyFloat_FromDouble(get(self->mat, i, j)));
                count++;
            }
        }
    } else {  // if 2D, print as nested list
        py_lst = PyList_New(rows);
        for (int i = 0; i < rows; i++) {
            PyList_SetItem(py_lst, i, PyList_New(cols));
            PyObject *curr_row = PyList_GetItem(py_lst, i);
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(curr_row, j, PyFloat_FromDouble(get(self->mat, i, j)));
            }
        }
    }
    return py_lst;
}

PyObject *Matrix61c_class_to_list(Matrix61c *self, PyObject *args) {
    PyObject *mat = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 1, &mat)) {
        if (!PyObject_TypeCheck(mat, &Matrix61cType)) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
            return NULL;
        }
        Matrix61c* mat61c = (Matrix61c*)mat;
        return Matrix61c_to_list(mat61c);
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
}

/*
 * Add class methods
 */
PyMethodDef Matrix61c_class_methods[] = {
    {"to_list", (PyCFunction)Matrix61c_class_to_list, METH_VARARGS, "Returns a list representation of numc.Matrix"},
    {NULL, NULL, 0, NULL}
};

/*
 * Matrix61c string representation. For printing purposes.
 */
PyObject *Matrix61c_repr(PyObject *self) {
    PyObject *py_lst = Matrix61c_to_list((Matrix61c *)self);
    return PyObject_Repr(py_lst);
}

/* NUMBER METHODS */

/*
 * Add the second numc.Matrix (Matrix61c) object to the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_add(Matrix61c* self, PyObject* args) {
    /* TODO: test */
    if (!PyObject_TypeCheck(args,&Matrix61cType))
    {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }

    Matrix61c *second = (Matrix61c*) args;

    PyObject *self_r = NULL; //= PyTuple_GetItem(self->shape, 0);
    PyObject *self_c = NULL; //= PyTuple_GetItem(self->shape, 1);
    PyObject *secd_r = NULL; // = PyTuple_GetItem(second->shape, 0);
    PyObject *secd_c = NULL; //= PyTuple_GetItem(second->shape, 1);

    int parsed0 = PyArg_UnpackTuple(self->shape, "args", 2, 2, &self_r, &self_c);
    if (parsed0==-1) {
        // parse failed
        return NULL;
    }
    int parsed1 = PyArg_UnpackTuple(self->shape, "args", 2, 2, &secd_r, &secd_c);
    if (parsed1==-1) {
        // parse failed
        return NULL;
    }
    int r1 = PyLong_AsLong(self_r);
    int r2 = PyLong_AsLong(secd_r);
    int c1 = PyLong_AsLong(self_c);
    int c2 = PyLong_AsLong(secd_c);
    if (r1 == -1|| r2 == -1 || r1 != r2 || c1 == -1 || c2==-1 || c1 != c2){ // self->shape != second->shape { //may need to do individual comparison or build comparator?
        PyErr_SetString(PyExc_ValueError, "Shapes are not the same.");
        return NULL;
    }


    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    rv->shape = get_shape(r1, c1);
    matrix *new_mat;

    int alloc_failed = allocate_matrix(&new_mat, r1, c1);
    if (alloc_failed != 0) {
        // TODO: failed to allocate
        PyErr_SetString(PyExc_RuntimeError, "(Add): failed to allocate matrix.");
        return NULL;
    }
    rv->mat = new_mat;
    // init_fill((PyObject*)rv, self->mat->rows, self->mat->cols, 0);
    int pass = add_matrix(rv->mat, self->mat, second->mat);

    if (pass != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot perform matrix add.");
        return NULL;
    }
    PyObject *retval = (PyObject*) rv;
    // printf("%s\n", "HELLO WORKED PAST line 352, casting to pyobject.");
    return retval;
}

/*
 * Substract the second numc.Matrix (Matrix61c) object from the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_sub(Matrix61c* self, PyObject* args) {
    /* TODO: test */
    if (PyObject_TypeCheck(args,&Matrix61cType) == 0)
    {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }

    Matrix61c *second = (Matrix61c*) args;
    PyObject *self_r = PyTuple_GetItem(self->shape, 0);
    PyObject *self_c = PyTuple_GetItem(self->shape, 1);
    PyObject *secd_r = PyTuple_GetItem(second->shape, 0);
    PyObject *secd_c = PyTuple_GetItem(second->shape, 1);
    if (PyLong_AsLong(self_r) != PyLong_AsLong(secd_r) || PyLong_AsLong(self_c)!= PyLong_AsLong(secd_c)) {// self->shape != second->shape { //may need to do individual comparison or build comparator?
        PyErr_SetString(PyExc_ValueError, "Shapes are not the same.");
        return NULL;
    }
    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    rv->shape = self->shape;
    matrix *new_mat;
    int pass = allocate_matrix(&new_mat, PyLong_AsLong(self_r), PyLong_AsLong(self_c));
    if (pass != 0) {
        // TODO: failed to allocate
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate matrix.");
        return NULL;
    }
    rv->mat = new_mat;
    pass = sub_matrix(rv->mat, self->mat, second->mat);
    if (pass != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot perform matrix add.");
        return NULL;
    }
    return (PyObject *) rv;

}

/*
 * NOT element-wise multiplication. The first operand is self, and the second operand
 * can be obtained by casting `args`.
 */
PyObject *Matrix61c_multiply(Matrix61c* self, PyObject *args) {

    if (PyObject_TypeCheck(args,&Matrix61cType) == 0)
    {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }

    Matrix61c *second = (Matrix61c*) args;

    PyObject *self_r = PyTuple_GetItem(self->shape, 0);
    PyObject *secd_c = PyTuple_GetItem(second->shape, 1);
    PyObject *self_c = PyTuple_GetItem(self->shape, 1);
    PyObject *secd_r = PyTuple_GetItem(second->shape, 0);
    if (PyLong_AsLong(self_c) != PyLong_AsLong(secd_r)){ // self->shape != second->shape { //may need to do individual comparison or build comparator?
        PyErr_SetString(PyExc_ValueError, "Shapes are not the same.");
        return NULL;
    }
    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    rv->shape = get_shape(PyLong_AsLong(self_r), PyLong_AsLong(secd_c));
    matrix *new_mat;
    int pass = allocate_matrix(&new_mat, PyLong_AsLong(self_r), PyLong_AsLong(secd_c));
    if (pass != 0) {
        // TODO: failed to allocate
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate matrix.");
        return NULL;
    }
    rv->mat = new_mat;
    pass = mul_matrix(rv->mat, self->mat, second->mat);
    if (pass != 0) {
        PyErr_SetString(PyExc_ValueError, "Cannot perform matrix mul if # cols A != # cols B.");
        return NULL;
    }
    return (PyObject *) rv;
}

/*
 * Negates the given numc.Matrix.
 */
PyObject *Matrix61c_neg(Matrix61c* self) {
    /* TODO: test */
    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    rv->shape = self->shape;
    matrix *new_mat;
    int pass = allocate_matrix(&new_mat, self->mat->rows, self->mat->cols);
    if (pass != 0) {
        // TODO: failed to allocate
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate matrix.");
        return NULL;
    }
    rv->mat = new_mat;
    pass = neg_matrix(rv->mat, self->mat);
    if (pass != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot perform matrix neg.");
        return NULL;
    }
    return (PyObject *) rv;
}

/*
 * Take the element-wise absolute value of this numc.Matrix.
 */
PyObject *Matrix61c_abs(Matrix61c *self) {
    /* TODO: YOUR CODE HERE */
    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    rv->shape = self->shape;
    matrix *new_mat;
    int pass = allocate_matrix(&new_mat, self->mat->rows, self->mat->cols);
    if (pass != 0) {
        // TODO: failed to allocate
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate matrix.");
        return NULL;
    }
    rv->mat = new_mat;
    pass = abs_matrix(rv->mat, self->mat);
    if (pass != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot perform matrix abs.");
        return NULL;
    }
    return (PyObject *) rv;
}

/*
 * Raise numc.Matrix (Matrix61c) to the `pow`th power. You can ignore the argument `optional`.
 */
PyObject *Matrix61c_pow(Matrix61c *self, PyObject *pow, PyObject *optional) {
    if(!PyLong_Check(pow)){
        PyErr_SetString(PyExc_TypeError, "Invalid arguments.");
        return NULL;
    }
    int p = PyLong_AsLong(pow);
    if (p < 0) {
        PyErr_SetString(PyExc_ValueError, "Power cannot be negative.");
        return NULL;
    }
    PyObject *self_r = PyTuple_GetItem(self->shape, 0);
    PyObject *self_c = PyTuple_GetItem(self->shape, 1);
    if (PyLong_AsLong(self_r) != PyLong_AsLong(self_c)){
        PyErr_SetString(PyExc_ValueError, "Not a square matrix.");
        return NULL;
    }
    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    rv->shape = self->shape;
    matrix *new_mat;
    int pass = allocate_matrix(&new_mat, PyLong_AsLong(self_r), PyLong_AsLong(self_c));
    if (pass != 0) {
        // TODO: failed to allocate
        PyErr_SetString(PyExc_RuntimeError, "failed to allocate matrix.");
        return NULL;
    }
    rv->mat = new_mat;
    pass = pow_matrix(rv->mat, self->mat, p);
    if (pass != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot perform matrix pow.");
        return NULL;
    }
    return (PyObject *) rv;
}

/*
 * Create a PyNumberMethods struct for overloading operators with all the number methods you have
 * define. You might find this link helpful: https://docs.python.org/3.6/c-api/typeobj.html
 */
PyNumberMethods Matrix61c_as_number = {
    /* TODO: test*/
        .nb_add = (binaryfunc) Matrix61c_add,
        .nb_subtract = (binaryfunc) Matrix61c_sub,
        .nb_multiply = (binaryfunc) Matrix61c_multiply,
        .nb_remainder= 0,
        .nb_divmod = 0,
        .nb_power = (ternaryfunc) Matrix61c_pow,
        .nb_negative = (unaryfunc) Matrix61c_neg,
        .nb_positive = 0,
        .nb_absolute = (unaryfunc) Matrix61c_abs,
        .nb_bool= 0,
        .nb_invert= 0,
        .nb_lshift= 0,
        .nb_rshift= 0,
        .nb_and= 0,
        .nb_xor= 0,
        .nb_or= 0,
        .nb_int= 0,
        .nb_reserved= 0,
        .nb_float= 0,
        .nb_inplace_add= 0,
        .nb_inplace_subtract= 0,
        .nb_inplace_multiply= 0,
        .nb_inplace_remainder= 0,
        .nb_inplace_power= 0,
        .nb_inplace_lshift= 0,
        .nb_inplace_rshift= 0,
        .nb_inplace_and= 0,
        .nb_inplace_xor= 0,
        .nb_inplace_or= 0,
        .nb_floor_divide= 0,
        .nb_true_divide= 0,
        .nb_inplace_floor_divide= 0,
        .nb_inplace_true_divide= 0,
        .nb_index= 0,
        .nb_matrix_multiply = (binaryfunc) Matrix61c_multiply, // multiply
        .nb_inplace_matrix_multiply= 0
};


/* INSTANCE METHODS */

/*
 * Given a numc.Matrix self, parse `args` to (int) row, (int) col, and (double/int) val.
 * Return None in Python (this is different from returning null).
 */
PyObject *Matrix61c_set_value(Matrix61c *self, PyObject* args) {
    // parse args
    PyObject *row;
    PyObject *col;
    PyObject *value;
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    PyObject *rs;
    PyObject *cs;
    if (PyArg_UnpackTuple(args, "args", 3, 3, &arg1, &arg2, &arg3)) {
        if(PyTuple_Size(self->shape) == 2) {
            row = arg1;
            col = arg2;
            value = arg3;
            rs = PyTuple_GetItem(self->shape, 0);
            cs = PyTuple_GetItem(self->shape, 1);
        }
        else if (PyTuple_Size(self->shape) == 1) {
            row = arg1;
            if (!arg3) { // arg3 not inputed
                value = arg2;
                col = PyLong_FromLong(0);
            } else {
                col = arg2;
                value = arg3;
            }
            rs = PyTuple_GetItem(self->shape, 0);
            cs = PyLong_FromLong(1);
        }
        else {
            // throw error bc shape should be 2 only
            PyErr_SetString(PyExc_TypeError, "Invalid Matrix input.");
            return NULL;
        } // ERR
        // arguments parsed

        int rows = PyLong_AsLong(rs);
        int cols = PyLong_AsLong(cs);
        if (PyLong_AsLong(row) >= rows || PyLong_AsLong(col) >= cols || PyLong_AsLong(row) < 0 || PyLong_AsLong(col) < 0) {
            PyErr_SetString(PyExc_IndexError, "Index out of range.");
            return NULL;
        }
        if (row && col && value && PyLong_Check(row) && PyLong_Check(col) && PyLong_Check(value))
        {
            // int value
            set(self->mat, PyLong_AsLong(row), PyLong_AsLong(col), PyLong_AsDouble(value));
        } else if (row && col && value && PyLong_Check(row) && PyLong_Check(col) && PyFloat_Check(value)){
            // float vaPyObject *rslue
            set(self->mat, PyLong_AsLong(row), PyLong_AsLong(col), PyFloat_AsDouble(value));
        } else {
            // incorrect value
            PyErr_SetString(PyExc_TypeError, "Value input was not float/int.");
            return NULL;
        }
    } else {
        //throw type error
        PyErr_SetString(PyExc_TypeError, "Incorrect number of elements: does not have 3 (row, col, value).");
        return NULL;
    }

    // replace value using matrix.h's set(matrix *mat, int row, int col, double val)
    return Py_None;
}

/*
 * Given a numc.Matrix `self`, parse `args` to (int) row and (int) col.
 * Return the value at the `row`th row and `col`th column, which is a Python
 * float/int.
 */
PyObject *Matrix61c_get_value(Matrix61c *self, PyObject* args) {
    //parse args
    PyObject *r;
    PyObject *c;
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    //printf("Start");
    if (PyArg_UnpackTuple(args, "args", 2, 2, &arg1, &arg2)) { // change min args to 1 to handle 1-d one index case!!
        // arguments parsed
        PyObject *self_r;
        PyObject *self_c;

        if(PyTuple_Size(self->shape) == 1) {
            self_r = PyTuple_GetItem(self->shape, 0);
            self_c = PyLong_FromLong(1);
            r = arg1;
            if(!arg2) {
                c = PyLong_FromLong(0);
            } else {
                c = arg2;
            }
        } else if (PyTuple_Size(self->shape) == 2 && arg1 && arg2) {
            self_r = PyTuple_GetItem(self->shape, 0);
            self_c = PyTuple_GetItem(self->shape, 1);
            r = arg1;
            c = arg2;
        } else {
            // Input error, incorrect number of args
            PyErr_SetString(PyExc_TypeError,"Incorrect Number of Args");
            return NULL;
        }

        int rows = PyLong_AsLong(self_r);
        int cols = PyLong_AsLong(self_c);
        if (r && c && PyLong_Check(r) && PyLong_Check(c)) {
            if (PyLong_AsLong(r) >= rows || PyLong_AsLong(c) >= cols || PyLong_AsLong(r) < 0 ||
                PyLong_AsLong(c) < 0) {
                PyErr_SetString(PyExc_IndexError, "Index out of range.");
                return NULL;
            }

            double val = get(self->mat, PyLong_AsLong(r), PyLong_AsLong(c));
            PyObject *res = PyFloat_FromDouble(val);
            return res;
        } else {
            PyErr_SetString(PyExc_TypeError, "Input i,j values were not parsed/not ints.");
            return NULL;
        }
    } else {
        //throw type error //TODO : may not be type-- could be value?
        PyErr_SetString(PyExc_TypeError, "the number of arguments parsed from args is not 1 or 2.");
        return NULL;
    }
}

/*
 * Create an array of PyMethodDef structs to hold the instance methods.
 * Name the python function corresponding to Matrix61c_get_value as "get" and Matrix61c_set_value
 * as "set"
 * You might find this link helpful: https://docs.python.org/3.6/c-api/structures.html
 */
PyMethodDef Matrix61c_methods[] = {
    {"get", (PyCFunction)Matrix61c_get_value, METH_VARARGS, "get a row,col indexed value from matrix"},
    {"set", (PyCFunction)Matrix61c_set_value, METH_VARARGS, "set a row,col indexed valued from matrix" },
    {NULL, NULL, 0, NULL}
};

/* INDEXING */

/*
 * Given a numc.Matrix `self`, index into it with `key`. Return the indexed result.
 */
PyObject *Matrix61c_subscript(Matrix61c* self, PyObject* key) {
    // Tuple Pack has size, rows, (cols)
    PyObject *shape = self->shape;
    if (PyTuple_Size(shape) == 1) {
        // for 1d matrix: key can be an integer or single slice
        // unsupported slices should throw a type error
        int len = PyLong_AsLong(PyTuple_GetItem(shape, 0));
        //printf("LENGTH GOTTEN\n");
        if (PyLong_Check(key)) {// INTEGER CASE
            int index = PyLong_AsLong(key);
            //printf("INDEX GOTTEN\n");
            if (index < 0 || index >= len) {
                PyErr_SetString(PyExc_IndexError, "key index out of range. (1-d int)");
                return NULL;
            }
            // return the ith value of self->matrix
            PyObject *args = PyTuple_New(2);
            if (args == NULL) {
                PyErr_SetString(PyExc_RuntimeError, "cannot make new tuple.");
                return NULL;
            }
            //PyTuple_SetItem(args, 0, PyLong_FromLong(index));
            //printf("HELO1");
            //PyTuple_SetItem(args, 1, PyLong_FromLong(0));
            //printf("HELO2");
            //PyObject *rv = Matrix61c_get_value(self, args);
            //printf("HELO3");
            //return (PyObject*) rv;
            return PyFloat_FromDouble(get(self->mat, 0, index));
        } // INT CASE
        else if (PySlice_Check(key) == 1) { // SLICE CASE
            Py_ssize_t stop;
            Py_ssize_t start;
            Py_ssize_t step;
            Py_ssize_t slicelength;
            int len = PyLong_AsLong(PyTuple_GetItem(shape, 0));
            int pass = PySlice_GetIndicesEx(key, len, &start, &stop, &step, &slicelength);
            if (pass == -1) {
                PyErr_SetString(PyExc_RuntimeError, "indices not able to be parsed from slice");
                return NULL;
            }

            // TODO: do we need to check start stop index errors?
            if (step != 1) {
                PyErr_SetString(PyExc_ValueError, "step size != 1");
                return NULL;
            }

            if (slicelength < 1 || slicelength > len) {
                PyErr_SetString(PyExc_ValueError, "slice length < 1");
                return NULL;
            }

            Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            rv->shape = get_shape(slicelength, 1);
            matrix *new_mat;
            // TODO: may need fixin
            int alloc_ref_failed = allocate_matrix_ref(&new_mat, self->mat, start, 0, slicelength, 1);

            if (alloc_ref_failed) {
                PyErr_SetString(PyExc_RuntimeError, "Cannot malloc ref matrix. (1-d int)");
                return NULL;
            }
            rv->mat = new_mat;

            // TODO: may need to do the same thing as slice in 2-d
            return (PyObject*) rv;
        } // 1-d slice
        else {
            PyErr_SetString(PyExc_TypeError, "Key not valid! (Key was not an int or slice) (1-d)");
            return NULL;
        } // ERR
    }
    else if (PyTuple_Size(shape) == 2) { // 2d matrix
        // for 2d matrix: key can be a integer, slice, or tuple of two slices/ints
        if (PyLong_Check(key)) {// INTEGER CASE
            // GET i'TH ROW OF MATRIX.
            int row = PyLong_AsLong(key);
            if (row >= PyLong_AsLong(PyTuple_GetItem(shape, 0)) || row < 0) {
                PyErr_SetString(PyExc_RuntimeError, "Int index is out of range.");
                return NULL;
            }
            int cols = PyLong_AsLong(PyTuple_GetItem(shape, 1));

            Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            rv->shape = get_shape(1, cols);
            matrix *new_mat;
            int alloc_ref_failed = allocate_matrix_ref(&new_mat, self->mat, row, 0, 1, cols);
            if (alloc_ref_failed) {
                PyErr_SetString(PyExc_RuntimeError, "Cannot malloc ref matrix. (1-d int)");
                return NULL;
            }
            rv->mat = new_mat;
            //PyObject *new_key = PyTuple_New(2);
            //PyTuple_SetItem(new_key, 0, PyLong_FromLong(0));
            //PyObject * new_slice = PySlice_New(PyLong_FromLong(0), PyLong_FromLong(cols), PyLong_FromLong(1));
            //PyTuple_SetItem(new_key, 1, new_slice);
            //return (PyObject*) Matrix61c_subscript(self, new_key);
            return (PyObject *) rv;
        } // INT CASE
        else if (PySlice_Check(key)) { // 2-d SLICE CASE
            // GET SLICE ROWS FROM MATRIX
            Py_ssize_t stop;
            Py_ssize_t start;
            Py_ssize_t step;
            Py_ssize_t slicelength;
            PyObject *slice = key;
            int len = PyLong_AsLong(PyTuple_GetItem(shape, 0));
            int pass = PySlice_GetIndicesEx(slice, len, &start, &stop, &step, &slicelength);
            if (pass == -1) {
                PyErr_SetString(PyExc_RuntimeError, "indices not able to be parsed from slice");
                return NULL;
            }
            if (start < 0 || start >= len || stop <= 0 || stop > len) {
                PyErr_SetString(PyExc_ValueError, "slice indices out of bounds");
                return NULL;
            }

            if (step != 1) {
                PyErr_SetString(PyExc_ValueError, "step size != 1");
                return NULL;
            }

            if (slicelength < 1 || slicelength > len) {
                PyErr_SetString(PyExc_ValueError, "slice length < 1");
                return NULL;
            }
            int cols = PyLong_AsLong(PyTuple_GetItem(self->shape, 1));
            Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            rv->shape = get_shape(slicelength, cols);
            matrix *new_mat;
            int alloc_ref_failed = allocate_matrix_ref(&new_mat, self->mat, start, 0, slicelength, cols);
            if (alloc_ref_failed) {
                PyErr_SetString(PyExc_RuntimeError, "Cannot malloc ref matrix. (1-d int)");
                return NULL;
            }
            rv->mat = new_mat;
            return (PyObject*) rv;
        } // SLICE CASE
        else if (PyTuple_Check(key)){ // TUPLE CASE
            PyObject *key0 = PyTuple_GetItem(key, 0);
            PyObject *key1 = PyTuple_GetItem(key, 1);

            int rows = PyLong_AsLong(PyTuple_GetItem(shape, 0));
            int cols = PyLong_AsLong(PyTuple_GetItem(shape, 1));

            if (PyLong_Check(key0)) {
                int row = PyLong_AsLong(key0);
                if (row < 0 || row >= rows) {
                    PyErr_SetString(PyExc_IndexError, "key0 index out of bounds (2-d tuple)");
                    return NULL;
                }
                if (PyLong_Check(key1)) { // int, int
                    // int int case
                    int col = PyLong_AsLong(key1);
                    if (col < 0 || col >= cols) {
                        PyErr_SetString(PyExc_IndexError, "key1 index out of bounds (2-d int, int)");
                        return NULL;
                    }
                    PyObject *args = PyTuple_New(2);
                    if (args == NULL) {
                        PyErr_SetString(PyExc_RuntimeError, "failed to allocate new tuple");
                        return NULL;
                    }
                    PyTuple_SetItem(args, 0, PyLong_FromLong(row));
                    PyTuple_SetItem(args, 1, PyLong_FromLong(col));

                    PyObject *rv = Matrix61c_get_value(self, args);
                    return (PyObject *) rv;

                } // int, int
                else if (PySlice_Check(key1)) {
                    // int, slice case
                    Py_ssize_t stop1;
                    Py_ssize_t start1;
                    Py_ssize_t step1;
                    Py_ssize_t slicelength1;

                    int pass1 = PySlice_GetIndicesEx(key1, cols, &start1, &stop1, &step1, &slicelength1);
                    if (pass1 == -1) {
                        // throw error
                        PyErr_SetString(PyExc_RuntimeError, "indices not able to be parsed from key1");
                        return NULL;
                    } if (step1 != 1) {
                        PyErr_SetString(PyExc_ValueError, "step size != 1");
                        return NULL;
                    } if (slicelength1 < 1 || slicelength1 > rows) {
                        PyErr_SetString(PyExc_ValueError, "slicelength1 < 1 or > rows");
                        return NULL;
                    }

                    if (slicelength1 == 1) {
                        PyObject *rv = PyFloat_FromDouble(get(self->mat, row, start1));
                        return rv;
                    }
                    else {
                        Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                        rv->shape = get_shape(1, slicelength1);
                        // if slice is 0:1, get just the integer value.
                        matrix *new_mat;
                        int alloc_ref_failed = allocate_matrix_ref(&new_mat, self->mat, row, start1, 1, slicelength1);
                        if (alloc_ref_failed) {
                            PyErr_SetString(PyExc_RuntimeError, "Cannot malloc ref matrix. (2-d int)");
                            return NULL;
                        }
                        rv->mat = new_mat;
                        return (PyObject *) rv;
                    }



                } // int, slice
                else {
                    PyErr_SetString(PyExc_TypeError, "key1 not int or slice (2-d tuple)");
                    return NULL;
                } // err
            } //key0 = int
            else if (PySlice_Check(key0)) {
                Py_ssize_t stop0;
                Py_ssize_t start0;
                Py_ssize_t step0;
                Py_ssize_t slicelength0;

                int pass0 = PySlice_GetIndicesEx(key0, rows, &start0, &stop0, &step0, &slicelength0);
                if (pass0 == -1) {
                    // throw error
                    PyErr_SetString(PyExc_RuntimeError, "indices not able to be parsed from key0");
                    return NULL;
                } if (step0 != 1) {
                    PyErr_SetString(PyExc_ValueError, "step size != 1");
                    return NULL;
                } if (slicelength0 < 1 || slicelength0 > rows) {
                    PyErr_SetString(PyExc_ValueError, "slicelength0 < 1 or > rows");
                    return NULL;
                }

                if (PyLong_Check(key1)) {
                    // slice, int case
                    int c = PyLong_AsLong(key1);
                    if (c < 0 || c >= cols) {
                        PyErr_SetString(PyExc_IndexError, "Index out of bounds for key 1. (2-d slice, int)");
                        return NULL;
                    }
                    if (slicelength0 == 1) {
                        PyObject *rv = PyFloat_FromDouble(get(self->mat, start0, c));
                        return rv;
                    } else {
                        Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                        rv->shape = get_shape(slicelength0, 1);
                        matrix *new_mat;
                        int alloc_ref_failed = allocate_matrix_ref(&new_mat, self->mat, start0, c, slicelength0, 1);
                        if (alloc_ref_failed) {
                            PyErr_SetString(PyExc_RuntimeError, "Cannot malloc ref matrix. (2-d int)");
                            return NULL;
                        }
                        rv->mat = new_mat;
                        return (PyObject *) rv;
                    }
                } //slice int
                else if (PySlice_Check(key1)) { // slice, slice
                    Py_ssize_t stop1;
                    Py_ssize_t start1;
                    Py_ssize_t step1;
                    Py_ssize_t slicelength1;

                    int pass1 = PySlice_GetIndicesEx(key1, cols, &start1, &stop1, &step1, &slicelength1);
                    if (pass1 == -1) {
                        // throw error
                        PyErr_SetString(PyExc_RuntimeError, "indices not able to be parsed from key1");
                        return NULL;
                    } if (step1 != 1) {
                        PyErr_SetString(PyExc_ValueError, "step size != 1");
                        return NULL;
                    } if (slicelength1 < 1 || slicelength1 > rows) {
                        PyErr_SetString(PyExc_ValueError, "slicelength1 < 1 or > rows");
                        return NULL;
                    }

                    if(slicelength0 == 1 && slicelength1 == 1) {
                        return (PyObject *) PyFloat_FromDouble(get(self->mat, start0, start1));
                    }
                    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                    rv->shape = get_shape(slicelength0, slicelength1);
                    matrix *new_mat;
                    int alloc_ref_failed = allocate_matrix_ref(&new_mat, self->mat, start0, start1, slicelength0,
                                                               slicelength1);
                    if (alloc_ref_failed) {
                        PyErr_SetString(PyExc_RuntimeError, "Cannot malloc ref matrix. (2-d int)");
                        return NULL;
                    }
                    rv->mat = new_mat;
                    return (PyObject *) rv;
                } // slice slice
                else {
                    PyErr_SetString(PyExc_TypeError, "not a valid key1 input: not int or slice");
                    return NULL;
                }
            } //key0 = slice
            else {
                PyErr_SetString(PyExc_TypeError, "not a valid key0 input: not int or slice");
                return NULL;
            }
        } // TUPLE CASE
        else {
            PyErr_SetString(PyExc_TypeError, "not a valid key input: not int, slice, or tuple.");
            return NULL;
        }
    } // 2-d cases

    PyErr_SetString(PyExc_TypeError, "Invalid arguments: Not a valid matrix.");
    return NULL;
    // value error if key is a slice or tuple that has one or more slices, and atleast one
    // of these slices has a step size not equal to 1 or if slice len is <1
}


/*
 * Given a numc.Matrix `self`, index into it with `key`, and set the indexed result to `v`.
 */
int Matrix61c_set_subscript(Matrix61c* self, PyObject *key, PyObject *v) {
    // check shape of mat

    if (PyTuple_Size(self->shape) == 1) { // matrix is 1-d
        int len = PyLong_AsLong(PyTuple_GetItem(self->shape, 0)); // len of 1-d matrix
        if (PyLong_Check(key)) { // key is an integer
            int i = PyLong_AsLong(key); // convert key to int value.
            if (i < 0 || i >= len) {
                // THROW VALUE ERROR: Index out of bounds
                PyErr_SetString(PyExc_IndexError, "Index Out of Bounds for 1-d Array");
                return -1; // return -1 on failure
            }
            if (!PyLong_Check(v) && !PyFloat_Check(v)) {
                // for 1-d matrix, the value needs to be a single float or int
                PyErr_SetString(PyExc_TypeError, "Value is not a float/int.");
                return -1; // return -1 on failure
            }
            if (PyLong_Check(v)) {
                set(self->mat, 0, i, PyLong_AsDouble(v));
            } else if (PyFloat_Check(v)) {
                set(self->mat, 0, i, PyFloat_AsDouble(v));
            }
            return 0; // zero on success
        } // key is an int
        else if (PySlice_Check(key)) { // key is a slice
            Py_ssize_t stop;
            Py_ssize_t start;
            Py_ssize_t step;
            Py_ssize_t slicelength;
            int sliced = PySlice_GetIndicesEx(key, len, &start, &stop, &step, &slicelength);
            if (sliced == -1) {
                // THROW RUNTIME ERROR: Failed to slice
                PyErr_SetString(PyExc_RuntimeError, "Failed to Slice 1-d slice key.");
                return -1;
            }
            if (step != 1) {
                // THROW ERROR: Slice not valid.
                PyErr_SetString(PyExc_ValueError, "step size != 1 for 1-d slice");
                return -1;
            }
            if (slicelength > len || slicelength < 1) {
                // not valid slice.
                PyErr_SetString(PyExc_ValueError, "slicelength not valid for 1-d slice");
                return -1;
            }
            if (start < 0 || start >= len || stop < 0 || stop > len) {
                PyErr_SetString(PyExc_ValueError, "start and stop are out of bounds for length of 1-d array.");
                return -1;
            }
            if (slicelength == 1 && !PyLong_Check(v) && !PyFloat_Check(v)) {
                // Type Error: value must be a float or int.
                PyErr_SetString(PyExc_TypeError, "input value must be a float or int. (1x1 1-d slice)");
                return -1;
            }

            if (PyLong_Check(v)) { // VALUE IS AN INT OR FLOAT
                double val = PyLong_AsDouble(v);
                for (int i = start; i < stop; i += step) {
                    set(self->mat, 0, i, val);
                }
            } else if (PyFloat_Check(v)) {
                double val = PyFloat_AsDouble(v);
                for (int i = start; i < stop; i += step) {
                    set(self->mat, 0, i, val);
                }
            } else if (PyList_Check(v)) { // VALUE IS A LIST
                if (PyList_Size(v) != slicelength && slicelength > 1) {
                    // VALUE ERROR
                    PyErr_SetString(PyExc_ValueError, "value must be a list of same size as 1-d slice.");
                    return -1;
                }
                int count = 0;
                for (int i = start; i < stop; i += step) {
                    PyObject *item = PyList_GetItem(v, count);
                    if (PyLong_Check(item)) {
                        set(self->mat, 0, i, PyLong_AsDouble(item));
                    } else if (PyFloat_Check(item)) {
                        set(self->mat, 0, i, PyFloat_AsDouble(item));
                    } else {
                        // VALUE ERROR: Item in v is not float or int
                        PyErr_SetString(PyExc_ValueError, "item in value is not float or int.");
                        return -1;
                    }
                    count = count + 1;
                }
            } else {
                // TYPE ERROR: Slice Length is > 1 but v is not a list.
                PyErr_SetString(PyExc_TypeError, "value must be a list for slice length > 1");
                return -1;
            }
        } // key is a slice
        else {
            // THROW TYPE ERROR
            PyErr_SetString(PyExc_TypeError, "Key is not a slice or int (1-d)");
            return -1;
        } // ERR
    } // MATRIX is 1-d
    else if (PyTuple_Size(self->shape) == 2) { // matrix is 2-d
        // TODO: DO 2-d version!!!
        // get rows and cols from self
        PyObject *rs = PyTuple_GetItem(self->shape, 0);
        PyObject *cs = PyTuple_GetItem(self->shape, 1);
        int rows = PyLong_AsLong(rs);
        int cols = PyLong_AsLong(cs);
        int pass;
        if (PyLong_Check(key)) {// if key is int
            int i = PyLong_AsLong(key);
            // set subscript for the key-th row (1-d matrix) of self->mat
            // make sure that v is an int or list of equal length as # cols
            if (i < 0 || i >= rows) {
                PyErr_SetString(PyExc_IndexError, "key was out of bounds (2-d) row int");
            }
            if (PyLong_Check(v)) {
                // TODO: maybe need to throw error instead>>>
                for (int j = 0; j < cols; j++) {
                    set(self->mat, i, j, PyLong_AsDouble(v));
                }
            } else if (PyFloat_Check(v)) {
                // TODO: maybe not valid.
                for (int j = 0; j < cols; j++) {
                    set(self->mat, i, j, PyFloat_AsDouble(v));
                }
            } else if (PyList_Check(v)) {
                // check list is equal length as cols
                if (PyList_Size(v) != cols) {
                    // incorrect array len
                    PyErr_SetString(PyExc_TypeError, "value list is not the correct size.");
                    return -1;
                }
                for (int j = 0; j < cols; j++) {
                    PyObject *item = PyList_GetItem(v, j);
                    if (PyLong_Check(item)) {
                        set(self->mat, i, j, PyLong_AsDouble(item));
                    } else if (PyFloat_Check(item)) {
                        set(self->mat, i, j, PyFloat_AsDouble(item));
                    } else {
                        // VALUE ERR: v is not a float or int
                        PyErr_SetString(PyExc_ValueError, "item in v is not float/int");
                        return -1;
                    }
                }
            } else {
                // THROW INPUT ERROR;
                PyErr_SetString(PyExc_TypeError, "value not list/int/float. (2-d single key)");
                return -1;
            }
        } // key is an int
        else if (PySlice_Check(key)) { // if key is slice
            // --> this is a slice for the rows of the matrix ie. 0:2 is rows 0, 1
            // make sure v is list with length = cols or int/float
            // for each row of slice call Matrix_set_subscript(mat, row_index, val);
            Py_ssize_t stop;
            Py_ssize_t start;
            Py_ssize_t step;
            Py_ssize_t slicelength;
            pass = PySlice_GetIndicesEx(key, rows, &start, &stop, &step, &slicelength);
            if (pass == -1) {
                PyErr_SetString(PyExc_RuntimeError, "slice not found. (2-d single key).");
                return -1;
            }
            if (stop > rows || stop < 0) {
                //printf("stop is : %ld\n", stop);
                PyErr_SetString(PyExc_ValueError, "stop out of bounds");
                return -1;
            }
            if (start >= rows || start < 0) {
                PyErr_SetString(PyExc_ValueError, "start out of bounds");
                return -1;
            }
            if (step != 1) {
                PyErr_SetString(PyExc_ValueError, "step size != 1 (2-d single slice)");
                return -1;
            }
            if (slicelength < 1 || slicelength > rows) {
                PyErr_SetString(PyExc_ValueError, "Not a valid slice; slicelength < 0 (2-d single slice)");
                return -1;
            }
            if (PyLong_Check(v)) {
                for (int i = start; i < stop; i++) {
                    for (int j = 0; j < cols; j++) {
                        set(self->mat, i, j, PyLong_AsDouble(v));
                    }
                }
            } else if (PyFloat_Check(v)) {
                for (int i = start; i < stop; i++) {
                    for (int j = 0; j < cols; j++) {
                        set(self->mat, i, j, PyFloat_AsDouble(v));
                    }
                }
            } else if (PyList_Check(v)) {

                if (PyList_Size(v) != cols) {
                    PyErr_SetString(PyExc_ValueError, "v is not a list with correct length.");
                    return -1;
                }
                for (int i = start; i < stop; i++) {
                    //printf("i: %d\n", i);
                    for (int j = 0; j < cols; j++) {
                        //printf("j: %d\n", j);
                        PyObject *val = PyList_GetItem(v, j);
                        //printf("pylist v size: %ld\n", PyList_Size(v));
                        //printf("val begotten!\n");
                        if (PyList_Check(val)) {
                            PyErr_SetString(PyExc_ValueError, "for 2-d single splice, v cannot be a matrix.");
                            return -1;
                        } else {
                            if (PyLong_Check(val)) {
                                double value = PyLong_AsDouble(val);
                                //printf("value: %f\n", value);
                                set(self->mat, i, j, value);
                            } else if (PyFloat_Check(val)) {
                                double value = PyFloat_AsDouble(val);
                                //printf("value: %f\n", value);
                                set(self->mat, i, j, value);
                            } else {
                                PyErr_SetString(PyExc_ValueError, "value in v is not float/int.");
                                return -1;
                            }
                        }
                    }
                }
            }
        }  // key is a slice
        else if (PyTuple_Check(key)) {
            PyObject *key0 = PyTuple_GetItem(key, 0);
            PyObject *key1 = PyTuple_GetItem(key, 1);
            if (PyLong_Check(key0)) {
                int k0 = PyLong_AsLong(key0);
                if (k0 < 0 || k0 > rows) {
                    PyErr_SetString(PyExc_IndexError, "k0 out of index. (2-d int,int)");
                    return -1;
                }
                if (PyLong_Check(key1)) {
                    // int, int
                    int k1 = PyLong_AsLong(key1);
                    if (k1 < 0 || k1 > cols) {
                        PyErr_SetString(PyExc_IndexError, "k1 out of index. (2-d int,int)");
                        return -1;
                    }
                    if (PyLong_Check(v)) {
                        set(self->mat, k0, k1, PyLong_AsDouble(v));
                    } else if (PyFloat_Check(v)) {
                        set(self->mat, k0, k1, PyFloat_AsDouble(v));
                    } else {
                        // VALUE ERR? TYPE ERR? value is not an int or float
                        PyErr_SetString(PyExc_ValueError, "v is not int/float (2-d int,int)");
                        return -1;
                    }

                } // int, int case
                else if (PySlice_Check(key1)) {
                    // int, slice
                    Py_ssize_t stop;
                    Py_ssize_t start;
                    Py_ssize_t step;
                    Py_ssize_t slicelength;
                    pass = PySlice_GetIndicesEx(key1, cols, &start, &stop, &step, &slicelength);
                    if (pass == -1) {
                        PyErr_SetString(PyExc_RuntimeError, "slice not found (2-d int, slice)");
                        return -1;
                    }
                    if (stop > cols || stop < 0) {
                        PyErr_SetString(PyExc_ValueError, "stop out of bounds");
                        return -1;
                    }
                    if (start >= cols || start < 0) {
                        PyErr_SetString(PyExc_ValueError, "start out of bounds");
                        return -1;
                    }
                    if (step != 1) {
                        PyErr_SetString(PyExc_ValueError, "step size != 1 (2-d int, slice)");
                        return -1;
                    }
                    if (slicelength < 1 || slicelength > cols) {
                        PyErr_SetString(PyExc_ValueError, "slice length < 1 or > cols(2-d int, slice)");
                        return -1;
                    }

                    if (PyLong_Check(v)) {
                        for (int j = start; j < stop; j++) {
                            set(self->mat, k0, j, PyLong_AsDouble(v));
                        }
                    } else if (PyFloat_Check(v)) {
                        for (int j = start; j < stop; j++) {
                            set(self->mat, k0, j, PyFloat_AsDouble(v));
                        }
                    } else if (PyList_Check(v)) {
                        //printf("just to see::: slen = %ld, rows = %d, listsize = %ld\n", slicelength, rows, PyList_Size(v));
                        if (PyList_Size(v) != rows) {
                            PyErr_SetString(PyExc_ValueError, "list size is not same size as rows (2-d int, slice)");
                            return -1;
                        }
                        for (int j = start; j < stop; j++) {
                            if (PyFloat_Check(PyList_GetItem(v, j))) {
                                set(self->mat, k0, j, PyFloat_AsDouble(PyList_GetItem(v, j)));
                            } else if (PyLong_Check(PyList_GetItem(v, j))) {
                                set(self->mat, k0, j, PyLong_AsDouble(PyList_GetItem(v, j)));
                            } else {
                                // jth value of value list 'v' is not an int or float
                                PyErr_SetString(PyExc_ValueError, "value in v is not int/float (2-d int, slice)");
                                return -1;
                            }
                        }
                    }
                } // int, slice case
            } // key0 is an int
            else if (PySlice_Check(key0)) {
                Py_ssize_t stop;
                Py_ssize_t start;
                Py_ssize_t step;
                Py_ssize_t slicelength;
                pass = PySlice_GetIndicesEx(key0, rows, &start, &stop, &step, &slicelength);
                if (pass == -1) {
                    PyErr_SetString(PyExc_RuntimeError, "slice not found. (2-d slice, int)");
                    return -1;
                }
                if (stop > rows || stop < 0) {
                    PyErr_SetString(PyExc_ValueError, "stop out of bounds (2-d slice, int)");
                    return -1;
                }
                if (start >= rows || start < 0) {
                    PyErr_SetString(PyExc_ValueError, "start out of bounds (2-d slice, int)");
                    return -1;
                }
                if (step != 1) {
                    PyErr_SetString(PyExc_ValueError, "step != 1 (2-d slice, int)");
                    return -1;
                }
                if (slicelength < 1 || slicelength > rows) {
                    PyErr_SetString(PyExc_ValueError, "slicelength out of bounds (2-d slice, int)");
                    return -1;
                }
                if (PyLong_Check(key1)) {
                    // slice, int
                    int k1 = PyLong_AsLong(key1);
                    if (k1 < 0 || k1 > cols) {
                        PyErr_SetString(PyExc_IndexError, "k1 out of index. (2-d int,int)");
                        return -1;
                    }
                    if (PyList_Check(v)) {
                        if (slicelength != PyList_Size(v) || slicelength == 1) {
                            // length of v must equal slice length
                            PyErr_SetString(PyExc_TypeError, "length of v must equal slice len (2-d slice, int)");
                            return -1;
                        }
                        for (int i = start; i < stop; i++) {//for each row in sliced rows:
                            // set key1'th value as the count'th value of list 'v'
                            if (PyLong_Check(PyList_GetItem(v, i))) {
                                set(self->mat, i, k1, PyLong_AsDouble(PyList_GetItem(v, i)));
                            } else if (PyFloat_Check(PyList_GetItem(v, i))) {
                                set(self->mat, i, k1, PyFloat_AsDouble(PyList_GetItem(v, i)));
                            } else {
                                // error: val of list 'v' is not an int or float.
                                PyErr_SetString(PyExc_ValueError, "value in v is not int/float (2-d slice, int)");
                                return -1;
                            }
                        }
                    } else if (PyLong_Check(v)) {
                        for (int i = start; i < stop; i++) {//for each row in sliced rows:
                            // set key1'th value as the int val 'v'
                            set(self->mat, i, k1, PyLong_AsDouble(v));
                        }
                    } else if (PyFloat_Check(v)) {
                        for (int i = start; i < stop; i++) {//for each row in sliced rows:
                            // set key1'th value as the int val 'v'
                            set(self->mat, i, k1, PyFloat_AsDouble(v));
                        }
                    } else {
                        // 'v' is not list, or int/float.
                        PyErr_SetString(PyExc_ValueError, "v is not list/int/float (2-d slice, int)");
                        return -1;
                    }

                } // slice, int
                else if (PySlice_Check(key1)) {
                    // slice, slice
                    Py_ssize_t stop1;
                    Py_ssize_t start1;
                    Py_ssize_t step1;
                    Py_ssize_t slicelength1;
                    int pass1 = PySlice_GetIndicesEx(key1, cols, &start1, &stop1, &step1, &slicelength1);
                    if (pass1 == -1) {
                        PyErr_SetString(PyExc_RuntimeError, "slice not found (2-d slice, slice)");
                        return -1;
                    }
                    if (stop1 > cols || stop1 < 0) { // slicelength1, slicelength2 must equal the shape of 'v'if (stop1 >= cols|| stop < 0) {
                        PyErr_SetString(PyExc_ValueError, "stop1 out of bounds (2-d slice, slice)");
                        return -1;
                    }
                    if (start1 >= cols || start < 0) {
                        PyErr_SetString(PyExc_ValueError, "start1 out of bounds (2-d slice, slice)");
                        return -1;
                    }
                    if (step1 != 1) {
                        PyErr_SetString(PyExc_ValueError, "step1 != 1 (2-d slice, slice)");
                        return -1;
                    }
                    if (slicelength1 < 1 || slicelength1 > cols) {
                        PyErr_SetString(PyExc_ValueError, "slicelength1 out of bounds (2-d slice, slice)");
                        return -1;
                    }


                    if (PyList_Check(PyList_GetItem(v, 0)) == 0) {
                        // create a new 'v' list with duplicate v's
                        //last_coded
                        //printf("oh look an ouija bouard-- v size= %ld\n", PyList_Size(v));
                        //printf("HELO\n");
                        PyObject *temp_v = PyList_New(slicelength1);
                        if (temp_v == NULL) {
                            PyErr_SetString(PyExc_RuntimeError, "failed to make new temp_v (2-d slice, slice)");
                            return -1;
                        }
                        //printf("How are you?\n");
                        int set;
                        for (int tempi = 0; tempi < slicelength1; tempi++) {
                            set = PyList_SetItem(temp_v, tempi, v);
                            if (set != 0) {
                                PyErr_SetString(PyExc_RuntimeError, "failed to set temp_v with v (2-d slice, slice)");
                                return -1;
                            }
                           // printf("give us a sign\n");
                        }
                        v = temp_v;
                        //printf("*cue light flicker*\n");
                    }

                    Py_ssize_t list_rows = PyList_Size(v);
                    //printf("your sanity is so low-- list_rows = %ld\n", list_rows);
                    Py_ssize_t list_cols = PyList_Size(PyList_GetItem(v, 0));
                    //printf("%ld, %ld, %ld, %ld", list_rows, list_cols, slicelength, slicelength1);
                    if (list_rows <= 0 || list_cols <= 0 || list_rows != slicelength || list_cols != slicelength1) {
                        // nested list 'v' does not have same shape or # inputs as needed for 2d splice.
                        PyErr_SetString(PyExc_ValueError, "v is not in correct shape (2-d slice, slice)");
                        return -1;
                    }
                    int lr = 0;
                    int lc = 0;
                    for (int i = start; i < stop; i++) {
                        for (int j = start1; j < stop1; j++) {
                            PyObject *val = PyList_GetItem(PyList_GetItem(v, lr), lc);
                            if (PyLong_Check(val)) {
                                set(self->mat, i, j, PyLong_AsDouble(val));
                            } else if (PyFloat_Check(val)) {
                                set(self->mat, i, j, PyFloat_AsDouble(val));
                            } else {
                                // VALUE ERROR
                                PyErr_SetString(PyExc_ValueError, "value in v is not int/float (2-d slice, slice)");
                                return -1;
                            }
                            lc = lc + 1;
                        }
                        lc = 0;
                        lr = lr + 1;
                    }
                } // slice, slice
                else {
                    PyErr_SetString(PyExc_TypeError, "key1 is not int or slice (2-d tuple)");
                    return -1;
                } // ERR
            } // key0 is a slice
            else {
                PyErr_SetString(PyExc_TypeError, "key0 is not int/slice (2-d tuple)");
                return -1;
            } // ERR
        }  // key is a tuple
        else {
            // THROW VALUE ERROR --> too many args in shape.
            PyErr_SetString(PyExc_TypeError, "key is not int/slice/tuple (2-d)");
            return -1;
        }
    } // MATRIX is 2-d
    else {
        // THROW VALUE ERROR --> too many args in shape.
        PyErr_SetString(PyExc_ValueError, "too many args in self->shape");
        return -1;
    } // self->shape not valid :: invalid matrix error::

    return 0; // ON SUCCESS!
}

PyMappingMethods Matrix61c_mapping = {
    NULL,
    (binaryfunc) Matrix61c_subscript,
    (objobjargproc) Matrix61c_set_subscript,
};

/* INSTANCE ATTRIBUTES*/
PyMemberDef Matrix61c_members[] = {
    {
        "shape", T_OBJECT_EX, offsetof(Matrix61c, shape), 0,
        "(rows, cols)"
    },
    {NULL}  /* Sentinel */
};

PyTypeObject Matrix61cType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "numc.Matrix",
    .tp_basicsize = sizeof(Matrix61c),
    .tp_dealloc = (destructor)Matrix61c_dealloc,
    .tp_repr = (reprfunc)Matrix61c_repr,
    .tp_as_number = &Matrix61c_as_number,
    .tp_flags = Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_BASETYPE,
    .tp_doc = "numc.Matrix objects",
    .tp_methods = Matrix61c_methods,
    .tp_members = Matrix61c_members,
    .tp_as_mapping = &Matrix61c_mapping,
    .tp_init = (initproc)Matrix61c_init,
    .tp_new = Matrix61c_new
};


struct PyModuleDef numcmodule = {
    PyModuleDef_HEAD_INIT,
    "numc",
    "Numc matrix operations",
    -1,
    Matrix61c_class_methods
};

/* Initialize the numc module */
PyMODINIT_FUNC PyInit_numc(void) {
    PyObject* m;

    if (PyType_Ready(&Matrix61cType) < 0)
        return NULL;

    m = PyModule_Create(&numcmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&Matrix61cType);
    PyModule_AddObject(m, "Matrix", (PyObject *)&Matrix61cType);
    printf("CS61C Fall 2020 Project 4: numc imported!\n");
    fflush(stdout);
    return m;
}