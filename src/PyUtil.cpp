#include <cstdlib>
#include "PyUtil.h"

void print_py_obj(PyObject *pyobj)
{
    PyObject_Print(pyobj, stdout, Py_PRINT_RAW);
}

PyObject* load_npz(std::string fname)
{
    PyObject *pName, *pModule, *pFunc, *pArgs, *pValue;
    const char *module = "numpy";
    const char *func = "load";
    pName = PyUnicode_FromString(module);
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL)
    {
	pFunc = PyObject_GetAttrString(pModule, func);
	if (pFunc && PyCallable_Check(pFunc))
	{
	    pValue = PyUnicode_FromString(fname.c_str());
	    pArgs = PyTuple_New(1);
	    PyTuple_SetItem(pArgs, 0, pValue);
	    pValue = PyObject_CallObject(pFunc, pArgs);
	    Py_DECREF(pArgs);
	    Py_DECREF(pFunc);
	    Py_DECREF(pModule);
	    return pValue;
	}
	else
	{
	    if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "failed to call function \"%s\"\n", func);
	}
	Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else
    {
	PyErr_Print();
        fprintf(stderr, "failed to load \"%s\"\n", module);
	return NULL;
    }
    return NULL;
}

Eigen::MatrixXd npa2mat(PyObject *npa)
{
    Eigen::MatrixXd mat;
    PyObject *pValue;
    int n, m;

    PyObject *pShape = PyObject_GetAttrString(npa, "shape");
    //print_py_obj(pShape);
    pValue = PyTuple_GetItem(pShape, 0);
    n = PyLong_AsSize_t(pValue);    // number of rows
    //Py_DECREF(pValue);
    pValue = PyTuple_GetItem(pShape, 1);
    m = PyLong_AsSize_t(pValue);    // number of columns
    //Py_DECREF(pValue);
    Py_DECREF(pShape);

    mat = Eigen::MatrixXd(n, m);
    for (int i = 0; i < n; ++i)
    {
	for (int j = 0; j < m; ++j)
	{
	    pValue = PyObject_CallMethod(npa, "item", "ii", i, j);
	    double d = PyFloat_AsDouble(pValue);
	    Py_DECREF(pValue);
	    mat(i,j) = d;
	}
    }

    return mat;
}

std::vector<Eigen::MatrixXd> npa2mat_vector(PyObject *npa)
{
    PyObject *pValue;
    int n, m, l;

    PyObject *pShape = PyObject_GetAttrString(npa, "shape");
    //print_py_obj(pShape);
    pValue = PyTuple_GetItem(pShape, 0);
    l = PyLong_AsSize_t(pValue);    // number of matrices
    pValue = PyTuple_GetItem(pShape, 1);
    n = PyLong_AsSize_t(pValue);    // number of rows
    //Py_DECREF(pValue);
    pValue = PyTuple_GetItem(pShape, 2);
    m = PyLong_AsSize_t(pValue);    // number of columns
    //Py_DECREF(pValue);
    Py_DECREF(pShape);

    std::vector<Eigen::MatrixXd> list;
    for (int k = 0; k < l; ++k)
    {
	Eigen::MatrixXd mat = Eigen::MatrixXd(n, m);
	for (int i = 0; i < n; ++i)
	{
	    for (int j = 0; j < m; ++j)
	    {
		pValue = PyObject_CallMethod(npa, "item", "iii", k, i, j);
		double d = PyFloat_AsDouble(pValue);
		Py_DECREF(pValue);
		mat(i,j) = d;
	    }
	}
	list.push_back(mat);
    }

    return list;
}
