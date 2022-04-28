#ifndef PYUTIL_H
#define PYUTIL_H

#include <string>
#include <Eigen/Core>
#include <Python.h>


void print_py_obj(PyObject* pyobj);
PyObject* load_npz(std::string fname);
Eigen::MatrixXd npa2mat(PyObject* npa);

#endif
