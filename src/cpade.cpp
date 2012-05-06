/* 
 * File:   expmat.cpp
 * Author: Robert Gantner
 *
 * Created on December 01, 2011
 */

// Eigen Matrix Exponential in:
//   /usr/include/eigen3/unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h
// Python Pade implementation in:
//   /usr/lib/python2.7/site-packages/scipy/linalg/matfuncs.py

#ifndef EXPM_HPP
#define EXPM_HPP

#include "globalinc.hpp"
//#include <boost/shared_ptr.hpp>
#include <unsupported/Eigen/MatrixFunctions>

#ifdef PYTHONMODULE
    #include <Python.h>
    #include <boost/python.hpp>
    #include <numpy/arrayobject.h>
#endif

#include <Eigen/Core>
using namespace Eigen;

// delete this when done testing
#include <iostream>
using namespace std;

/**
 * Pade approximation of matrix exponential. Calls built-in Eigen function
 * @param A nxn matrix to exponentiate
 * @param expA nxn matrix into which result is to be stored
 */
template<class Derived>
void expm(const MatrixBase<Derived>& A, MatrixBase<Derived>& expA){
    expA = A.exp();
}

#ifdef PYTHONMODULE
/**
 * wrapper function that creates Eigen matrices of the correct type and
 * calls the expm function of the corresponding type.
 * @param A_in numpy array of dimension nxn
 * @param v_in numpy array of dimension nx1
 * @param k number of iterations
 * @param V_out numpy array of dimension nxk
 * @param H_out numpy array of dimension kxk
 */
template<class T>
void call_expm(PyObject* A_in, PyObject* expA_out) {
    npy_intp* shape = PyArray_DIMS(A_in);
    int n = shape[0];
    Map< Eigen::Matrix< T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor > > _A_in((T *) PyArray_DATA(A_in), n,n);
    Map< Eigen::Matrix< T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor > > _expA_out((T *) PyArray_DATA(expA_out), n,n);
    // call arnoldi algorithm. This function is also templated!
    expm(_A_in, _expA_out);
}

void expm_py(PyObject* A_in, PyObject* expA_out) {
    // sanity checks
    bool error = false;
    npy_intp* shape;
    int d;
    // A_in
    d = PyArray_NDIM(A_in);
    shape = PyArray_DIMS(A_in);
    if (d != 2 || shape[0] != shape[1]){
        cout << "arnoldi_py error: A is not a matrix in R^(n x n)" << endl;
        error = true;
    }
    int n = shape[0];
    // expA_out
    d = PyArray_NDIM(expA_out);
    shape = PyArray_DIMS(expA_out);
    if (d == 0 || shape[0] != n || shape[1] != n){
        cout << "arnoldi_py error: expA is not a matrix in R^(n x n), n=A.shape[0]" << endl;
        error = true;
    }
    
    if (error) return;
    
    // determine numpy data type
    PyArray_Descr* descr = PyArray_DESCR(A_in);
    switch (descr->type) {
    // float32 (single precision)
    case 'f':
        call_expm<float>(A_in, expA_out);
    break;
    
    // float64 (double precision)
    case 'd':
        call_expm<double>(A_in, expA_out);
    break;
    
    // complex64 (2x single precision)
    case 'F':
        call_expm<std::complex<float> >(A_in, expA_out);
    break;
    
    // complex128 (2x double precision)
    case 'D':
        call_expm<std::complex<double> >(A_in, expA_out);
    break;
        
//// unsupported types
//    case 'e': // float 16
//    case 'g': // floa128 (quadruple precision)
//    case 'G': // complex256 (2x quadruple precision)
        
        // default
    default:
        cout << "expm_py error: unknown type: " << descr->type << endl;
        return;
    }
}

namespace bp = boost::python;
#include "boost/python.hpp"
BOOST_PYTHON_MODULE(cpade) {
    bp::numeric::array::set_module_and_type("numpy", "ndarray");
    
    boost::python::def("expm", expm_py, "Compute the matrix exponential using Pade approximation.\n    Parameters\n        A : array, shape(M,M)\n            Matrix to be exponentiated\n\n        expA : array, shape(M,M)\n            Matrix exponential of A (overwritten by function)");
}
#endif


int main (int argc, char const *argv[])
{
    size_t n = 4;
    
    // initialize A
    Types<double>::Matrix A(n,n);
    for (unsigned int i=0; i<n; ++i) {
        A(i,i) = 2;
    }
    for (unsigned int i=1; i<n; ++i) {
        A(i,i-1) = -1;
        A(i-1,i) = -1;
    }
    
    // calculate exponential with expm function
    Types<double>::Matrix expA;
    expm(A,expA);
    
    // output
    cout << "A:\n" << A << endl;
    cout << "A.exp():\n" << A.exp() << endl;
    cout << "expm:\n" << expA << endl;
    
    return 0;
}


#endif    /* EXPM_HPP */
