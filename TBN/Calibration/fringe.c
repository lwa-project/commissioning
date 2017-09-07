#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <cblas.h>
#include <stdlib.h>
#include <complex.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "numpy/arrayobject.h"

#define PI 3.1415926535898
#define imaginary _Complex_I

/*
  Simple Fringing of Complex data
*/

static PyObject *Simple(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *visOut;
	PyArrayObject *data, *vis;
	
	int refX, refY;
	float clipLevel;
	long i, j, k, nStand, nSamps;
	
	if(!PyArg_ParseTuple(args, "Oiif", &signals, &refX, &refY, &clipLevel)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 2-D complex64");
		return NULL;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[1];
	dims[0] = (npy_intp) nStand;
	vis = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_COMPLEX64, 0);
	if(vis == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	int ref;
	
	Py_BEGIN_ALLOW_THREADS
	
	// Pointers
	float complex *a;
	float complex *b;
	a = (float complex *) PyArray_DATA(data);
	b = (float complex *) PyArray_DATA(vis);
	double complex tempV;
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(j, k, ref, tempV)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(i=0; i<nStand; i++) {
			// Which reference to use based on polarization ordering
			if( i % 2 == 0 ) {
				ref = refX;
			} else {
				ref = refY;
			}
			
			// Go!
			k = 0;
			tempV = 0.0;
			for(j=0; j<nSamps; j++) {
				// Bad input value
				if(cabs(*(a + (nSamps*i + j))) >= clipLevel) {
					continue;
				}
				
				// Bad reference input value
				if(cabs(*(a + (nSamps*ref + j))) >= clipLevel) {
					continue;
				}
				
				tempV += *(a + (nSamps*i + j)) * conj(*(a + (nSamps*ref + j)));
				k++;
			}
			
			// Average
			if( k > 0 ) {
				*(b + i) = tempV / (float) k;
			} else {
				*(b + i) = 0.0;
			}
		}
	}
	
	Py_END_ALLOW_THREADS
	
	visOut = Py_BuildValue("O", PyArray_Return(vis));
	
	Py_XDECREF(data);
	Py_XDECREF(vis);
	
	return visOut;
}

PyDoc_STRVAR(Simple_doc, \
"Parallel function to do simple fringing of complex data.");

/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef FringeMethods[] = {
	{"Simple",  (PyCFunction) Simple,  METH_VARARGS, Simple_doc}, 
	{NULL,      NULL,                  0,            NULL      }
};

PyDoc_STRVAR(fringe_doc, \
"Test extension to for fringing quickly.");

/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC initfringe(void) {
	PyObject *m;

	// Module definitions and functions
	m = Py_InitModule3("fringe", FringeMethods, fringe_doc);
	import_array();
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.1"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev: 1 $"));
	
}
