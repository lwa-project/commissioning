#include "Python.h"
#include <cmath>
#include <complex>

#ifdef _OPENMP
	#include <omp.h>
	
	// OpenMP scheduling method
	#ifndef OMP_SCHEDULER
	#define OMP_SCHEDULER dynamic
	#endif
#endif

#include "numpy/arrayobject.h"

/*
  Complex type definitions
*/

typedef std::complex<float> Complex32;
typedef std::complex<double> Complex64;

/*
  Simple Fringing of Complex data
*/

template<typename InType, typename OutType>
void compute_fringe(long nStand,
                    long nSamp,
									  float clipLevel,
										int refX,
										int refY,
									  InType const* data,
									  OutType* vis) {
	long i, j, k;
	int ref;
	InType temp, tempRef;
	Complex64 tempV;
	
	Py_BEGIN_ALLOW_THREADS
  
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(j, k, ref, temp, tempRef, tempV)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
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
			for(j=0; j<nSamp; j++) {
				// Load
				temp = *(data + (nSamp*i + j));
				tempRef = *(data + (nSamp*ref + j));
				
				// Bad input value
				if( abs(temp) >= clipLevel ) {
					continue;
				}
				
				// Bad reference input value
				if( abs(tempRef) >= clipLevel ) {
					continue;
				}
				
				tempV += temp * conj(tempRef);
				k++;
			}
			
			// Average
			if( k > 0 ) {
				*(vis + i) = tempV / (Complex64) k;
			} else {
				*(vis + i) = 0.0;
			}
		}
	}

	Py_END_ALLOW_THREADS
}

static PyObject *Simple(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *visOut;
	PyArrayObject *data=NULL, *vis=NULL;
	
	int refX, refY;
	float clipLevel;
	long nStand, nSamps;
	
	if(!PyArg_ParseTuple(args, "Oiif", &signals, &refX, &refY, &clipLevel)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, 
																											  PyArray_TYPE((PyArrayObject *) signals), 
																											  2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signals as a 2-D array");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[1];
	dims[0] = (npy_intp) nStand;
	vis = (PyArrayObject*) PyArray_EMPTY(1, dims, NPY_COMPLEX64, 0);
	if(vis == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
#define LAUNCH_FRINGE(IterType) \
        compute_fringe<IterType>(nStand, nSamps, clipLevel, refX, refY, \
					                       (IterType*) PyArray_DATA(data), \
                                 (Complex32*) PyArray_DATA(vis))
	switch( PyArray_TYPE(data) ){
		case( NPY_COMPLEX64  ): LAUNCH_FRINGE(Complex32); break;
		case( NPY_COMPLEX128 ): LAUNCH_FRINGE(Complex64); break;
		default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
}
		
#undef LAUNCH_FRINGE
	
	visOut = Py_BuildValue("O", PyArray_Return(vis));
	
	Py_XDECREF(data);
	Py_XDECREF(vis);
	
	return visOut;

fail:
	Py_XDECREF(data);
	Py_XDECREF(vis);
	
	return NULL;
}

PyDoc_STRVAR(Simple_doc, \
"Parallel function to do simple fringing of complex data.");

/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef fringe_methods[] = {
	{"Simple",  (PyCFunction) Simple,  METH_VARARGS, Simple_doc}, 
	{NULL,      NULL,                  0,            NULL      }
};

PyDoc_STRVAR(fringe_doc, \
"Test extension to for fringing quickly.");

/*
  Module Setup - Initialization
*/

static int fringe_exec(PyObject* module) {
		import_array();
		
		// Version and revision information
		PyModule_AddObject(module, "__version__", PyUnicode_FromString("0.1"));
		
		// Function listings
		PyObject* all = PyList_New(0);
		PyList_Append(all, PyUnicode_FromString("Simple"));
		PyModule_AddObject(module, "__all__", all);
		return 0;
}

static PyModuleDef_Slot fringe_slots[] = {
    {Py_mod_exec, (void *)&fringe_exec},
    {0,           NULL}
};

static PyModuleDef fringe_def = {
    PyModuleDef_HEAD_INIT,    /* m_base */
    "fringe",                 /* m_name */
    fringe_doc,               /* m_doc */
    0,                        /* m_size */
    fringe_methods,           /* m_methods */
    fringe_slots,             /* m_slots */
    NULL,                     /* m_traverse */
    NULL,                     /* m_clear */
    NULL,                     /* m_free */
};

PyMODINIT_FUNC PyInit_fringe(void) {
	  return PyModuleDef_Init(&fringe_def);
}
