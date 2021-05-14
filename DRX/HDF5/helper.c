#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#ifdef _OPENMP
	#include <omp.h>
	
	// OpenMP scheduling method
	#ifndef OMP_SCHEDULER
	#define OMP_SCHEDULER dynamic
	#endif
#endif

#include "numpy/arrayobject.h"

#include "py3_compat.h"


#define IS_NOT_NAN(X) (X == X)


static PyObject *FastAxis1Mean(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long i, j, k, ik, nStand, nSamps, nChans, jPrime;
	
	if(!PyArg_ParseTuple(args, "O", &spectra)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChans;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Pointers
	float *a, *b;
	double tempV;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataF);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, jk, tempV, jPrime)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(ik=0; ik<nStand*nChans; ik++) {
			i = ik / nChans;
			k = ik % nChans;
			
			tempV = 0.0;
			
			jPrime = 0;
			for(j=0; j<nSamps; j++) {
				if( IS_NOT_NAN(*(a + nSamps*nChans*i + nChans*j + k)) ) {
					tempV += (double) *(a + nSamps*nChans*i + nChans*j + k);
					jPrime++;
				}
			}
			
			if( jPrime > 0 ) {
				tempV /= (double) jPrime;
			} else {
				tempV = 1.0;
			}
			*(b + nChans*i + k) = (float) tempV;
		}
	}
	
	Py_END_ALLOW_THREADS
	
	spectraF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return spectraF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return NULL;
}

PyDoc_STRVAR(FastAxis1Mean_doc, \
"Given a 3-D numpy.float32 array, compute the mean along the first axis\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (stands by time by channels) array of data\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * mean: 2-D numpy.float32 (stands by channels) of time-averaged spectra\n\
");


static PyObject *FastAxis0MinMax(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long i, j, k, nStand, nSamps, nChans;
	long int chanMin, chanMax;
	chanMin = 0;
	chanMax = -1;
	
	static char *kwlist[] = {"spectra", "chanMin", "chanMax", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|ll", kwlist, &spectra, &chanMin, &chanMax)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	if( chanMax < chanMin ) {
		chanMax = nChans;
	}
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) 2;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Pointers
	float *a, *b;
	float tempMin, tempMax;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataF);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, tempMin, tempMax)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(i=0; i<nStand; i++) {
			tempMin = 1e200;
			tempMax = -tempMin;
			
			for(j=0; j<nSamps; j++) {
				for(k=chanMin; k<chanMax; k++) {
					if( (float) *(a + nSamps*nChans*i + nChans*j + k) < tempMin ) {
						tempMin = (float) *(a + nSamps*nChans*i + nChans*j + k);
					} else if( (float) *(a + nSamps*nChans*i + nChans*j + k) > tempMax ) {
						tempMax = (float) *(a + nSamps*nChans*i + nChans*j + k);
					}
				}
			}
			
			*(b + 2*i + 0) = (float) tempMin;
			*(b + 2*i + 1) = (float) tempMax;
		}
	}
	
	Py_END_ALLOW_THREADS
	
	spectraF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return spectraF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return NULL;
}

PyDoc_STRVAR(FastAxis0MinMax_doc, \
"Given a 3-D numpy.float32 array, compute the minimum and maximum values along\n\
the zeroth axis\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (stands by time by channels) array of data\n\
\n\
Input keywords are:\n\
 * chanMin: Channel to start at (default = 0)\n\
 * chanMax: Channel to stop at (default = -1 => maximum channel number\n\
\n\
Outputs:\n\
 * minmax: 2-D numpy.float32 (stands by min/max) of the spectra\n\
");


static PyObject *FastAxis1Bandpass(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *bandpass, *spectraF;
	PyArrayObject *data=NULL, *dataB=NULL;
	
	long i, j, k, ik, nStand, nSamps, nChans;
	
	if(!PyArg_ParseTuple(args, "OO", &spectra, &bandpass)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	dataB = (PyArrayObject *) PyArray_ContiguousFromObject(bandpass, NPY_FLOAT32, 2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		goto fail;
	}
	if( dataB == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input bandpass array to 2-D float32");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	
	Py_BEGIN_ALLOW_THREADS
	
	// Pointers
	float *a, *b;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataB);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(ik=0; ik<nStand*nChans; ik++) {
			i = ik / nChans;
			k = ik % nChans;
			
			for(j=0; j<nSamps; j++) {
				*(a + nSamps*nChans*i + nChans*j + k) /= *(b + nChans*i + k);
			}
		}
	}
	
	Py_END_ALLOW_THREADS
	
	spectraF = Py_BuildValue("i", 1);
	
	Py_XDECREF(data);
	Py_XDECREF(dataB);
	
	return spectraF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(dataB);
	
	return NULL;
}

PyDoc_STRVAR(FastAxis1Bandpass_doc, \
"Given a 3-D numpy.float32 array of spectra and a 2-D numpy.float32 of\n\
bandpass shapes, apply a bandpass correction to the spectra.\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (stands by time by channels) array of data\n\
 * bandpass: 2-D numpy.float32 (stands by channels) array of bandpass shapes\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * bandpassed: 3-D numpy.float32 (stands by time by channels) of bandpass-\n\
               corrected spectra\n\
");


int cmpfloat(const void *a, const void *b) {
	/*
	 * cmpfloat - Comparison function for qsort-ing an array of float values.
	 */
	if( *(const float *) a < *(const float *) b ) {
		return -1;
	}
	return *(const float *) a > *(const float *) b;
}


static PyObject *FastAxis1Median(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long i, j, k, nStand, nSamps, nChans, jPrime;
	
	if(!PyArg_ParseTuple(args, "O", &spectra)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChans;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Pointers
	float *a, *b;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataF);
	
	float *tempV;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, tempV, jPrime)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(k=0; k<nChans; k++) {
			tempV = (float *) malloc(nSamps*sizeof(float));
			
			for(i=0; i<nStand; i++) {
				
				jPrime = 0;
				for(j=0; j<nSamps; j++) {
					*(tempV + jPrime) = *(a + nSamps*nChans*i + nChans*j + k);
					if( IS_NOT_NAN(*(tempV + jPrime)) ) {
						jPrime++;
					}
				}
				
				if( jPrime > 0 ) {
					qsort(tempV, jPrime, sizeof(float), cmpfloat);
					
					*(b + nChans*i + k) = *(tempV + jPrime/2);
				} else {
					*(b + nChans*i + k) = 1.0;
				}
			}
			
			free(tempV);
		}
	}
	
	Py_END_ALLOW_THREADS
	
	spectraF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return spectraF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return NULL;
}

PyDoc_STRVAR(FastAxis1Median_doc, \
"Given a 3-D numpy.float32 array, compute the approximate median along the\n\
first axis\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (stands by time by channels) array of data\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * median: 2-D numpy.float32 (stands by channels) of the median spectra\n\
\n\
.. note::\n\
\tThis function uses a median-of-medians method and may give slightly\n\
\tdifferent results relative to numpy.median() for arrays with even\n\
\tnumbers of elements.\n\
");


static PyObject *FastAxis0Percentiles5And99(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data=NULL;
	
	long i, j, k, nStand, nSamps, nChans, jPrime;
	long int stand, chanMin, chanMax;
	stand = 0;
	chanMin = 0;
	chanMax = -1;
	
	static char *kwlist[] = {"spectra", "stand", "chanMin", "chanMax", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "Ol|ll", kwlist, &spectra, &stand, &chanMin, &chanMax)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	if( chanMax < chanMin ) {
		chanMax = nChans;
	}
	
	// Pointers
	float *temp5, *temp99;
	
	Py_BEGIN_ALLOW_THREADS
	
	// Pointers
	float *a, *tempV;
	a = (float *) PyArray_DATA(data);
	
	temp5  = (float *) malloc((chanMax-chanMin)*sizeof(float));
	temp99 = (float *) malloc((chanMax-chanMin)*sizeof(float));
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, tempV, jPrime)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(k=chanMin; k<chanMax; k++) {
			i = stand;
			
			tempV = (float *) malloc(nSamps*sizeof(float));
			
			jPrime = 0;
			for(j=0; j<nSamps; j++) {
				*(tempV + jPrime) = *(a + nSamps*nChans*i + nChans*j + k);
				if( IS_NOT_NAN(*(tempV + jPrime)) ){
					jPrime++;
				}
			}
			
			qsort(tempV, jPrime, sizeof(float), cmpfloat);
			
			if( jPrime > 0 ) {
				*(temp5  + k - chanMin) = *(tempV + (int) (jPrime * 0.05));
				*(temp99 + k - chanMin) = *(tempV + (int) (jPrime * 0.99));
			} else {
				*(temp5  + k - chanMin) = 0.1;
				*(temp99 + k - chanMin) = 0.2;
			}
			
			free(tempV);
		}
	}
	
	qsort(temp5,  (chanMax-chanMin), sizeof(float), cmpfloat);
	qsort(temp99, (chanMax-chanMin), sizeof(float), cmpfloat);
	
	Py_END_ALLOW_THREADS
	
	spectraF = Py_BuildValue("ff", *(temp5 + (int) ((chanMax-chanMin)*0.05)), *(temp99 + (int) ((chanMax-chanMin) * 0.99)));
	
	free(temp5);
	free(temp99);
	Py_XDECREF(data);
	
	return spectraF;
	
fail:
	Py_XDECREF(data);
	
	return NULL;
}

PyDoc_STRVAR(FastAxis0Percentiles5And99_doc, \
"Given a 3-D numpy.float32 array, compute the approximate fifth and 99th \n\
percentiles for a particular index along the zeroth axis.\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (stands by time by channels) array of data\n\
 * stand: index along the zeroth axis to compute these values\n\
\n\
Input keywords are:\n\
 * chanMin: Channel to start at (default = 0)\n\
 * chanMax: Channel to stop at (default = -1 => maximum channel number\n\
\n\
Outputs:\n\
 * percentiles: two-element tuple of the fifth and 99th percentile values\n\
\n\
.. note::\n\
\tThis function uses percentile-of-percentile method to compute\n\
\tthe returned values.  These should be *close* to those returned\n\
\tby the scipy.stats.scoreatpercentile function for many cases.\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef HelperMethods[] = {
	{"FastAxis1Mean",              (PyCFunction) FastAxis1Mean,              METH_VARARGS,               FastAxis1Mean_doc}, 
	{"FastAxis0MinMax",            (PyCFunction) FastAxis0MinMax,            METH_VARARGS|METH_KEYWORDS, FastAxis0MinMax_doc}, 
	{"FastAxis1Bandpass",          (PyCFunction) FastAxis1Bandpass,          METH_VARARGS,               FastAxis1Bandpass_doc}, 
	{"FastAxis1Median",            (PyCFunction) FastAxis1Median,            METH_VARARGS,               FastAxis1Median_doc}, 
	{"FastAxis0Percentiles5And99", (PyCFunction) FastAxis0Percentiles5And99, METH_VARARGS|METH_KEYWORDS, FastAxis0Percentiles5And99_doc},
	{NULL,                         NULL,                                     0,                          NULL}
};

PyDoc_STRVAR(helper_doc, \
"HELPER - The Heuristic Extension aLlowing Parallel Efficiency in Rendering\n\
\n\
Extension to help plotHDF.py along by running the data computations in\n\
parallel.  See the individual functions for more details.");


/*
  Module Setup - Initialization
*/

MOD_INIT(_helper) {
        PyObject *m, *all;

        // Module definitions and functions
        MOD_DEF(m, "_helper", HelperMethods, helper_doc);
        if( m == NULL ) {
        return MOD_ERROR_VAL;
    }
        import_array();
        
        // Version information
        PyModule_AddObject(m, "__version__", PyString_FromString("0.2"));
        
        // Function listings
        all = PyList_New(0);
        PyList_Append(all, PyString_FromString("FastAxis1Mean"));
        PyList_Append(all, PyString_FromString("FastAxis0MinMax"));
        PyList_Append(all, PyString_FromString("FastAxis1Bandpass"));
        PyList_Append(all, PyString_FromString("FastAxis1Median"));
        PyList_Append(all, PyString_FromString("FastAxis0Percentiles5And99"));
        PyModule_AddObject(m, "__all__", all);
        
        return MOD_SUCCESS_VAL(m);
}
