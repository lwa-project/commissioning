#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "numpy/arrayobject.h"


static PyObject *FastAxis0Mean(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data, *dataF;
	
	long i, j, k, jk, nStand, nSamps, nChans, iPrime;
	
	if(!PyArg_ParseTuple(args, "O", &spectra)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	
	// Get the properties of the data
	nSamps = (long) data->dimensions[0];
	nStand = (long) data->dimensions[1];
	nChans = (long) data->dimensions[2];
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChans;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_FLOAT32);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Pointers
	float *a, *b;
	float tempV;
	a = (float *) data->data;
	b = (float *) dataF->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, jk, tempV, iPrime)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(jk=0; jk<nStand*nChans; jk++) {
			j = jk / nChans;
			k = jk % nChans;
			
			tempV = 0.0;
			
			iPrime = 0;
			for(i=0; i<nSamps; i++) {
				if( *(a + nStand*nChans*i + nChans*j + k) == *(a + nStand*nChans*i + nChans*j + k) ) {
					tempV += (float) *(a + nStand*nChans*i + nChans*j + k);
					iPrime++;
				}
			}
			
			if( iPrime > 0 ) {
				tempV /= (float) iPrime;
			} else {
				tempV = 1.0;
			}
			*(b + nChans*j + k) = tempV;
		}
	}
	
	Py_XDECREF(data);
	
	spectraF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);
	
	return spectraF;
}

PyDoc_STRVAR(FastAxis0Mean_doc, \
"Given a 3-D numpy.float32 array, compute the mean along the zeroth axis\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * mean: 2-D numpy.float32 (stands by channels) of time-averaged spectra\n\
");


static PyObject *FastAxis1MinMax(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data, *dataF;
	
	long i, j, k, nStand, nSamps, nChans;
	long int chanMin, chanMax;
	chanMin = 0;
	chanMax = -1;
	
	static char *kwlist[] = {"spectra", "chanMin", "chanMax", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|ll", kwlist, &spectra, &chanMin, &chanMax)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	
	// Get the properties of the data
	nSamps = (long) data->dimensions[0];
	nStand = (long) data->dimensions[1];
	nChans = (long) data->dimensions[2];
	if( chanMax < chanMin ) {
		chanMax = nChans;
	}
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) 2;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_FLOAT32);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Pointers
	float *a, *b;
	float tempMin, tempMax;
	a = (float *) data->data;
	b = (float *) dataF->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, tempMin, tempMax)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(j=0; j<nStand; j++) {
			tempMin = 1e200;
			tempMax = -tempMin;
			
			for(k=chanMin; k<chanMax; k++) {
				for(i=0; i<nSamps; i++) {
					if( (float) *(a + nStand*nChans*i + nChans*j + k) < tempMin ) {
						tempMin = (float) *(a + nStand*nChans*i + nChans*j + k);
					} else if( (float) *(a + nStand*nChans*i + nChans*j + k) > tempMax ) {
						tempMax = (float) *(a + nStand*nChans*i + nChans*j + k);
					}
				}
			}
			
			*(b + 2*j + 0) = (float) tempMin;
			*(b + 2*j + 1) = (float) tempMax;
		}
	}
	
	Py_XDECREF(data);
	
	spectraF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);
	
	return spectraF;
}

PyDoc_STRVAR(FastAxis1MinMax_doc, \
"Given a 3-D numpy.float32 array, compute the minimum and maximum values along\n\
the first axis\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
\n\
Input keywords are:\n\
 * chanMin: Channel to start at (default = 0)\n\
 * chanMax: Channel to stop at (default = -1 => maximum channel number\n\
\n\
Outputs:\n\
 * minmax: 2-D numpy.float32 (stands by min/max) of the spectra\n\
");


static PyObject *FastAxis0Bandpass(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *bandpass, *spectraF;
	PyArrayObject *data, *dataB;
	
	long i, j, k, jk, nStand, nSamps, nChans;
	
	if(!PyArg_ParseTuple(args, "OO", &spectra, &bandpass)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	dataB = (PyArrayObject *) PyArray_ContiguousFromObject(bandpass, NPY_FLOAT32, 2, 2);
	// Get the properties of the data
	nSamps = (long) data->dimensions[0];
	nStand = (long) data->dimensions[1];
	nChans = (long) data->dimensions[2];
	
	// Pointers
	float *a, *b;
	a = (float *) data->data;
	b = (float *) dataB->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(jk=0; jk<nStand*nChans; jk++) {
			j = jk / nChans;
			k = jk % nChans;
			
			for(i=0; i<nSamps; i++) {
				*(a + nStand*nChans*i + nChans*j + k) /= *(b + nChans*j + k);
			}
		}
	}
	
	Py_XDECREF(data);
	Py_XDECREF(dataB);
	
	spectraF = Py_BuildValue("i", 1);
	return spectraF;
}

PyDoc_STRVAR(FastAxis0Bandpass_doc, \
"Given a 3-D numpy.float32 array of spectra and a 2-D numpy.float32 of\n\
bandpass shapes, apply a bandpass correction to the spectra.\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
 * bandpass: 2-D numpy.float32 (stands by channels) array of bandpass shapes\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * bandpassed: 3-D numpy.float32 (time by stands by channels) of bandpass-\n\
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


static PyObject *FastAxis0Median(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data, *dataF;
	
	long i, j, k, nStand, nSamps, nChans, iPrime;
	
	if(!PyArg_ParseTuple(args, "O", &spectra)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	
	// Get the properties of the data
	nSamps = (long) data->dimensions[0];
	nStand = (long) data->dimensions[1];
	nChans = (long) data->dimensions[2];
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChans;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_FLOAT32);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Pointers
	float *a, *b;
	a = (float *) data->data;
	b = (float *) dataF->data;
	
	float *tempV;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, tempV, iPrime)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(k=0; k<nChans; k++) {
			tempV = (float *) malloc(nSamps*sizeof(float));
			
			for(j=0; j<nStand; j++) {
				
				iPrime = 0;
				for(i=0; i<nSamps; i++) {
					*(tempV + iPrime) = *(a + nStand*nChans*i + nChans*j + k);
					if( *(tempV + iPrime) == *(tempV + iPrime) ) {
						iPrime++;
					}
				}
				
				if( iPrime > 0 ) {
					qsort(tempV, iPrime, sizeof(float), cmpfloat);
					
					*(b + nChans*j + k) = *(tempV + iPrime/2);
				} else {
					*(b + nChans*j + k) = 1.0;
				}
			}
			
			free(tempV);
		}
	}
	
	Py_XDECREF(data);
	
	spectraF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);
	
	return spectraF;
}

PyDoc_STRVAR(FastAxis0Median_doc, \
"Given a 3-D numpy.float32 array, compute the approximate median along the\n\
zeroth axis\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
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


static PyObject *FastAxis1Percentiles5And99(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data;
	
	long i, j, k, nStand, nSamps, nChans, iPrime;
	long int stand, chanMin, chanMax;
	stand = 0;
	chanMin = 0;
	chanMax = -1;
	
	static char *kwlist[] = {"spectra", "stand", "chanMin", "chanMax", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "Ol|ll", kwlist, &spectra, &stand, &chanMin, &chanMax)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	
	// Get the properties of the data
	nSamps = (long) data->dimensions[0];
	nStand = (long) data->dimensions[1];
	nChans = (long) data->dimensions[2];
	if( chanMax < chanMin ) {
		chanMax = nChans;
	}
	
	// Pointers
	float *a;
	a = (float *) data->data;
	
	float *tempV, *temp5, *temp99;
	
	temp5  = (float *) malloc((chanMax-chanMin)*sizeof(float));
	temp99 = (float *) malloc((chanMax-chanMin)*sizeof(float));
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, tempV, iPrime)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(k=chanMin; k<chanMax; k++) {
			j = stand;
			
			tempV = (float *) malloc(nSamps*sizeof(float));
			
			iPrime = 0;
			for(i=0; i<nSamps; i++) {
				*(tempV + iPrime) = *(a + nStand*nChans*i + nChans*j + k);
				if( *(tempV + iPrime) == *(tempV + iPrime) ){
					iPrime++;
				}
			}
			
			qsort(tempV, iPrime, sizeof(float), cmpfloat);
			
			if( iPrime > 0 ) {
				*(temp5  + k - chanMin) = *(tempV + (int) (iPrime * 0.05));
				*(temp99 + k - chanMin) = *(tempV + (int) (iPrime * 0.99));
			} else {
				*(temp5  + k - chanMin) = 0.1;
				*(temp99 + k - chanMin) = 0.2;
			}
			
			free(tempV);
		}
	}
	
	Py_XDECREF(data);
	
	qsort(temp5,  (chanMax-chanMin), sizeof(float), cmpfloat);
	qsort(temp99, (chanMax-chanMin), sizeof(float), cmpfloat);
	
	spectraF = Py_BuildValue("ff", *(temp5 + (int) ((chanMax-chanMin)*0.05)), *(temp99 + (int) ((chanMax-chanMin) * 0.99)));
	
	free(temp5);
	free(temp99);
	
	return spectraF;
}

PyDoc_STRVAR(FastAxis1Percentiles5And99_doc, \
"Given a 3-D numpy.float32 array, compute the approximate fifth and 99th \n\
percentiles for a particular index along the first axis.\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
 * stand: index along the first axis to compute these values\n\
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
	{"FastAxis0Mean",              (PyCFunction) FastAxis0Mean,              METH_VARARGS,               FastAxis0Mean_doc}, 
	{"FastAxis1MinMax",            (PyCFunction) FastAxis1MinMax,            METH_VARARGS|METH_KEYWORDS, FastAxis1MinMax_doc}, 
	{"FastAxis0Bandpass",          (PyCFunction) FastAxis0Bandpass,          METH_VARARGS,               FastAxis0Bandpass_doc}, 
	{"FastAxis0Median",            (PyCFunction) FastAxis0Median,            METH_VARARGS,               FastAxis0Median_doc}, 
	{"FastAxis1Percentiles5And99", (PyCFunction) FastAxis1Percentiles5And99, METH_VARARGS|METH_KEYWORDS, FastAxis1Percentiles5And99_doc},
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

PyMODINIT_FUNC init_helper(void) {
	PyObject *m;

	// Module definitions and functions
	m = Py_InitModule3("_helper", HelperMethods, helper_doc);
	import_array();
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.1"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
}
