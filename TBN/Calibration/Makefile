CFLAGS = $(shell python-config --cflags) \
         $(shell python -c "from __future__ import print_function; import numpy; print('-I' + numpy.get_include())") \
         -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

LDFLAGS = $(shell python-config --ldflags)

finge.so: fringe.o
	$(CC) -o fringe.so fringe.o $(LDFLAGS) -lm -shared -fopenmp -march=native

fringe.o: fringe.c
	$(CC) -c $(CFLAGS) -fPIC -o fringe.o fringe.c -funroll-loops -O3 -fopenmp -march=native

clean:
	rm -rf fringe.o fringe.so

