finge.so: fringe.o
	gcc fringe.o -o fringe.so -fopenmp -lm -shared

fringe.o:
	gcc -c fringe.c -o fringe.o -I/usr/include/python2.6 -I/usr/lib/python2.6/dist-packages/numpy/core/include -fopenmp -O3 -fPIC -shared

