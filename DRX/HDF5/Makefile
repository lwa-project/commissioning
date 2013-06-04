CFLAGS = $(shell python-config --includes) $(shell python -c "import numpy; print '-I' + numpy.__path__[0] + '/core/include/numpy/'")

_helper.so: helper.o
	gcc $(CFLAGS) -o _helper.so helper.o -lm -shared -fopenmp -lfftw3
	
helper.o:
	gcc -c $(CFLAGS) -fPIC -o helper.o helper.c -funroll-loops -O3 -fopenmp
	
clean:
	rm -rf helper.o _helper.so
