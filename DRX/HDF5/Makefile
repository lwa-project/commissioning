CFLAGS = $(shell python-config --cflags) $(shell python -c "import numpy; print '-I' + numpy.get_include()")

LDFLAGS = $(shell python-config --ldflags)

_helper.so: helper.o
	$(CC) $(LDFLAGS) -o _helper.so helper.o -lm -shared -fopenmp
	
helper.o:
	$(CC) -c $(CFLAGS) -fPIC -o helper.o helper.c -funroll-loops -O3 -fopenmp
	
clean:
	rm -rf helper.o _helper.so
