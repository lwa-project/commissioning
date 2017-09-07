CFLAGS = $(shell python-config --cflags) \
         $(shell python -c "import numpy; print '-I' + numpy.get_include()") \
         -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

LDFLAGS = $(shell python-config --ldflags)

_helper.so: helper.o
	$(CC) -o _helper.so helper.o $(LDFLAGS) -lm -shared -fopenmp -march=native
	
helper.o: helper.c
	$(CC) -c $(CFLAGS) -fPIC -o helper.o helper.c -funroll-loops -O3 -fopenmp -march=native
	
clean:
	rm -rf helper.o _helper.so
