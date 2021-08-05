.PHONY: all
all: _helper.so

_helper.so: helper.c setup.py
	python setup.py build
	mv build/lib*/*.so .
	rm -rf build
	
clean:
	rm -rf _helper.so
