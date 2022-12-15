PYTHON ?= python3
ifeq (, $(shell which $(PYTHON)))
	PYTHON = python
endif

.PHONY: all
all: fringe.so

fringe.so: fringe.c setup.py
	$(PYTHON) setup.py build
	mv build/lib*/*.so .
	rm -rf build
	
clean:
	rm -rf fringe.*so fringe.*dylib
