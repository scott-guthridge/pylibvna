.PHONY: examples

all: build
	( cd build && meson compile )

build:
	meson setup build

examples:
	( cd examples && ${MAKE} )

install:
	pip install .

fast_install:
	python3 -m pip install --no-build-isolation --editable .

test: all
	( cd build && meson test )
