all:
	python3 -m build

install: all
	pip install dist/libvna*.whl

fast_install:
	python3 -m pip install --no-build-isolation --editable .

test:
	python -m unittest discover tests
