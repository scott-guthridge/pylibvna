all:
	python3 -m build

install: all
	pip install dist/libvna*.whl

fast_install:
	pip install -e . --user

test:
	python -m unittest discover tests
