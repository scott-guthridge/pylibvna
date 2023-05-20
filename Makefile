all:
	python3 -m build

install: all
	pip install dist/vna*.whl

fast_install:
	pip install -e . --user

test:
	pytest tests
