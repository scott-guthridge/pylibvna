all: install

slow_build:
	python3 -m build

install:
	pip install --user -e .

test:
	pytest tests
