all:
	pip wheel -e .

slow_build:
	python3 -m build


install:
	pip install -e . --user

sdist:
	python3 setup.py sdist

test:
	pytest tests
