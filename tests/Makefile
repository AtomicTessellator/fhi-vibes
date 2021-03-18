test: test_coverage clean

test_basic:
	flake8 && flake8 ..
	pytest --ignore aims --ignore fireworks

test_coverage:
	flake8 && flake8 ..
	pytest --cov --cov-config=../setup.cfg --ignore fireworks -x

test_all:
	flake8 && flake8 ..
	pytest

clean:
	./clean.py
