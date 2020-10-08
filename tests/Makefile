test: test_coverage clean

test_basic:
	pytest --ignore fireworks

test_coverage:
	flake8 .
	flake8 ../vibes
	pytest --cov --cov-config=../setup.cfg --ignore fireworks

test_all:
	pytest

clean:
	./clean.py
