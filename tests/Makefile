test: test_coverage clean

test_basic:
	pytest --ignore fireworks

test_coverage:
	pytest --cov --cov-config=../setup.cfg --ignore fireworks

test_all:
	pytest

clean:
	./clean.py
