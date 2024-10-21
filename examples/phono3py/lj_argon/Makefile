nq=6

all: run output

run:
	vibes run phono3py

output:
	vibes output phono3py phono3py/trajectory.son
	cd phono3py/output ; phono3py-load phono3py_params.yaml --mesh $(nq) $(nq) $(nq) --br

clean:
	rm -r phono3py

