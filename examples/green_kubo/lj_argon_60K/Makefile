all: gk

gk: trajectory.nc
	vibes output gk --interpolate $<

trajectory.nc: md/trajectory.son
	vibes output md -fc FORCE_CONSTANTS md/trajectory.son

md/trajectory.son:
	vibes run md md.in

clean:
	rm -r md
	rm trajectory.nc
	rm greenkubo.nc
