NATOMS=100
T=60

all: md post

geometry.in.primitive:
	ln -sf ../relaxation_02/relaxation/geometry.in.next_step $@

geometry.in.supercell: geometry.in.primitive
	vibes utils make-supercell $< -n $(NATOMS) -o $@

geometry.in: geometry.in.supercell
	vibes utils create-samples $< -T $(T)
	mv geometry.in.supercell.* $@

md: geometry.in md.in.template
	cp md.in.template md.in
	sed -i -e 's/<temp>/$(T)/g' md.in
	vibes run md

post: md/trajectory.son
	vibes output md md/trajectory.son
	vibes info md trajectory.nc

clean:
	rm geometry.in* md.in trajectory.nc
	rm -r md
