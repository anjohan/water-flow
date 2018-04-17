data:
	mkdir -p data

lammps = $(shell find ~ -name lmp_*mpi* 2> /dev/null)

log.system: system.in
	$(lammps) -in system.in

data/amorphous.data: data make_amorphous.in system.in
	mpirun $(lammps) -in make_amorphous.in
data/unpassivated.data: data/amorphous.data data system.in common_regions_groups.in make_setup.in data
	mpirun $(lammps) -in make_setup.in
data/unpassivatedwithH.data: data/unpassivated.data fix_setup_addH.py data
	python fix_setup_addH.py
data/passivated.xyz: data/unpassivated.data system.in common_regions_groups.in passivate.in data
	mpirun -np 1 $(lammps) -in passivate.in
water.inp: packmol_template.in water.py log.system
	python water.py
data/withwater.xyz: data/passivated.xyz water.inp water.xyz
	packmol < water.inp
data/withwater.data: data/withwater.xyz withwater_xyztodata.py log.system
	python withwater_xyztodata.py
