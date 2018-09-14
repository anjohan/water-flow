lammps = $(shell find ~ -name lmp_*mpi* 2> /dev/null)

log.system: system.in
	mpirun $(lammps) -in system.in

data/amorphous.data: make_amorphous.in system.in
	mkdir -p data
	date
	mpirun $(lammps) -in make_amorphous.in  # 30 min on 128 cores
	date
data/unpassivated.data: data/amorphous.data system.in common_regions_groups.in make_setup.in
	date
	mpirun $(lammps) -in make_setup.in
	date
data/unpassivatedwithH.data: data/unpassivated.data fix_setup_addH.py
	date
	python fix_setup_addH.py
	date
data/passivated.xyz: data/unpassivatedwithH.data system.in common_regions_groups.in passivate.in
	date
	mpirun -np 1 $(lammps) -in passivate.in
	date
data/water.inp: packmol_template.inp water.py data/unpassivatedwithH.data
	date
	python water.py
	date
data/withwater.xyz: data/passivated.xyz data/water.inp water.xyz
	date
	packmol < data/water.inp
	date
data/withwater.data: data/withwater.xyz withwater_xyztodata.py
	date
	python3 withwater_xyztodata.py
	date
data/water_%.bin: data/withwater.data water.in system.in common_regions_groups.in
	date
	mpirun $(lammps) -var FORCE $* -in water.in
	date
