lammps = $(shell find ~ -name lmp 2> /dev/null)

# all: data/water_0.0001.bin
all: $(foreach force,$(shell seq 0.00001 0.00001 0.0001), data/water_$(force).bin)

log.system: system.in
	mpirun $(lammps) -in system.in

data/amorphous.data: make_amorphous.in system.in
	mkdir -p data
	mpirun $(lammps) -in make_amorphous.in  # 30 min on 128 cores

data/unpassivated.data: data/amorphous.data system.in common_regions_groups.in make_setup.in
	mpirun $(lammps) -in make_setup.in

data/unpassivatedwithH.data: data/unpassivated.data fix_setup_addH.py
	python3 fix_setup_addH.py

data/passivated.xyz: data/unpassivatedwithH.data system.in common_regions_groups.in passivate.in
	mpirun -np 1 $(lammps) -in passivate.in

data/water.inp: packmol_template.inp water.py data/unpassivatedwithH.data
	python3 water.py

data/withwater.xyz: data/passivated.xyz data/water.inp water.xyz
	packmol < data/water.inp

data/withwater.data: data/withwater.xyz withwater_xyztodata.py
	python3 withwater_xyztodata.py

data/water_%.bin: data/withwater.data water.in system.in common_regions_groups.in
	mpirun $(lammps) -var FORCE $* -in water.in
