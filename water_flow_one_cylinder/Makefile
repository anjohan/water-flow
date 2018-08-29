SHELL := /bin/bash
indir = .
outdir = data
lammps = $(shell find ~ -name lmp_*mpi* 2> /dev/null) -var INDIR $(indir) -var OUTDIR $(outdir)

$(shell mkdir -p $(outdir))

all: $(outdir)/done

$(outdir)/data.amorphous: $(indir)/in.make_amorphous $(indir)/in.common_variables
	mpirun $(lammps) -in $<

$(outdir)/data.amorphous_withH: $(indir)/data_addH.py $(outdir)/data.amorphous
	python $^ $@

$(outdir)/restart.SiO2: $(indir)/in.setup $(outdir)/data.amorphous_withH
	mpirun $(lammps) -in $<

$(outdir)/add_water.inp: $(indir)/in.write_packmol $(indir)/in.common_variables
	mpirun $(lammps) -in $<

$(outdir)/xyz.water: $(outdir)/add_water.inp
	packmol < $<

$(outdir)/dump.water: $(outdir)/xyz.water $(indir)/xyz2dump.py $(outdir)/restart.SiO2
	python $(indir)/xyz2dump.py $$(cat $(outdir)/num_atoms.dat) $< $@

$(outdir)/restart.water: $(indir)/in.add_water $(outdir)/dump.water $(outdir)/restart.SiO2
	mpirun $(lammps) -in $<

$(outdir)/done: $(indir)/in.flow $(outdir)/restart.water
	mpirun $(lammps) -in $<
