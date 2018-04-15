#!/bin/bash
set -e -x

source regn.sh

if [ ! -d data ]
then
    mkdir data
fi

cp job.sh data/job.sh
cp ./job_passivering.pbs data/job_passivering.pbs

source system.sh

echo "log log.passivering
    #/atom 1 silicon
    #/atom 2 oxygen
    #/atom 3 hydrogen
    #/bond 1 2 2.6
    #/bond 2 3 1.4

    units metal
    boundary p p p
    atom_style	molecular

    read_data klartilpassivering.data

    pair_style  usc
    pair_coeff  * * ../SiOH2O.vashishta Si O H
    pair_modify coord 2 1 2.0 0.3
    pair_modify coord 2 3 1.4 0.3
    mass            1 28.08
    mass            2 15.9994
    mass            3 1.00794

    thermo 10
    thermo_style custom step time temp press pzz etotal cpuremain

    fix passivate all passivate
    run 1
    unfix passivate

    write_restart passivert.restart
    write_dump all xyz passivert.xyz
    " > data/passivering.in
