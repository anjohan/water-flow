#!/bin/bash
set -e -x

source regn.sh

if [ ! -d data ]
then
    mkdir data
fi

cp job.sh data/
cp job_vann.pbs data/job_vann.pbs
cp job_fram_vann.sh data/job_fram_vann.sh

source system.sh

tid=3000
sluttid=100
ant=$(py "${tid}/${dt}")
sluttant=$(py "${sluttid}/${dt}")

dumpperiode=10000
langdumpperiode=$(py "${dumpperiode}*10")
dx=$(py "${x}/30")
Vperchunk=$(py "${dx}**3")

echo "log log.vann
    #/atom 1 silicon
    #/atom 2 oxygen
    #/atom 3 hydrogen
    #/bond 1 2 2.6
    #/bond 2 3 1.4

    units metal
    boundary p p p
    atom_style atomic

    read_data medvann.data

    pair_style  usc
    pair_coeff  * * ../SiOH2O.vashishta Si O H
    pair_modify coord 2 1 2.0 0.3
    pair_modify coord 2 3 1.4 0.3
    mass            1 28.08
    mass            2 15.9994
    mass            3 1.00794

    region helesystemet block EDGE EDGE EDGE EDGE EDGE EDGE

    region tomromnederst block EDGE EDGE EDGE EDGE EDGE ${nederst}
    region tomromz block EDGE EDGE EDGE EDGE ${nederst} EDGE side out
    region mykbunn1 block EDGE EDGE EDGE EDGE ${nederst} ${midtenavbunn1}
    region hardbunn block EDGE EDGE EDGE EDGE ${midtenavbunn1} ${midtenavbunn2}
    region mykbunn2 block EDGE EDGE EDGE EDGE ${midtenavbunn2} ${toppavbunn}
    region bunn union 3 hardbunn mykbunn1 mykbunn2
    region sylinder cylinder z ${xmid} ${xmid} ${sylinderradius} ${nederst} EDGE

    group bunn region bunn
    group hardbunn region hardbunn
    group sylinder region sylinder

    group kanbevegeseg subtract all hardbunn

    group vann subtract kanbevegeseg bunn sylinder

    timestep ${dt}
    minimize 1e-6 1e-6 1000 1000
    write_data 01_medvann_minimert.data

    compute potperatom all pe/atom
    compute inndeling all chunk/atom bin/3d x lower ${dx} y lower ${dx} z lower ${dx}

    compute stress_x_V_per_atom all stress/atom NULL ke pair
    variable trykk_per_chunk atom '-(c_stress_x_V_per_atom[1]+c_stress_x_V_per_atom[2]+c_stress_x_V_per_atom[3])/(3*${Vperchunk})'

    fix temp_per_chunk all ave/chunk 100 100 10000 inndeling temp file temp_per_chunk.profile
    fix trykk_per_chunk all ave/chunk 100 100 10000 inndeling v_trykk_per_chunk norm sample file trykk_per_chunk.profile

    fix termostat kanbevegeseg nvt temp $T $T 1.0
    velocity kanbevegeseg create 600 277385 mom yes loop geom

    compute vacf vann vacf
    thermo 100
    thermo_style custom step time temp ke pe etotal press pzz spcpu cpuremain c_vacf[4]
    dump lagring all custom ${dumpperiode} vann.in.bin id type x y z xu yu zu c_potperatom
    dump sjeldenlagring all custom ${langdumpperiode} vannlavfrekv.in.bin id type x y z xu yu zu c_potperatom

    restart 100000 vann.*.restart

    run ${ant}

    write_data 02_medvann_${tid}ps.data

    " > data/vann.in
