#!/bin/bash
set -e -x

source regn.sh

if [ ! -d data ]
then
    mkdir data
fi

cp job.sh data/
cp job_vann.pbs data/job_vann.pbs

source system.sh

tid=1000
sluttid=100
ant=$(py "${tid}/${dt}")
sluttant=$(py "${sluttid}/${dt}")

dumpperiode=10000
langdumpperiode=$(py "${dumpperiode}*10")

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

    thermo 100
    thermo_style custom step time temp ke pe etotal press pzz spcpu cpuremain
    dump lagring all custom ${dumpperiode} vann.in.bin id type x y z
    dump sjeldenlagring all custom ${langdumpperiode} vannlavfrekv.in.bin id type x y z

    minimize 1e-6 1e-6 1000 1000

    write_data 01_medvann_minimert.data

    fix termostat kanbevegeseg nvt temp $T $T 1.0
    velocity kanbevegeseg create 600 277385 mom yes loop geom

    group silisium type 1
    group oksygen type 2

    timestep ${dt}

    restart 100000 vann.*.restart

    run ${ant}

    dump hyppiglagring all custom 1000 vannhoyfrekv.in.bin id type x y z
    dump veldighyppiglagring all custom 100 vannveldighoyfrekv.in.bin id type x y z

    run ${sluttant}

    write_data 02_medvann_1.1ns.data

    " > data/vann.in
