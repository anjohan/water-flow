#!/bin/bash
set -e -x

source regn.sh

if [ ! -d data ]
then
    mkdir data
fi

cp job.sh data/job.sh

source system.sh

echo $nederst > nederst.dat
echo $x > x.dat
echo $y > y.dat

echo "log log.passivering
    #/atom 1 silicon
    #/atom 2 oxygen
    #/atom 3 hydrogen
    #/bond 1 2 2.6
    #/bond 2 3 1.4

    units metal
    boundary p p p
    atom_style	molecular

    read_data ../betacristobalite.data
    replicate $antx $anty 25

    pair_style  usc
    pair_coeff  * * ../SiOH2O.vashishta Si O H
    pair_modify coord 2 1 2.0 0.3
    pair_modify coord 2 3 1.4 0.3
    mass            1 28.08
    mass            2 15.9994
    mass            3 1.00794

    region helesystemet block EDGE EDGE EDGE EDGE EDGE EDGE
    region ytreblokk block EDGE EDGE EDGE EDGE ${midtenavbunn} ${midtenavtopp} side out
    region bunnhard block EDGE EDGE EDGE EDGE ${nederst} ${midtenavbunn}
    region boks2 block EDGE EDGE EDGE EDGE ${midtenavbunn} ${toppavbunn}
    region kule sphere ${xmid} ${ymid} ${zkule} ${R}
    region boks1 block EDGE EDGE EDGE EDGE ${bunnavtopp} ${toppavtopp}
    region slett_boks block EDGE EDGE EDGE EDGE ${toppavtopp} EDGE

    region R_topp union 2 kule boks1
    region ALT union 4 kule boks1 boks2 bunnhard side out
    region bunn union 2 bunnhard boks2
    region sylinder cylinder z ${xmid} ${xmid} ${sylinderradius} EDGE EDGE


    group G_topp region R_topp
    group hard_topp region boks1
    group frosten region bunnhard
    group G_bunn region bunn
    group kanbevegeseg subtract all frosten
    group kanbevegesegutentoppen subtract kanbevegeseg hard_topp

    velocity kanbevegesegutentoppen create $T 277385 mom yes loop geom
    delete_atoms region slett_boks
    delete_atoms region ALT

    compute mass_center G_topp com
    compute normalkraft hard_topp reduce sum fz

    thermo 10
    thermo_style custom step time temp press pzz etotal c_mass_center[1] c_mass_center[3] c_normalkraft cpuremain

    group sio2 region helesystemet

    fix passivate all passivate
    run 1
    unfix passivate

    group newAtoms subtract all sio2
    group ytreblokk region ytreblokk
    group OHytreblokk intersect newAtoms ytreblokk
    delete_atoms group OHytreblokk compress yes
    group kanbevegesegutentoppen union kanbevegesegutentoppen newAtoms
    fix konstant2 kanbevegesegutentoppen nvt temp $T $T 1.0

    region nedrehalvdel block EDGE EDGE EDGE EDGE EDGE ${midten}
    group nedrehalvdel region nedrehalvdel
    write_dump nedrehalvdel xyz nedrehalvdel.xyz
    delete_atoms group nedrehalvdel

    write_restart ovrehalvdel.restart
    " > data/passivering.in
