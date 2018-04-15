#!/bin/bash
set -e -x

source regn.sh

for i in $(cat i.txt)
do
    if [ ! -d data/ilik$i ]
    then
        mkdir data/ilik$i
    fi
    if [ ! -f data/passivert.restart ]
    then
        echo "Passiver først!"
        exit 1
    else
        cp data/passivert.restart data/ilik$i/
    fi
    cp job.sh data/ilik$i/job.sh

    source system.sh
    nedtid=100000
    indent=$(py "3*$a")
    dt=0.0005
    v=$(py "(${indent}+${kuleklaring})/(${nedtid}*${dt})")
    antvann=50000

    echo "
        log log.ned
        #/atom 1 silicon
        #/atom 2 oxygen
        #/atom 3 hydrogen
        #/bond 1 2 2.6
        #/bond 2 3 1.4

        units metal
        boundary p p p
        atom_style	molecular
        read_restart passivert.restart

        pair_style  usc
        pair_coeff  * * ../../SiOH2O.vashishta Si O H
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

        compute mass_center G_topp com
        compute normalkraft hard_topp reduce sum fz

        thermo 10
        thermo_style custom step time temp press pzz etotal c_mass_center[1] c_mass_center[3] c_normalkraft cpuremain

        region waterRegion union 5 kule boks1 boks2 bunnhard ytreblokk side out
        molecule water ../../water.mol
        group altunntattvann region helesystemet
        create_atoms 0 random ${antvann} 23452 waterRegion mol water 1353
        group vann subtract all altunntattvann
        group kanbevegesegutentoppen union kanbevegesegutentoppen vann
        minimize 1e-6 1e-6 100 100
        velocity kanbevegesegutentoppen create $T 277385 mom yes loop geom
        dump lagring all custom 2000 ned.in.bin id type x y z
        dump litenlagring all custom 40000 liten_ned.in.bin id type x y z
        fix nedover hard_topp move linear 0 0 -$v
        fix konstant2 kanbevegesegutentoppen nvt temp $T $T 1.0

        timestep ${dt}

        run ${nedtid}

        unfix konstant2
        fix konstant2 kanbevegesegutentoppen nvt temp 300 300 1.0

        unfix nedover
        fix nedover hard_topp move linear 0 0 0
        run 100000
        write_restart nede.restart
    " > data/ilik$i/ned.in
done
