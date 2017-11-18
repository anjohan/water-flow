#!/bin/bash
set -e -x

source regn.sh

if [ ! -d data ]
then
    mkdir data
fi

cp job_oppsett.pbs data/job_oppsett.pbs

source system.sh

echo $toppavbunn > toppavbunn.dat
echo $nederst > nederst.dat
echo $x > x.dat
echo $y > y.dat
echo $toppavsylinder > toppavsylinder.dat
echo $sylinderradius > sylinderradius.dat
echo $sylinderhoyde > sylinderhoyde.dat

dt=0.002
tid=100
ant=$(py "${tid}/${dt}")

antdumps=100
dumpfrekv=$(py "${ant}/${antdumps}")

echo "log log.oppsett
    #/atom 1 silicon
    #/atom 2 oxygen
    #/atom 3 hydrogen
    #/bond 1 2 2.6
    #/bond 2 3 1.4

    units metal
    boundary p p p
    atom_style molecular

    read_data ../amorf.data

    # pair_style  usc
    # pair_coeff  * * ../SiOH2O.vashishta Si O H
    # pair_modify coord 2 1 2.0 0.3
    # pair_modify coord 2 3 1.4 0.3
    pair_style vashishta
    pair_coeff	* * ../SiO2.vashishta Si O
    mass            1 28.08
    mass            2 15.9994
    #mass            3 1.00794

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
    group beholdes union bunn sylinder
    group slettes subtract all beholdes

    delete_atoms group slettes

    group kanbevegeseg subtract all hardbunn


    thermo 10
    thermo_style custom step time temp press pzz etotal cpuremain
    dump lagring all custom ${dumpfrekv} oppsett.in.bin id type x y z

    group silisium type 1
    group oksygen type 2

    variable antsilisium equal count(silisium)
    variable antoksygen equal count(oksygen)
    variable antsilisiumx2 equal '2*count(silisium)'

    region skalfylles union 3 tomromz bunn sylinder side out
    group forstokiometrisering region helesystemet


    if '\${antoksygen} > \${antsilisiumx2}' then &
        \"variable antsilisiumsommangler equal '(v_antoksygen - v_antsilisiumx2)/2'\" &
        \"create_atoms 1 random \${antsilisiumsommangler} 142857 skalfylles\" &
    elif '\${antoksygen} < \${antsilisiumx2}' &
        \"variable antoksygensommangler equal 'v_antsilisiumx2 - v_antoksygen'\" &
        \"create_atoms 2 random \${antoksygensommangler} 142857 skalfylles\"

    group nyeatomer subtract all forstokiometrisering
    group kanbevegeseg union nyeatomer

    fix termostat kanbevegeseg nvt temp $T $T 1.0
    velocity kanbevegeseg create $T 277385 mom yes loop geom

    group silisium type 1
    group oksygen type 2

    timestep ${dt}
    run ${ant}


    # group sio2 region helesystemet
    # group OHgrupper subtract all sio2

    # group kanbevegeseg union kanbevegeseg OHgrupper

    write_restart klartilpassivering.restart
    " > data/oppsett.in
