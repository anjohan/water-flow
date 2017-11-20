#!/bin/bash
set -e -x

source regn.sh

if [ ! -d data ]
then
    mkdir data
fi

cp job_amorf.pbs data/job_amorf.pbs

source system.sh
startT=$T
toppT=4000
dt=0.002
topp=30
antopp=$(py "${topp}/${dt}")
ttermalisering=20
anttermalisering=$(py "${ttermalisering}/${dt}")
tned=120
antned=$(py "${tned}/${dt}")
ttermaliseringslutt=50
anttermaliseringslutt=$(py "${ttermaliseringslutt}/${dt}")
ttot=$(py "${topp}+${ttermalisering}+${tned}+${ttermaliseringslutt}")
antsteg=$(py "${ttot}/${dt}")
antdumps=100
dumpfrekv=$(py "${antsteg}/${antdumps}")

echo "log log.amorf
    #/atom 1 silicon
    #/atom 2 oxygen
    #/atom 3 hydrogen
    #/bond 1 2 2.6
    #/bond 2 3 1.4

    #package gpu 2
    #suffix gpu

    units metal
    boundary p p p
    atom_style molecular

    read_data ../betacristobalite_utenH.data
    replicate $antx $anty $systemhoyde

    pair_style vashishta
    pair_coeff	* * ../SiO2.vashishta Si O
    mass            1 28.08
    mass            2 15.9994

    neigh_modify every 1 delay 0 check yes
    timestep ${dt}

    thermo 100
    thermo_style custom step time temp press pzz etotal spcpu cpuremain
    dump lagring all custom ${dumpfrekv} amorf.in.bin id type x y z

    velocity all create ${startT} 277385 mom yes loop geom

    fix termostat all nvt temp $startT $toppT 1.0
    run ${antopp}

    fix termostat all nvt temp $toppT $toppT 1.0
    run ${anttermalisering}

    fix termostat all nvt temp $toppT $T 1.0
    run ${antned}

    fix termostat all nvt temp $T $T 1.0
    run ${anttermaliseringslutt}

    write_data amorf.data
    write_restart amorf.restart
    " > data/amorf.in
