#! /bin/bash

SOURCEDIR=$1
NEVENTS=$2
SEED=$3
ENERGY=$4
PDF=$5
BORNKT=$6
BORNSUPP=$7
WITHNEGWEIGHT=$8

export PATH=$PATH:/afs/cern.ch/work/m/mfasel/alice/alibuild/bin/
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mfasel/alice/alibuild/lib/python3.6/site-packages/
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mfasel/alice/alibuild/lib64/python3.6/site-packages/
eval `/afs/cern.ch/work/m/mfasel/alice/alibuild/bin/alienv -w /afs/cern.ch/work/m/mfasel/alice/sw --no-refresh printenv POWHEG/latest-ali-master-user-next-root6`

echo "Running central powheg prediction"
cmd=$(printf "python3 %s/create_powheginput.py -n %d -s %d -e %f -p %d -k %f -b %f" $SOURCEDIR $NEVENTS $SEED $ENERGY $PDF $BORNKT $BORNSUPP)
if [ $WITHNEGWEIGHT -gt 0 ]; then 
    cmd=$(printf "%s -w" "$cmd")
fi
eval $cmd

pwhg_main_dijet &> powheg_main.log
mv powheg.input powheg.input.main

# Generate configs for scale variation
echo "Running scale variations"
MUR=(0.5 0.5 0.5 1. 1. 1. 2. 2. 2.)
MUF=(0.5 1. 2. 0.5 1. 2. 0.5 1. 2.)
SCALEID=("rlfl" "rlfm" "rlfh" "rmfl" "rmfm" "rmfh" "rhfl" "rhfm" "rhfh")
for conf in `seq 0 8`; do
    mymuf=${MUF[$conf]}
    mymur=${MUR[$conf]}
    myscaleid=${SCALEID[$conf]}
    echo "Doing variation $myscaleid: MUF $mymuf, MUR $mymur"
    varcmd=$(printf "python3 %s/create_powheginput.py -n %d -s %d -e %f -p %d -k %f -b %f -f %f -r %f -c %s" $SOURCEDIR $NEVENTS $SEED $ENERGY $PDF $BORNKT $BORNSUPP $mymuf $mymur $myscaleid)
    if [ $WITHNEGWEIGHT -gt -0 ]; then
        varcmd=$(printf "%s -w" "$varcmd")
    fi
    eval $varcmd
    
    pwhg_main_dijet &> powheg_scale_$myscaleid.log
    mv powheg.input powheg.input.scale.$myscaleid
    mv pwgevents-rwgt.lhe pwgevents.lhe
done

# Pack output
OUTPUTFILES=(FlavRegList pwg-btlgrid.top pwgborngrid.top pwgcounters.dat pwgevents.lhe pwghistnorms.top pwgxgrid.dat realequiv virtequiv bornequiv pwg-stat.dat pwgboundviolations.dat pwggrid.dat pwgubound.dat pwhg_checklimits realequivregions)
zipcmd="zip powheg.zip"
for o in ${OUTPUTFILES[@]}; do
    zipcmd=$(printf "%s %s" "$zipcmd" $o)
done

INPUTFILES=($(ls -1 | grep input))
for i in ${INPUTFILES[@]}; do
    zipcmd=$(printf "%s %s" "$zipcmd" $i)
done

LOGFILES=($(ls -1 | grep log))
for l in ${LOGFILES[@]}; do
    zipcmd=$(printf "%s %s" "$zipcmd" $l)
done
eval $zipcmd