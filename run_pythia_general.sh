#! /bin/bash
SOURCEDIR=$1
VARIATION=$2
SEED=$3

source $HOME/alice_setenv
PACKAGES=(pythia fastjet ROOT lhapdf-pdfsets)
for pack in ${PACKAGES[@]}; do
    eval `/usr/bin/alienv --no-refresh printenv $pack/latest-ali-master-next-root6` 
done
# Additional PDF sets locally installed
source $HOME/lhapdf_data_setenv

export CONFIG_SEED=$SEED
OUTPUTFILE=$(printf "POWHEGPYTHIA_%s.root" $VARIATION)
LOGFILE=$(printf "pythia8_%s.log" $VARIATION)
cmd=$(printf "root -l -b -q \'%s/RunPythia8.C(\"%s\", \"%s\")\' &> %s" $SOURCEDIR $OUTPUTFILE $VARIATION $LOGFILE)
eval $cmd
