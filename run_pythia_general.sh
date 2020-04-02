#! /bin/bash
SOURCEDIR=$1
VARIATION=$2
SEED=$3
PDFSET=$4
PYTHIATUNE=$5

source $CONF/alice_setenv
PACKAGES=(pythia fastjet ROOT lhapdf-pdfsets)
DEFAULTS=root6
ALIENV=/usr/local/bin/alienv
for pack in ${PACKAGES[@]}; do
    eval `$ALIENV --no-refresh printenv $pack/latest-ali-master-$DEFAULTS` 
done
# Additional PDF sets locally installed
source $CONF/lhapdf_data_setenv

export CONFIG_SEED=$SEED
LOGFILE=$(printf "pythia8_%s_%s_%s.log" $PDFSET $PYTHIATUNE $VARIATION)
cmd=$(printf "root -l -b -q \'%s/RunPythia8.C(\"%s\", \"%s\", \"%s\")\' &> %s" $SOURCEDIR $VARIATION $PDFSET $PYTHIATUNE $LOGFILE)
eval $cmd
