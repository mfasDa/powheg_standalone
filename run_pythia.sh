#! /bin/bash

SOURCEDIR=$1
WORKDIR=$2
SEED=$3

TMPDIR=/tmp/slurm_job_$SLURM_JOBID
if [ ! -d $TMPDIR ]; then mkdir -p $TMPDIR; fi
cd $TMPDIR

ALIBUILD_WORK_DIR=/clusterfs1/markus/alice/sw
PACKAGES=(pythia fastjet ROOT lhapdf-pdfsets)
for pack in ${PACKAGES[@]}; do
    eval `/usr/local/bin/alienv -w $ALIBUILD_WORK_DIR --no-refresh printenv $pack/latest-ali-master-root6` 
done
# Additional PDF sets locally installed
export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/clusterfs1/markus/alice/PDFSETS

export CONFIG_SEED=$SEED
cp $WORKDIR/powheg.zip $PWD
unzip powheg.zip

VARIATIONS=("main" "rlfl" "rlfm" "rlfh" "rmfl" "rmfm" "rmfh" "rhfl" "rhfm" "rhfh")
LOGFILES=()
ROOTFILES=()
for VAR in ${VARIATIONS[@]}; do
    echo "Doing variation $VAR"
    OUTPUTFILE=$(printf "POWHEGPYTHIA_%s.root" $VAR)
    LOGFILE=$(printf "pythia8_%s.log" $VAR)
    cmd=$(printf "root -l -b -q \'%s/RunPythia8.C(\"%s\", \"%s\")\' &> %s" $SOURCEDIR $OUTPUTFILE $VAR $LOGFILE)
    eval $cmd
    LOGFILES+=($LOGFILE)
    ROOTFILE+=($ROOTFILE)
done

# Copy back output
echo "Packing and copying back output"
rootcmd="zip pythia_roots.log "
for ROOTFILE in ${ROOTFILES[@]}; do
    rootcmd=$(printf "%s %s" "$rootcmd" $ROOTFILE)
done
echo "Command for root package: $rootcmd"
eval $rootcmd
cp -v pythia_roots.zip $WORKDIR
logcmd="zip pythia_logs.zip "
for LOGFILE in ${LOGFILES[@]}; do
    logcmd=$(printf "%s %s" "$logcmd" $LOGFILE)
done
echo "Command for log package: $logcmd"
eval $logcmd
cp -v pythia_logs.zip $WORKDIR

cd $WORKDIR
rm -rf $TMPDIR
echo "Done ..."