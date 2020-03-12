#! /bin/bash

SOURCEDIR=$1
WORKDIR=$2
SEED=$3
NEVENTS=$4
ENERGY=$5
PDF=$6
BORNKT=$7
BORNSUPP=$8

if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
cd $WORKDIR

eval `/usr/local/bin/alienv -w /clusterfs1/markus/alice/sw --no-refresh printenv POWHEG/latest-ali-master-root6`

echo "Running central powheg prediction"
$SOURCEDIR/create_powheginput.py -n $NEVENTS -s $SEED -e $ENERGY -p $PDF -k $BORNKT -b $BORNSUPP
pwhg_main_dijet &> powheg_main.log
mv powheg.input powheg.input.default

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
    $SOURCEDIR/create_powheginput.py -n $NEVENTS -s $SEED -e $ENERGY -p $PDF -k $BORNKT -b $BORNSUPP -f $mymuf -r $mymur -c $myscaleid
    pwhg_main_dijet &> powheg_scale_$myscaleid.log
    mv powheg.input powheg.input.scale.$myscaleid
    mv pwgevents-rwgt.lhe pwgevents.lhe
done

# Generate configs for pdf variation
PDFSETS=(13100 13101 13102 13103 13104 13105 13106 13107 13108 13109 13110 13111 13112 13113 13114 13115 13116 13117 13118 13119 13120 13121 13122 13123 13124 13125 13126 13127 13128 13129 13130 13131 13132 13133 13134 13135 13136 13137 13138 13139 13140 13141 13142 13143 13144 13145 13146 13147 13148 13149 13150 13151 13152 13153 13154 13155 13156 )
for pdf in ${PDFSETS[@]}; do
    echo "Doing variation pdf set: $pdf"
    $SOURCEDIR/create_powheginput.py -n $NEVENTS -s $SEED -e $ENERGY -p $pdf -k $BORNKT -b $BORNSUPP -c $pdf
    pwhg_main_dijet &> powheg_pdf_$pdf.log
    mv powheg.input powheg.input.pdf$pdf
    mv pwgevents-rwgt.lhe pwgevents.lhe
done
