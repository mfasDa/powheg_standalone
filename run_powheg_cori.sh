#! /bin/bash

source /global/homes/m/mfasel/alice_setenv
eval `/usr/bin/alienv --no-refresh printenv POWHEG/latest`
/usr/bin/alienv list
export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/global/cfs/projectdirs/alice/mfasel/alice/PDFSETS

pwhg_main_dijet &> powheg_main.log