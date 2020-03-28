#! /bin/bash

TAG=$1
source $HOME/alice_setenv
eval `/usr/bin/alienv --no-refresh printenv POWHEG/latest`
/usr/bin/alienv list
source $HOME/lhapdf_data_setenv

pwhg_main_dijet &> powheg_$TAG.log
