#! /usr/bin/env python3

"""
Example powheg.input

!randomseed 352345 ! uncomment to set the random seed to a value of your choice.
                   ! It generates the call RM48IN(352345,0,0) (see the RM48 manual).
                   ! THIS MAY ONLY AFFECTS THE GENERATION OF POWHEG EVENTS!
                   ! If POWHEG is interfaced to a shower MC, refer to the shower MC
                   ! documentation to set its seed.

!iseed 1

numevts 20000            ! number of events to be generated

ih1  1               ! hadron 1 (1 for protons, -1 for antiprotons)
ih2  1               ! hadron 2 (1 for protons, -1 for antiprotons)
ebeam1 6500d0        ! energy of beam 1
ebeam2 6500d0        ! energy of beam 2

bornktmin 1d0        ! (default 0d0) Generation cut: minimum kt in underlying Born
bornsuppfact 70d0    ! (default 0d0) Mass parameter for Born suppression factor.
                     ! If < 0 suppfact = 1.


! To be set only if using internal (mlm) pdfs
! 131 cteq6m
! ndns1 131            ! pdf set for hadron 1 (mlm numbering)
! ndns2 131            ! pdf set for hadron 2 (mlm numbering)

! To be set only if using LHA pdfs
! 10050  cteq6m
! 13100  CT14nlo
! 10550  cteq66-0
! 247000 NNPDF23_lo_as_0130_qed
 lhans1  10050      ! pdf set for hadron 1 (LHA numbering)
 lhans2  10050      ! pdf set for hadron 2 (LHA numbering)
! lhans1  13100      ! pdf set for hadron 1 (LHA numbering)
! lhans2  13100      ! pdf set for hadron 2 (LHA numbering)	
! lhans1  10550      ! pdf set for hadron 1 (LHA numbering)
! lhans2  10550      ! pdf set for hadron 2 (LHA numbering)	
! lhans1  247000       ! pdf set for hadron 1 (LHA numbering)	
! lhans2  247000       ! pdf set for hadron 2 (LHA numbering)	

! To be set only if using different pdf sets for the two incoming hadrons
#QCDLambda5  0.25    ! for not equal pdf sets 

#renscfact  0.5d0      ! (default 1d0) ren scale factor: muren  = muref * renscfact 
#facscfact  0.5d0      ! (default 1d0) fac scale factor: mufact = muref * facscfact 
renscfact  1.0      ! (default 1d0) ren scale factor: muren  = muref * renscfact 
facscfact  1.0      ! (default 1d0) fac scale factor: mufact = muref * facscfact 

! Parameters to allow or not the use of stored data
use-old-grid    1    ! If 1 use old grid if file pwggrids.dat is present (<> 1 regenerate)
use-old-ubound  1    ! If 1 use norm of upper bounding function stored
                     ! in pwgubound.dat, if present; <> 1 regenerate

! A typical call uses 1/1400 seconds (1400 calls per second)
ncall1 20000         ! No. calls for the construction of the importance sampling grid
itmx1 5              ! No. iterations for grid: total 100000 calls ~ 70 seconds
ncall2 20000         ! No. calls for the computation of the upper bounding
                     ! envelope for the generation of radiation
itmx2 5              ! No. iterations for the above

! Notice: the total number of calls is ncall2*itmx2*foldcsi*foldy*foldphi
! these folding numbers yield a negative fraction of 0.5% with bornktmin=10 GeV.
! With these settings: ncall2*itmx2*foldcsi*foldy*foldphi=5M, 60 minutes
foldcsi 5            ! No. folds on csi integration
foldy   5            ! No. folds on  y  integration
foldphi 2            ! No. folds on phi integration

nubound 500000       ! No. calls to set up the upper bounding norms for radiation.
                     ! This is performed using only the Born cross section (fast)

! OPTIONAL PARAMETERS

withnegweights 1    ! (default 0). If 1 use negative weights.
#bornonly 1          ! (default 0). If 1 compute underlying Born using LO
                     ! cross section only.

#ptsqmin    0.8      ! (default 0.8 GeV) minimum pt for generation of radiation
#charmthr   1.5      ! (default 1.5 GeV) charm treshold for gluon splitting 
#bottomthr  5.0      ! (default 5.0 GeV) bottom treshold for gluon splitting
#testplots  1        ! (default 0, do not) do NLO and PWHG distributions
#charmthrpdf  1.5    ! (default 1.5 GeV) pdf charm treshold
#bottomthrpdf 5.0    ! (default 5.0 GeV) pdf bottom treshold

#xupbound 2d0        ! increase upper bound for radiation generation

#iseed    5421       ! Start the random number generator with seed iseed
#rand1     0         ! skipping  rand2*100000000+rand1 numbers (see RM48
#rand2     0         ! short writeup in CERNLIB).              
#manyseeds 1         ! Used to perform multiple runs with different random
                     ! seeds in the same directory.
                     ! If set to 1, the program asks for an integer j;
                     ! The file pwgseeds.dat at line j is read, and the
                     ! integer at line j is used to initialize the random
                     ! sequence for the generation of the event.
                     ! The event file is called pwgevents-'j'.lhe

doublefsr 1          ! Default 0; if 1 use new mechanism to generate regions
                     ! such that the emitted harder than the
                     ! emitter in FSR is suppressed. If doublefsr=0 this is
                     ! only the case for emitted gluons (old behaviour). If
                     ! 1 it is also applied to emitted quarks.
                     ! If set, it strongly reduces spikes on showered output.



#par_diexp 4         ! default is 2. With 4 there is a stronger separation
#par_dijexp 4        ! of regions, it may help to reduce spikes when generating
#par_2gsupp 4        ! weighted events.

#storeinfo_rwgt 1
compute_rwgt 1

lhrwgt_id 'lfr'
lhrwgt_descr 'estimate uncertainties'
lhrwgt_group_name 'errors'
lhrwgt_group_combine 'foo'
"""

import argparse
import os
import sys
import subprocess

def make_pwhgin(numevents, seed, ebeam, pdfset, mufact, muren, bornkt, bornsupp, withnegweights, weightid):
    with open("powheg.input", "w") as configwriter:
        configwriter.write("randomseed %d\n" %seed)
        configwriter.write("iseed %d\n" %seed)
        configwriter.write("numevts %d\n" %numevents)
        configwriter.write("ih1  1\n")
        configwriter.write("ih2  1\n")
        configwriter.write("ebeam1 %.1f\n" %ebeam)
        configwriter.write("ebeam2 %.1f\n" %ebeam)
        configwriter.write("lhans1  %d\n" %pdfset)
        configwriter.write("lhans2  %d\n" %pdfset)
        configwriter.write("bornktmin %.1f\n" %bornkt)
        configwriter.write("bornsuppfact %.1f\n" %bornsupp)
        configwriter.write("renscfact  %.1f\n" %muren)
        configwriter.write("facscfact  %.1f\n" %mufact)
        configwriter.write("use-old-grid    1\n")
        configwriter.write("use-old-ubound  1\n")
        configwriter.write("ncall1 20000\n")
        configwriter.write("itmx1 5\n");
        configwriter.write("ncall2 20000\n")
        configwriter.write("itmx2 5\n");
        configwriter.write("foldcsi 5\n")
        configwriter.write("foldy   5\n")
        configwriter.write("foldphi 2\n")
        configwriter.write("doublefsr 1\n")
        configwriter.write("nubound 500000\n")
        if withnegweights:
            configwriter.write("withnegweights 1\n")
        else:
            configwriter.write("withnegweights 0\n")
        if len(weightid):
            configwriter.write("compute_rwgt 1\n")
            configwriter.write("lhrwgt_id \'%s\'\n" %weightid)
        else:
            configwriter.write("lhrwgt_id \'main\'\n")
        configwriter.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser("create_powheginput.py", "tool generating powheg.input for various configurations")
    parser.add_argument("-n", "--numevents", metavar="NUMEVENTS", type=int, required=True, help="Number of events")
    parser.add_argument("-s", "--seed", metavar="SEED", type=int, required=True, help="Random seed")
    parser.add_argument("-e", "--ebeam", metavar="BEAMEMERGY", type=float, required=True, help="Beam energy (GeV)")
    parser.add_argument("-p", "--pdfset", metavar="PDFSET", type=int, default=10050, help="PDF set")
    parser.add_argument("-f", "--mufact", metavar="MUFACT", type=float, default=1., help="Factorization scale")
    parser.add_argument("-r", "--muren", metavar="MUREN", type=float, default=1., help="Renomalization scale")
    parser.add_argument("-k", "--bornkt", metavar="BORNKT", type=float, default=1., help="Born-kt (GeV)")
    parser.add_argument("-b", "--bornsupp", metavar="BORNSUPP", type=float, default=70., help="Born suppression (GeV)")
    parser.add_argument("-c", "--recalcweight", metavar="RECALCWEIGHT", type=str, default="", help="Recalculate weights")
    parser.add_argument("-w", "--withnegweights", action='store_true', help="Use events with negative weights")
    options = parser.parse_args()
    make_pwhgin(options.numevents, options.seed, options.ebeam, options.pdfset, options.mufact, options.muren, options.bornkt, options.bornsupp, options.withnegweights, options.recalcweight)
