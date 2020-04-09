#! /usr/bin/env python

from mpi4py import MPI
import logging
import os
import subprocess
import sys

if __name__ == "__main__":
    logging.basicConfig(format='[%(levelname)s]: %(message)s', level=logging.INFO)
    SOURCEDIR = os.path.dirname(os.path.abspath(sys.argv[0]))
    CURRENTSLOT = MPI.COMM_WORLD.Get_rank()
    MINSLOT = int(sys.argv[1])
    VARIATION = sys.argv[2]
    PDFSET = sys.argv[3]
    TUNE = sys.argv[4]
    FILE =  sys.argv[5]
    logging.info("Starting worker %d ...", CURRENTSLOT)
    GLOBALSLOT = MINSLOT+CURRENTSLOT
    slotdir = "%04d" %GLOBALSLOT
    os.chdir(slotdir)
    logging.info("Running job for slot %d in workdir %s", CURRENTSLOT, os.getcwd())
    content = os.listdir(os.getcwd())
    if not len(content):
        logging.info("Working directory %s empty ...", os.getcwd())
    else:
        contentstring = ""
        for f in content:
            if len(contentstring):
                contentstring += ", "
            contentstring += f
        logging.info("Content of working directory %s: %s", os.getcwd(), contentstring)
    content = os.listdir(os.getcwd())
    subprocess.call("shifter {} {} {} {} {} {} {} &> run_pythia_{}_{}_{}.log".format(os.path.join(SOURCEDIR, "run_pythia_general.sh"), SOURCEDIR, VARIATION, GLOBALSLOT, PDFSET, TUNE, FILE, PDFSET, TUNE, VARIATION), shell = True)
    content = os.listdir(os.getcwd())
    if not len(content):
        logging.info("Working directory %s empty after job execution ...", os.getcwd())
    else:
        contentstring = ""
        for f in content:
            if len(contentstring):
                contentstring += ", "
            contentstring += f
        logging.info("Content of working directory after job execution %s: %s", os.getcwd(), contentstring)
    logging.info("Worker %d done ...", CURRENTSLOT)
