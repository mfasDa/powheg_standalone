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
    logging.info("Starting worker %d ...", CURRENTSLOT)
    slotdir = "%04d" %CURRENTSLOT
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
    subprocess.call("shifter {} &> run_powheg.log".format(os.path.join(SOURCEDIR, "run_powheg_cori.sh")), shell = True)
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