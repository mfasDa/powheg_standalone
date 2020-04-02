#! /usr/bin/env python3
import argparse
import logging
import math
import os
import sys
import subprocess

sourcedir = os.path.dirname(os.path.abspath(sys.argv[0]))

def createJobscript(workdir, seed, variation, pdfset, tune):
    jobscriptname = os.path.join(workdir, "jobscript_{}_{}_{}.sh".format(pdfset, tune, variation))
    logfile = os.path.join(workdir, "joboutput_pythia_{}_{}_{}.log".format(pdfset, tune, variation))
    with open(os.path.join(workdir, jobscriptname), "w") as exewriter:
        exewriter.write("#! /bin/bash\n")
        exewriter.write("#SBATCH -N 1\n")
        exewriter.write("#SBATCH -n 1\n")
        exewriter.write("#SBATCH -c 1\n")
        exewriter.write("#SBATCH -J pythia\n")
        exewriter.write("#SBATCH --partition=long\n")
        exewriter.write("#SBATCH -o {}\n".format(logfile))
        exewriter.write("export CONF=/clusterfs1/markus/alice\n")
        exewriter.write("echo \"Using working directory {}\"\n".format(workdir))
        exewriter.write("cd {}\n".format(workdir))
        exewriter.write("echo \"Preparing working directories and configurations ...\"\n")
        exewriter.write("echo \"Running simulation ... \"\n")
        exewriter.write("{}/run_pythia_general.sh {} {} {} {} {}  &> run_pythia_{}_{}_{}.log\n".format(sourcedir, sourcedir, variation, seed, pdfset, tune, pdfset, tune, variation))
        exewriter.write("rm -v {}\n".format(jobscriptname))
        exewriter.write("echo Job done\n")
        exewriter.close() 
    return jobscriptname

if __name__ == "__main__":
    logging.basicConfig(format='[%(levelname)s]: %(message)s', level=logging.INFO)
    parser = argparse.ArgumentParser("submit_powheg.py", "Submitter for POWHEG dijet process")
    parser.add_argument("-i", "--inputdir", metavar="INPUTDIR", type=str, required=True, help="Output directory")
    parser.add_argument("-v", "--variation", metavar="VARIATION", type=str, default="main", help="Scale variation")
    parser.add_argument("-p", "--pdfset", metavar="PDFSET", type=str, default="CT14nlo", help="PDF set")
    parser.add_argument("-t", "--tune", metavar="TUNE", type=str, default="Monash2013", help="PYTHIA tune")
    args = parser.parse_args()
    seed = 0
    for slot in sorted(os.listdir(args.inputdir)):
        slotdir = os.path.join(args.inputdir, slot) 
        if not os.path.isdir(slotdir):
            continue
        if not os.path.exists(os.path.join(slotdir, "pwgevents.lhe")):
            logging.error("Directory %s does not contain pwgevents.lhe", slotdir)
            continue
        logging.info("Processing job %s", slot)
        jobscript = createJobscript(slotdir, seed, args.variation, args.pdfset, args.tune)
        subprocess.call(["sbatch", jobscript])
        seed += 1
