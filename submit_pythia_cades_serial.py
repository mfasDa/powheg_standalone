#! /usr/bin/env python3
import argparse
import logging
import math
import os
import sys
import subprocess

sourcedir = os.path.dirname(os.path.abspath(sys.argv[0]))

def createJobscript(workdir, maxtime, seed, variation):
    jobscriptname = os.path.join(workdir, "jobscript.sh")
    logfile = os.path.join(workdir, "joboutput_pythia_{}.log".format(variation))
    with open(os.path.join(workdir, "jobscript.sh"), "w") as exewriter:
        exewriter.write("#! /bin/bash\n")
        exewriter.write("#SBATCH -A birthright\n")
        exewriter.write("#SBATCH -N 1\n")
        exewriter.write("#SBATCH -n 1\n")
        exewriter.write("#SBATCH -c 1\n")
        exewriter.write("#SBATCH -p gpu\n")
        exewriter.write("#SBATCH -J pythia\n")
        exewriter.write("#SBATCH -o {}\n".format(logfile))
        exewriter.write("#SBATCH -t {}\n".format(maxtime))
        exewriter.write("#SBATCH --mem=2G\n")
        exewriter.write("module load PE-gnu\n")
        exewriter.write("module load singularity\n")
        exewriter.write("echo \"Using working directory {}\"\n".format(workdir))
        exewriter.write("cd {}\n".format(workdir))
        exewriter.write("echo \"Preparing working directories and configurations ...\"\n")
        exewriter.write("echo \"Running simulation ... \"\n")
        exewriter.write("singularity exec -B /nfs/home:/nfs/home -B /lustre:/lustre /home/mfasel_alice/mfasel_cc7_alice.simg {}/run_pythia_general.sh {} {} {}  &> run_pythia_{}.log\n".format(sourcedir, sourcedir, variation, seed, variation))
        exewriter.write("rm -v {}\n".format(jobscriptname))
        exewriter.write("echo Job done\n")
        exewriter.close() 
    return jobscriptname

if __name__ == "__main__":
    logging.basicConfig(format='[%(levelname)s]: %(message)s', level=logging.INFO)
    parser = argparse.ArgumentParser("submit_powheg.py", "Submitter for POWHEG dijet process")
    parser.add_argument("-i", "--inputdir", metavar="INPUTDIR", type=str, required=True, help="Output directory")
    parser.add_argument("-v", "--variation", metavar="VARIATION", type=str, default="main", help="Scale variation")
    parser.add_argument("-t", "--time", metavar="TIME", default="10:00:00", help="Max. time")
    args = parser.parse_args()
    seed = 0
    for slot in os.listdir(args.inputdir):
        logging.info("Processing job %s", slot)
        slotdir = os.path.join(args.inputdir, slot) 
        jobscript = createJobscript(slotdir, args.time, seed, args.variation)
        subprocess.call(["sbatch", jobscript])
        seed += 1
