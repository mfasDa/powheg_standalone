#! /usr/bin/env python3

import os
import sys
import subprocess

sourcedir = os.path.dirname(os.path.abspath(sys.argv[0]))

def create_jobscript(workdir, jobname):
    jobscriptname = os.path.join(workdir, "jobscript.sh")
    logfilename = os.path.join(workdir, "pythia_master.log")
    with open(jobscriptname, "w") as scriptwriter:
        scriptwriter.write("#! /bin/bash\n")
        scriptwriter.write("#SBATCH -N 1\n")
        scriptwriter.write("#SBATCH -o %s\n" %logfilename)
        scriptwriter.write("#SBATCH --job-name=%s\n" %jobname)
        scriptwriter.write("#SBATCH --partition=long\n")
        scriptwriter.write("%s/run_pythia.sh %s %s $SLURM_JOBID\n" %(sourcedir, sourcedir, workdir))
        scriptwriter.write("rm -rf %s\n" %jobscriptname)
        scriptwriter.write("echo Done\n")
        scriptwriter.close()
    return jobscriptname

if __name__ == "__main__":
    workdir = sys.argv[1]
    jobname = "powhegpythia"
    if len(sys.argv) > 2:
        jobname = sys.argv[2]

    for jdir in os.listdir(workdir):
        jobdir = os.path.join(workdir, jdir)
        if os.path.exists(os.path.join(jobdir, "powheg.zip")):
            subprocess.call(["sbatch", create_jobscript(jobdir, jobname)])