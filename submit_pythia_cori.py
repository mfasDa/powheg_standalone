#! /usr/bin/env python3
import argparse
import logging
import math
import os
import sys
import subprocess

sourcedir = os.path.dirname(os.path.abspath(sys.argv[0]))

def createMasterJobscript(workdir, platform, queue, maxtime, nslots, minslot, masterID, variation, pdfset):
    if not os.path.exists(workdir):
        os.makedirs(workdir, 0o755)
    jobscriptname = os.path.join(workdir, "jobscript_pythia_{}_{}_master{}.sh".format(pdfset, variation, masterID))
    logfile = os.path.join(workdir, "joboutput_pythia_{}_{}_master{}.log".format(pdfset, variation, masterID))
    jobname = "pythia8_{}_{}".format(pdfset, variation)
    with open(jobscriptname, "w") as exewriter:
        exewriter.write("#! /bin/bash\n")
        exewriter.write("#SBATCH -N 1\n")
        exewriter.write("#SBATCH --ntasks-per-node={}\n".format(nslots))
        exewriter.write("#SBATCH --qos={}\n".format(queue))
        exewriter.write("#SBATCH -J powheg\n")
        exewriter.write("#SBATCH -o {}\n".format(logfile))
        exewriter.write("#SBATCH -t {}\n".format(maxtime))
        exewriter.write("#SBATCH -C {}\n".format(platform))
        exewriter.write("#SBATCH --image=docker:mfasel/cc7-alice:latest\n")
        exewriter.write("module load cray-python/3.7.3.2\n")
        exewriter.write("WORKDIR={}\n".format(workdir))
        exewriter.write("echo \"Running pythia on existing POWHEG {}\"\n".format(workdir))
        exewriter.write("cd $WORKDIR\n")
        exewriter.write("echo \"Running showering on {} cores, starting from slot {} ... \"\n".format(nslots, minslot))
        exewriter.write("SECONDS=0\n")
        exewriter.write("srun -n {} python3 {}/mpiwrapper_pythia.py {} {} {}\n".format(nslots, sourcedir, minslot, variation, pdfset))
        exewriter.write("duration=$SECONDS\n")
        exewriter.write("cd {}\n".format(workdir))
        exewriter.write("echo Job done after $duration seconds\n")
        exewriter.close() 
    return jobscriptname

if __name__ == "__main__":
    logging.basicConfig(format='[%(levelname)s]: %(message)s', level=logging.INFO)
    parser = argparse.ArgumentParser("submit_powheg.py", "Submitter for POWHEG dijet process")
    parser.add_argument("-i", "--inputdir", metavar="INPUTDIR", type=str, required=True, help="Directory with POWHEG events")
    parser.add_argument("-p", "--pdfset", metavar="PDFSET", type=str, default="CT14nlo", help="PDF set")
    parser.add_argument("-v", "--variation", metavar="VARIATION", type=str, default="main", help="Scale variation")
    parser.add_argument("-q", "--qos", metavar="QOS", default="regular", help="Cori queue")
    parser.add_argument("-t", "--time", metavar="TIME", default="10:00:00", help="Max. time")
    parser.add_argument("-c", "--constraint", metavar="CONSTRAINT", default="knl", help="Platform (haswell or knl)")
    args = parser.parse_args()
    workdir = os.path.abspath(args.inputdir)
    slotspermaster = 68
    if args.constraint == "knl":
        slotspermaster = 68
    elif args.constraint == "haswell":
        slotspermaster = 32
    else:
        print("Unknown architecture, select either haswell or knl")
        sys.exit(1)
    maxtime = args.time
    if args.qos == "debug":
        maxtime = "00:30:00"
    dirs = os.listdir(workdir)
    njobs = 0
    for d in dirs:
        if not d.isdigit():
            continue
        njobs += 1
    logging.info("Running %d jobs", njobs)
    nmasters = math.floor(int(njobs) / slotspermaster)
    if njobs % slotspermaster != 0:
        nmasters += 1
    minslot = 0
    for masterID in range(0, nmasters):
        logging.info("Processing master %d", masterID)
        nslotsworker = njobs - minslot
        if nslotsworker > slotspermaster:
            nslotsworker = slotspermaster
        jobscript = createMasterJobscript(workdir, args.constraint, args.qos, maxtime, nslotsworker, minslot, masterID, args.variation, args.pdfset)
        subprocess.call(["sbatch", jobscript])
        minslot += nslotsworker
