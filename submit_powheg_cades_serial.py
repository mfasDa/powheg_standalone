#! /usr/bin/env python3
import argparse
import logging
import math
import os
import sys
import subprocess

sourcedir = os.path.dirname(os.path.abspath(sys.argv[0]))

def createJobscript(outputdir, maxtime, slot, nevents, ebeam, pdfset, muf, mur, bornkt, bornsupp, withnegweight):
    slotdir = "%04d" %slot
    workdir = os.path.join(outputdir, slotdir)
    if not os.path.exists(workdir):
        os.makedirs(workdir, 0o755)
    jobscriptname = os.path.join(workdir, "jobscript.sh")
    logfile = os.path.join(workdir, "joboutput.log")
    with open(os.path.join(workdir, "jobscript.sh"), "w") as exewriter:
        exewriter.write("#! /bin/bash\n")
        exewriter.write("#SBATCH -A birthright\n")
        exewriter.write("#SBATCH -N 1\n")
        exewriter.write("#SBATCH -n 1\n")
        exewriter.write("#SBATCH -c 1\n")
        exewriter.write("#SBATCH -p gpu\n")
        exewriter.write("#SBATCH -J powheg\n")
        exewriter.write("#SBATCH -o {}\n".format(logfile))
        exewriter.write("#SBATCH -t {}\n".format(maxtime))
        exewriter.write("#SBATCH --mem=2G\n")
        exewriter.write("module load python/3.6.3\n")
        exewriter.write("module load PE-gnu\n")
        exewriter.write("module load singularity\n")
        exewriter.write("echo \"Using working directory {}\"\n".format(workdir))
        exewriter.write("cd {}\n".format(workdir))
        exewriter.write("echo \"Preparing working directories and configurations ...\"\n")
        configcreator = "{}/create_powheginput.py -n {} -s {} -e {} -p {} -f {} -r {} -k {} -b {}".format(sourcedir, nevents, slot, ebeam, pdfset, muf, mur, bornkt, bornsupp) 
        if withnegweight:
            configcreator += " -w"
        exewriter.write("python3 {}\n".format(configcreator))
        exewriter.write("echo \"Running simulation ... \"\n")
        exewriter.write("singularity exec -B /nfs/home:/nfs/home -B /lustre:/lustre /home/mfasel_alice/mfasel_cc7_alice.simg {}/run_powheg_general.sh &> powheg.log\n".format(sourcedir))
        exewriter.write("rm -v {}\n".format(jobscriptname))
        exewriter.write("echo Job done\n")
        exewriter.close() 
    return jobscriptname

if __name__ == "__main__":
    logging.basicConfig(format='[%(levelname)s]: %(message)s', level=logging.INFO)
    parser = argparse.ArgumentParser("submit_powheg.py", "Submitter for POWHEG dijet process")
    parser.add_argument("-j", "--jobs", metavar="NJOBS", type=int, default = 10, help="Number of jobs")
    parser.add_argument("-n", "--nevents", metavar="NEVENTS", type=int, default=50000, help="Number of events")
    parser.add_argument("-e", "--ebeam", metavar="EBEAM", type=float, default=6500., help="Beam energy")
    parser.add_argument("-o", "--outputdir", metavar="OUTPUTDIR", type=str, required=True, help="Output directory")
    parser.add_argument("-p", "--pdfset", metavar="PDFSET", type=int, default=10550, help="PDF set")
    parser.add_argument("-f", "--mufact", metavar="MUFACT", type=float, default=1., help="Factorization scale")
    parser.add_argument("-r", "--muren", metavar="MUREN", type=float, default=1., help="Renomalization scale")
    parser.add_argument("-k", "--bornkt", metavar="BORNKT", type=float, default=1., help="Born kt (GeV)")
    parser.add_argument("-s", "--bornsupp", metavar="BORNSUPP", type=float, default=70., help="Born suppression")
    parser.add_argument("-g", "--withnegweight", action='store_true', help="Include negative weights")
    parser.add_argument("-t", "--time", metavar="TIME", default="10:00:00", help="Max. time")
    args = parser.parse_args()
    outputdir = os.path.abspath(args.outputdir)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir, 0o755)
    for jobID in range(0, args.jobs):
        logging.info("Processing job %d", jobID)
        jobscript = createJobscript(outputdir, args.time, jobID, args.nevents, args.ebeam, args.pdfset, args.mufact, args.muren, args.bornkt, args.bornsupp, args.withnegweight)
        subprocess.call(["sbatch", jobscript])
