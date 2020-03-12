#! /usr/bin/env python3

import argparse
import os
import sys
import subprocess

sourcedir = os.path.dirname(os.path.abspath(sys.argv[0]))

def createJobscript(workdir, nevents, ebeam, pdf, bornkt, bornsupp):
    jobscriptname = os.path.join(workdir, "jobscript.sh")
    with open(jobscriptname, "w") as writer:
        writer.write("#! /bin/bash\n")
        writer.write("#SBATCH -n 1\n")
        writer.write("#SBATCH -o %s/joboutput.log\n" %workdir)
        writer.write("echo Running on $HOSTNAME\n")
        writer.write("%s/run_powheg.sh %s %s $SLURM_JOBID %d %.1f %d %.1f %.1f\n" %(sourcedir, sourcedir, workdir, nevents, ebeam, pdf, bornkt, bornsupp))
        writer.write("rm -rf %s\n" %jobscriptname)
        writer.close()
    return jobscriptname

if __name__ == "__main__":
    parser = argparse.ArgumentParser("submit_powheg.py", "Submitter for POWHEG dijet process")
    parser.add_argument("-j", "--jobs", metavar="NJOBS", type=int, default = 10, help="Number of jobs")
    parser.add_argument("-n", "--nevents", metavar="NEVENTS", type=int, default=50000, help="Number of events")
    parser.add_argument("-e", "--ebeam", metavar="EBEAM", type=float, default=6500., help="Beam energy")
    parser.add_argument("-o", "--outputdir", metavar="OUTPUTDIR", type=str, required=True, help="Output directory")
    parser.add_argument("-p", "--pdfset", metavar="PDFSET", type=int, default=10550, help="PDF set")
    parser.add_argument("-k", "--bornkt", metavar="BORNKT", type=float, default=1., help="Born kt (GeV)")
    parser.add_argument("-s", "--bornsupp", metavar="BORNSUPP", type=float, default=70., help="Born suppression")
    args = parser.parse_args()
    outputdir = os.path.abspath(args.outputdir)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir, 0o755)
    for jid in range(0, args.jobs):
        workdir = os.path.join(outputdir, "%04d" %jid)
        if not os.path.exists(workdir):
            os.makedirs(workdir, 0o755)
        jobscript = createJobscript(workdir, args.nevents, args.ebeam, args.pdfset, args.bornkt, args.bornsupp)
        subprocess.call(["sbatch", jobscript])
