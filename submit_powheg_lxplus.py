#! /usr/bin/env python3
import argparse
import os
import sys
import subprocess

sourcedir = os.path.dirname(os.path.abspath(sys.argv[0]))

def createJobscript(workdir, outputdir, njobs, nevents, ebeam, pdf, bornkt, bornsupp, withnegweight):
    jobscriptname = os.path.join(workdir, "job.submit")
    logdir = os.path.join(workdir, "logs")
    if not os.path.exists(logdir):
        os.makedirs(logdir, 0o755)
    for j in range(0, njobs):
        joboutoutdir=os.path.join(outputdir, "%04d" %j)
        if not os.path.exists(joboutoutdir):
            os.makedirs(joboutoutdir, 0o755)
    with open(os.path.join(workdir, "jobscript.sh"), "w") as exewriter:
        exewriter.write("#! /bin/bash\n")
        exewriter.write("SLOT=$1\n")
        exewriter.write("echo Running on $HOSTNAME\n")
        exewriter.write("echo $WORKDIR\n")
        exewriter.write("%s/run_powheg_lxplus.sh %s %d $SLOT %.1f %d %.1f %.1f %d\n" %(sourcedir, sourcedir, nevents, ebeam, pdf, bornkt, bornsupp, 1 if withnegweight else 0))
        exewriter.write("OUTPUTDIR=$(printf \"{}/%04d\" $SLOT)\n".format(outputdir))
        exewriter.write("cp powheg.zip $OUTPUTDIR\n")
        exewriter.write("echo Job done\n")
        exewriter.close() 
    with open(jobscriptname, "w") as writer:
        writer.write("executable = %s/jobscript.sh\n" %workdir)
        writer.write("arguments = $(ProcId)\n")
        writer.write("log = %s/job_$(ClusterId)_$(ProcId).log\n" %logdir)
        writer.write("output = %s/job_$(ClusterId)_$(ProcId).out\n" %logdir)
        writer.write("error = %s/job_$(ClusterId)_$(ProcId).err\n" %logdir)
        writer.write("requirements = (OpSysAndVer =?= \"CentOS7\")\n")
        writer.write("+JobFlavour = \"tomorrow\"\n")
        writer.write("queue %d\n" %njobs)
        writer.close()
    return jobscriptname

if __name__ == "__main__":
    parser = argparse.ArgumentParser("submit_powheg.py", "Submitter for POWHEG dijet process")
    parser.add_argument("-j", "--jobs", metavar="NJOBS", type=int, default = 10, help="Number of jobs")
    parser.add_argument("-n", "--nevents", metavar="NEVENTS", type=int, default=50000, help="Number of events")
    parser.add_argument("-e", "--ebeam", metavar="EBEAM", type=float, default=6500., help="Beam energy")
    parser.add_argument("-o", "--outputdir", metavar="OUTPUTDIR", type=str, required=True, help="Output directory")
    parser.add_argument("-w", "--workdir", metavar="WORKDIR", type=str, required=True, help="Working directory (on AFS)")
    parser.add_argument("-p", "--pdfset", metavar="PDFSET", type=int, default=10550, help="PDF set")
    parser.add_argument("-k", "--bornkt", metavar="BORNKT", type=float, default=1., help="Born kt (GeV)")
    parser.add_argument("-s", "--bornsupp", metavar="BORNSUPP", type=float, default=70., help="Born suppression")
    parser.add_argument("-g", "--withnegweight", action='store_true', help="Include negative weights")
    args = parser.parse_args()
    outputdir = os.path.abspath(args.outputdir)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir, 0o755)
    jobscript = createJobscript(args.workdir, outputdir, args.jobs, args.nevents, args.ebeam, args.pdfset, args.bornkt, args.bornsupp, args.withnegweight)
    subprocess.call(["condor_submit", jobscript])
