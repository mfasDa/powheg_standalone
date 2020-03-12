#! /usr/bin/env python3
import argparse
import logging
import math
import os
import sys
import subprocess

from configurator_cades import configurator_cades
from configurator_nersc import configurator_nersc

sourcedir = os.path.dirname(os.path.abspath(sys.argv[0]))

def createMasterJobscript(workdir, outputdir, queue_configurator, nslots, minslot, nevents, ebeam, pdfset, muf, mur, bornkt, bornsupp, withnegweight):
    if not os.path.exists(workdir):
        os.makedirs(workdir, 0o755)
    jobscriptname = os.path.join(workdir, "jobscript.sh")
    logfile = os.path.join(workdir, "joboutput.log")
    with open(os.path.join(workdir, "jobscript.sh"), "w") as exewriter:
        exewriter.write("#! /bin/bash\n")
        queue_configurator.configure_batch(exewriter, logfile, nslots)
        queue_configurator.configure_modules(exewriter)
        queue_configurator.configure_workdir(exewriter)
        exewriter.write("echo \"Preparing working directories and configurations ...\"\n")
        for slot in range(0, nslots):
            slotdir = "%04d" %slot
            globalslot = slot + minslot
            exewriter.write("mkdir $WORKDIR/{}\n".format(slotdir))
            exewriter.write("cd $WORKDIR/{}\n".format(slotdir))
            configcreator = "{}/create_powheginput.py -n {} -s {} -e {} -p {} -f {} -r {} -k {} -b {}".format(sourcedir, nevents, globalslot, ebeam, pdfset, muf, mur, bornkt, bornsupp) 
            if withnegweight:
                configcreator += " -w"
            exewriter.write("python3 {}\n".format(configcreator))
        exewriter.write("cd $WORKDIR\n")
        exewriter.write("echo \"Running simulation on {} cores ... \"\n".format(nslots))
        exewriter.write("srun -n {} python3 {}/mpiwrapper_general.py {} {} \n".format(nslots, sourcedir, queue_configurator.get_name(), queue_configurator.get_image()))
        for slot in range(0, nslots):
            globalslot = slot + minslot
            slotdir = "%04d" %slot
            localdir = "%04d" %slot
            globaldir = "%04d" %globalslot
            exewriter.write("mv {} {}\n".format(localdir, globaldir))
            exewriter.write("zip -r {}.zip {}/ \n".format(globaldir, globaldir))
            exewriter.write("cp {}.zip {}\n".format(globaldir, outputdir))
        exewriter.write("cd {}\n".format(workdir))
        exewriter.write("rm -rf $WORKDIR\n")
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
    parser.add_argument("-q", "--qos", metavar="QOS", default="regular", help="Cori queue")
    parser.add_argument("-c", "--constraint", metavar="CONSTRAINT", default="knl", help="Platform (haswell or knl)")
    parser.add_argument("-t", "--time", metavar="TIME", default="10:00:00", help="Max. time")
    args = parser.parse_args()
    outputdir = os.path.abspath(args.outputdir)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir, 0o755)
    queue_configurator = None
    if os.getenv("NERSC_HOST") != None:
        queue_configurator = configurator_nersc(args.time, "docker:mfasel/cc7-alice:latest", args.constraint, args.queues)
    else:
        queue_configurator = configurator_cades(args.time, "/home/mfasel_alice/mfasel_cc7_alice.simg")
    slotspermaster = 0
    try:
        slotspermaster = queue_configurator.get_slots_per_master()
    except platform_exception as e:
        logging.error("%s", e)
        sys.exit(1)
    nmasters = math.floor(int(args.jobs) / slotspermaster)
    if args.jobs % slotspermaster != 0:
        nmasters += 1
    minslot = 0
    for masterID in range(0, nmasters):
        logging.info("Processing master %d", masterID)
        nslotsworker = args.jobs - minslot
        if nslotsworker > slotspermaster:
            nslotsworker = slotspermaster
        workdir = os.path.join(outputdir, "master%d" %masterID)
        jobscript = createMasterJobscript(workdir, outputdir, queue_configurator, nslotsworker, minslot, args.nevents, args.ebeam, args.pdfset, args.mufact, args.muren, args.bornkt, args.bornsupp, args.withnegweight)
        subprocess.call(["sbatch", jobscript])
        minslot += nslotsworker
