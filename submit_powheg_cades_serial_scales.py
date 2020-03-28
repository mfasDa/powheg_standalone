#! /usr/bin/env python3
import argparse
import logging
import math
import os
import sys
import subprocess

sourcedir = os.path.dirname(os.path.abspath(sys.argv[0]))

def createJobscript(workdir, maxtime):
    jobscriptname = os.path.join(workdir, "jobscript_powheg_variations.sh")
    logfile = os.path.join(workdir, "joboutput_powheg_variations.log")
    variations = {"rlfl": {"renscale": 0.5, "factscale": 0.5}, 
                  "rlfm": {"renscale": 0.5, "factscale": 1.},
                  "rmfl": {"renscale": 1., "factscale": 0.5},
                  "rhfm": {"renscale": 2., "factscale": 1.},
                  "rmfh": {"renscale": 1., "factscale": 2.},
                  "rhfh": {"renscale": 2., "factscale": 2.}}
    with open(jobscriptname, "w") as exewriter:
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
        for tag,scales in variations.items():
            configcreator = "{}/createScaleVariation.py -i powheg.input.main -t {} -r {} -f {}".format(sourcedir, tag, scales["renscale"], scales["factscale"]) 
            exewriter.write("python3 {}\n".format(configcreator))
            exewriter.write("echo \"Running simulation {} (muf {}, mur{}) ... \"\n".format(tag, scales["factscale"], scales["renscale"]))
            exewriter.write("singularity exec -B /nfs/home:/nfs/home -B /lustre:/lustre /home/mfasel_alice/mfasel_cc7_alice.simg {}/run_powheg_general.sh {} &> run_powheg_{}.log\n".format(sourcedir, tag, tag))
            exewriter.write("mv powheg.input powheg.input.scale.{}\n".format(tag))
            exewriter.write("mv pwgevents-rwgt.lhe pwgevents.lhe\n")
        exewriter.write("rm -v {}\n".format(jobscriptname))
        exewriter.write("echo Job done\n")
        exewriter.close() 
    return jobscriptname

if __name__ == "__main__":
    logging.basicConfig(format='[%(levelname)s]: %(message)s', level=logging.INFO)
    parser = argparse.ArgumentParser("submit_powheg.py", "Submitter for POWHEG dijet process")
    parser.add_argument("-i", "--inputdir", metavar="INPUTDIR", type=str, required=True, help="Output directory")
    parser.add_argument("-t", "--time", metavar="TIME", default="10:00:00", help="Max. time")
    args = parser.parse_args()
    for chunk in os.listdir(args.inputdir):
        logging.info("Processing job directory %s", chunk)
        jobscript = createJobscript(os.path.join(args.inputdir, chunk), args.time)
        subprocess.call(["sbatch", jobscript])
