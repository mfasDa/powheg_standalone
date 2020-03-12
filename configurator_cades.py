#! /usr/bin/env python3

import subprocess

class configurator_cades:

    def __init__(self, jobtime, image):
        self.__jobtime = jobtime
        self.__image = image

    def configure_batch(self, jobscript, logfile, nslots):
        jobscript.write("#SBATCH -N 1\n")
        jobscript.write("#SBATCH -A birthright\n")
        jobscript.write("#SBATCH --ntasks-per-node={}\n".format(nslots))
        jobscript.write("#SBATCH --partition=gpu\n")
        jobscript.write("#SBATCH --mem=32G\n")
        jobscript.write("#SBATCH -J powheg\n")
        jobscript.write("#SBATCH -o {}\n".format(logfile))
        jobscript.write("#SBATCH -t {}\n".format(self.__jobtime))

    def configure_workdir(self, jobscript):
        jobscript.write("SCRATCH=/lustre/or-hydra/cades-birthright/mfasel_alice\n")
        jobscript.write("WORKDIR=$SCRATCH/spool/$SLURM_JOBID\n")
        jobscript.write("echo \"Using working directory $WORKDIR\"\n")
        jobscript.write("if [ ! -d $WORKDIR ]; then mkdir -p $WORKDIR; fi\n")
        jobscript.write("cd $WORKDIR\n")

    def configure_modules(self, jobscript):
        jobscript.write("module load python/3.6.3\n")
        jobscript.write("module load PE-gnu\n")
        jobscript.write("module load singularity\n")

    def run_image(self, taskscript, logfile):
        subprocess.call("singularity exec -B /home:/home -B /lustre:/lustre exec {} {} &> {}".format(self.__image, taskscript, logfile), shell = True) 

    def get_slots_per_master(self):
        return 32

    def get_image(self):
        return self.__image

    @staticmethod
    def get_default_image():
        return "/home/mfasel_alice/mfasel_cc7_alice.simg"

    def get_name(self):
        return "cades"
        