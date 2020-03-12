#! /usr/bin/env python3

import subprocess

class configurator_cades:

    def __init__(self, jobtime, image, platform, queue):
        self.__jobtime = jobtime
        self.__image = image
        self.__platform = platform
        self.__queue = queue

    def configure_batch(self, jobscript, logfile, nslots):
        jobscript.write("#SBATCH -N 1\n")
        jobscript.write("#SBATCH -A birthright\n")
        jobscript.write("#SBATCH --ntasks-per-node={}\n".format(nslots))
        jobscript.write("#SBATCH --platform=gpu\n")
        jobscript.write("#SBATCH -J powheg\n")
        jobscript.write("#SBATCH -o {}\n".format(logfile))
        jobscript.write("#SBATCH -t {}\n".format(self.__jobtime))

    def configure_spool(self, jobscript):
        exewriter.write("WORKDIR=$SCRATCH/spool/$SLURM_JOBID\n")
        exewriter.write("echo \"Using working directory $WORKDIR\"\n")
        exewriter.write("if [ ! -d $WORKDIR ]; then mkdir -p $WORKDIR; fi\n")
        exewriter.write("cd $WORKDIR\n")

    def configure_modules(self, jobscript):
        jobscript.write("module load PE-gnu\n")
        jobscript.write("module load singularity")

    def run_image(self, taskscript, logfile):
        subprocess.call("singularity exec -B /home:/home -B /lustre:/lustre exec {} {} &> {}".format(self.__image, taskscript, logfile), shell = True) 

    def get_slots_per_master(self):
        return 32

    def get_image(self):
        return self.__image

    @staticmethod
    def get_default_image()
        return "/home/mfasel_alice/mfasel_cc7_alice.simg"

    def get_name():
        return "cades"
        