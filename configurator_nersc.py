#! /usr/bin/env python3
import subprocess

class platform_exception(Exception):

    def __init__(self, platform):
        self.__platform = platform

    def __str__(self):
        return "Unknown platform {}, select either \"haswell\" or \"knl\"".format(self.__platform)

class configurator_nersc:

    def __init__(self, jobtime, image, platform, queue):
        self.__jobtime = jobtime
        self.__image = image
        self.__platform = platform
        self.__queue = queue

    def configure_batch(self, jobscript, logfile, nslots):
        jobscript.write("#SBATCH -N 1\n")
        jobscript.write("#SBATCH --ntasks-per-node={}\n".format(nslots))
        jobscript.write("#SBATCH --qos={}\n".format(self.__queue))
        jobscript.write("#SBATCH -J powheg\n")
        jobscript.write("#SBATCH -o {}\n".format(logfile))
        jobscript.write("#SBATCH -t {}\n".format(self.__jobtime))
        jobscript.write("#SBATCH -C {}\n".format(self.__platform))
        jobscript.write("#SBATCH --image={}\n".format(self.__image))

    def configure_spool(self, jobscript):
        jobscript.write("export SCRATCH=/lustre/or-hydra/cades-birthright/mfasel_alice\n")
        jobscript.write("WORKDIR=$SCRATCH/spool/$SLURM_JOBID\n")
        jobscript.write("echo \"Using working directory $WORKDIR\"\n")
        jobscript.write("if [ ! -d $WORKDIR ]; then mkdir -p $WORKDIR; fi\n")
        jobscript.write("cd $WORKDIR\n")

    def configure_modules(self, jobscript):
        jobscript.write("module load cray-python/3.7.3.2\n"))

    def run_image(self, taskscript, logfile):
        subprocess.call("shifter {} &> {}".format(taskscript, logfile), shell = True)

    def get_slots_per_master(self):
        if self.__platform = "knl":
            return 68
        elif self.__platform == "haswell":
            return 32
        raise platform_exception(self.__platform)

    def get_image(self):
        return self.__image

    @staticmethod
    def get_default_image():
        return "docker:mfasel/cc7-alice:latest"

    def get_name():
        return "cori"