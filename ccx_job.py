# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Job to submit from CAE
"""


import os, logging, subprocess, time
import multiprocessing as mp
from PyQt5 import QtWidgets


class Job:


    # Create job object
    def __init__(self, settings, file_name):
        self.settings = settings
        self.rename(file_name)


    # Rename job
    def rename(self, file_name):
        self.path = os.path.abspath(file_name) # full path to INP file
        self.dir = os.path.dirname(file_name) # working directory
        self.name = os.path.basename(file_name) # INP file name


    # Submit job
    def submit(self):
        if os.path.isfile(self.settings.path_ccx):
            logging.info('Job submitted.')
            job_completed = mp.Event()

            # Run analysis in a detached process
            job = mp.Process(target=self.runJob, args=(job_completed,))
            job.start()

            # Wait for job completion
            while not job_completed.is_set():
                QtWidgets.qApp.processEvents() # do not freeze GUI
                time.sleep(0.1)

            # Terminate job process after finishing the job
            job.terminate()

            logging.info("Job completed.")
        else:
            logging.error('Wrong path to CCX: ' + self.settings.path_ccx)


    # Run analysis in a detached process
    def runJob(self, job_completed):

        # Enable multithreading with needed amount of cores
        os.environ['OMP_NUM_THREADS'] = str(mp.cpu_count())

        # Run job in CalculiX + write log
        os.chdir(self.dir)
        command = '{0} -i {1} > {1}.log'\
                    .format(self.settings.path_ccx, self.path[:-4])
        subprocess.run(command, shell=True)

        # Notify about job completeion
        job_completed.set()
