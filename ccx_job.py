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
    def __init__(self, settings, path):
        self.settings = settings
        self.name = self.settings.job_prefix
        self.path = '.'
        self.rename(path)


    # Rename job
    def rename(self, path):

        # Full path
        self.path = os.path.abspath(path)

        # Job name with prefix
        self.name = os.path.basename(path)
        if not self.name.startswith(self.settings.job_prefix):
            self.name = self.settings.job_prefix + self.name

        # Full path to job with prefix
        self.path = os.path.join(os.path.dirname(path), self.name)


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
        command = '{0} -i {1} > {1}.log'\
                    .format(self.settings.path_ccx, self.path[:-4])
        subprocess.run(command, shell=True)

        # Notify about job completeion
        job_completed.set()
