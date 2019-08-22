# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Job submition and convertion. Run a detached process and
    send messages to queue. This scheme doesn't freeze CAE
    while analysis is running or files are converting.
"""


import os, logging, subprocess, time, ccx_log, queue, threading
import multiprocessing as mp
from PyQt5 import QtWidgets


class Job:


    # Create job object
    def __init__(self, settings, file_name):
        self.settings = settings
        self.rename(file_name)
        self.home_dir = os.getcwd() # installation directory


    # Rename job
    def rename(self, file_name):
        self.dir = os.path.dirname(file_name) # working directory
        self.name = os.path.basename(file_name) # INP file name
        self.inp = os.path.abspath(file_name) # full path to INP file with extension
        self.path = self.inp[:-4] # full path to INP without extension
        self.frd = self.path + '.frd' # full path to job results file
        self.log = self.path + '.log' # full path to job log file
        self.unv = self.path + '.unv' # full path to imported UNV file


    # Submit job for CalculiX
    def submit(self):
        if os.path.isfile(self.settings.path_ccx):

            # Enable multithreading with needed amount of cores
            os.environ['OMP_NUM_THREADS'] = str(mp.cpu_count())

            command = [self.settings.path_ccx, '-i', self.path]
            self.run(command)
        else:
            logging.error('Wrong path to CCX: ' + self.settings.path_ccx)


    # Convert FRD to VTU
    def exportVTU(self):
        if os.path.isfile(self.frd):
            extension = '.exe' if os.name=='nt' else '' # file extension in OS
            converter_path = \
                os.path.join(self.home_dir, 'converters',
                    'ccx2paraview' + extension)
            command = [converter_path, self.frd, 'vtu']
            self.run(command)
        else:
            logging.error('File not found: ' + self.frd)


    def importUNV(self):
        extension = '.exe' if os.name=='nt' else '' # file extension in OS
        converter_path = os.path.join(self.home_dir,
                            'converters', 'unv2ccx' + extension)
        command = [converter_path, self.unv]
        self.run(command)


    # Run job command and log stdout without blocking GUI
    def run(self, command):
        start = time.perf_counter() # start time
        os.chdir(self.dir)

        # Run command
        p = subprocess.Popen(command,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        # Make daemon to enqueue stdout
        q = queue.Queue()
        t = threading.Thread(target=enqueue_output, args=(p.stdout, q))
        t.daemon = True # thread dies with the program
        t.start()

        # Read and log stdout without blocking GUI
        with open(self.log, 'w') as lf:
            while True:
                try:
                    line = q.get_nowait()
                    if line == 'END':
                        logging.info('Total {:.0f} seconds.'\
                            .format(time.perf_counter()-start)) # end time
                        break
                    ccx_log.logLine(line)
                    lf.write(line + '\n')
                except queue.Empty:
                    QtWidgets.qApp.processEvents() # do not block GUI

        os.chdir(self.home_dir)


# Put stdout lines to queue
def enqueue_output(stdout, queue):
    for line in iter(stdout.readline, b''):
        line = line.decode().strip()
        queue.put(line)
    queue.put('END') # mark to break while loop
    stdout.close()
