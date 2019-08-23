# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Job submition and convertion. Run a detached process and
    send messages to queue. This scheme doesn't freeze CAE
    while analysis is running or files are converting.
"""


import os, sys, logging, subprocess, time, ccx_log, queue, threading
import multiprocessing as mp
from PyQt5 import QtWidgets


class Job:


    # Create job object
    def __init__(self, settings, file_name):
        self.settings = settings
        if not len(file_name):
            file_name = 'job.inp'
        self.rename(file_name)
        self.home_dir = os.path.dirname(sys.argv[0]) # app. home directory
        self.extension = '.exe' if os.name=='nt' else '' # file extension in OS


    # Rename job
    def rename(self, file_name):
        self.dir = os.path.dirname(file_name) # working directory
        self.name = os.path.basename(file_name) # INP file name
        self.inp = os.path.abspath(file_name) # full path to INP file with extension
        self.path = self.inp[:-4] # full path to INP without extension
        self.frd = self.path + '.frd' # full path to job results file
        self.log = self.path + '.log' # full path to job log file
        self.unv = self.path + '.unv' # full path to imported UNV file


    # Convert UNV to INP
    def importUNV(self):
        converter_path = os.path.join(self.home_dir,
            'converters', 'unv2ccx' + self.extension)
        command = [converter_path, self.unv]
        self.run(command)


    # Submit INP to CalculiX
    def submit(self):
        if os.path.isfile(self.settings.path_ccx):
            os.environ['OMP_NUM_THREADS'] = str(mp.cpu_count()) # enable multithreading
            command = [self.settings.path_ccx, '-i', self.path]
            self.run(command)
        else:
            logging.error('Wrong path to CCX: ' + self.settings.path_ccx)


    # Open FRD in GraphiX
    def openCGX(self):
        if os.path.isfile(self.frd):
            if os.path.isfile(self.settings.path_cgx):
                command = [self.settings.path_cgx, '-o', self.frd]
                self.run(command)
            else:
                logging.error('Wrong path to CGX: ' + self.settings.path_cgx)
        else:
            logging.error('File not found: ' + self.frd)


    # Convert FRD to VTU
    def exportVTU(self):
        if os.path.isfile(self.frd):
            if os.path.isfile(self.settings.path_paraview):
                converter_path = os.path.join(self.home_dir,
                    'converters', 'ccx2paraview' + self.extension)
                command = [converter_path, self.frd, 'vtu']
                self.run(command)
            else:
                logging.error('Wrong path to Paraview: ' + self.settings.path_paraview)
        else:
            logging.error('File not found: ' + self.frd)


    # Open VTU in Paraview
    def openParaview(self):

        # Count result VTU files
        file_list = []
        for f in os.listdir(self.dir):
            f = os.path.basename(f)
            if f.lower() == self.name[:-4] + '.vtu':
                file_list = []
                break
            if f.lower().endswith('.vtu') and f.startswith(self.name[:-4]):
                file_list.append(f)
        if len(file_list) > 1:
            vtu_path = self.path + '...vtu'
        else:
            vtu_path = self.path + '.vtu'

        if os.path.isfile(self.settings.path_paraview):
            command = [self.settings.path_paraview, '--data=' + vtu_path]
            self.run(command)
        else:
            logging.error('Wrong path to Paraview: ' + self.settings.path_paraview)


    # Run commands and log stdout without blocking GUI
    def run(self, command):
        os.chdir(self.dir)
        start = time.perf_counter() # start time

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
