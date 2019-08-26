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
        self.home_dir = os.path.dirname(os.path.abspath(sys.argv[0])) # app. home directory
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
            'bin', 'unv2ccx' + self.extension)
        command = [converter_path, self.unv]
        self.run(command)


    # Open INP file in external text editor
    def editINP(self):
        if os.path.isfile(self.inp):
            os.system(self.settings.path_editor + ' ' + self.inp)
        else:
            logging.error('File not found: ' + self.inp)
            logging.error('Write input first.')


    # Dialog window to filter fortran subroutines
    def openSubroutine(self):
        file_name = QtWidgets.QFileDialog.getOpenFileName(None, \
            'Open a subroutine', \
            os.path.join(self.home_dir, 'src', 'ccx_free_form_fortran'), \
            'FORTRAN (*.f)')[0]
        if file_name:
            os.system(self.settings.path_editor + ' ' + file_name)


    # Recompile CalculiX sources with updated subroutines
    def rebuildCCX(self):
        if os.name=='nt':
            logging.warning('CalculiX rebuild not jet supported in fuck\'n Windows!')
            logging.warning('Please, format C: and install normal OS like Ubuntu.')
        else:
            p0 = os.path.join('src', 'ccx_free_form_fortran')
            p1 = os.path.join(p0, 'ccx_2.15_MT')
            p2 = os.path.join('bin', 'ccx')
            subprocess.run('make -f Makefile_MT -C ' + p0, shell=True)
            subprocess.run('cp ' + p1 + ' ' + p2, shell=True)

            # TODO Use self.run()
            # command = ['make', '-f', 'Makefile_MT', '-C', p0, ';cp', p1, p2]
            # self.run(command)

            # TODO Implement event notification to make logger wait for completion
            logging.info('Compilation OK!')


    # Submit INP to CalculiX
    def submit(self):
        if os.path.isfile(self.inp):
            if os.path.isfile(self.settings.path_ccx):
                os.environ['OMP_NUM_THREADS'] = str(mp.cpu_count()) # enable multithreading
                command = [self.settings.path_ccx, '-i', self.path]
                self.run(command)
            else:
                logging.error('Wrong path to CCX: ' + self.settings.path_ccx)
        else:
            logging.error('File not found: ' + self.inp)
            logging.error('Write input first.')


    # Open log file in external text editor
    def viewLog(self):
        if os.path.isfile(self.log):
            os.system(self.settings.path_editor + ' ' + self.log)
        else:
            logging.error('File not found: ' + self.log)
            logging.error('Submit analysis first.')


    # Open FRD in GraphiX
    def openCGX(self):
        if os.path.isfile(self.frd):
            if os.path.isfile(self.settings.path_cgx):
                command = self.settings.path_cgx + ' -o ' + self.frd
                os.system(command)
            else:
                logging.error('Wrong path to CGX: ' + self.settings.path_cgx)
        else:
            logging.error('File not found: ' + self.frd)
            logging.error('Submit analysis first.')


    # Convert FRD to VTU
    def exportVTU(self):
        if os.path.isfile(self.frd):
            if os.path.isfile(self.settings.path_paraview):
                converter_path = os.path.join(self.home_dir,
                    'bin', 'ccx2paraview' + self.extension)
                command = [converter_path, self.frd, 'vtu']
                self.run(command)
            else:
                logging.error('Wrong path to Paraview: ' + self.settings.path_paraview)
        else:
            logging.error('File not found: ' + self.frd)
            logging.error('Submit analysis first.')


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
        elif len(file_list) == 1:
            vtu_path = self.path + '.vtu'
        else:
            logging.error('VTU file not found.')
            logging.error('Export VTU results first.')
            return

        if os.path.isfile(self.settings.path_paraview):
            command = self.settings.path_paraview + ' --data=' + vtu_path
            os.system(command)
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
                    time.sleep(0.1) # reduce CPU usage

        os.chdir(self.home_dir)


# Put stdout lines to queue
def enqueue_output(stdout, queue):
    for line in iter(stdout.readline, b''):
        line = line.decode().strip()
        queue.put(line)
    queue.put('END') # mark to break while loop
    stdout.close()


# Test ccx sources recompilation
if __name__ == '__main__':
    import clean
    Job(None, '').rebuildCCX()
    clean.cache()
