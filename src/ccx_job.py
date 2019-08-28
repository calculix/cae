# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Job submition and convertion. Run a detached process and
    send messages to queue. This scheme doesn't freeze CAE
    while analysis is running or files are converting.
"""


import os, sys, logging, queue
from PyQt5 import QtWidgets
try:
    from .ccx_log import logLine
except:
    from ccx_log import logLine


class Job:


    # Create job object
    def __init__(self, settings, file_name):
        self.settings = settings
        if not len(file_name):
            file_name = 'job.inp'
        self.rename(file_name)
        self.home_dir = os.path.dirname(os.path.abspath(sys.argv[0])) # app. home directory
        self.extension = '.exe' if os.name=='nt' else '' # file extension in OS
        self.path_ccx = os.path.join(self.home_dir, 'bin', 'ccx_2.15_MT') + self.extension


    # Rename job
    def rename(self, file_name):
        self.dir = os.path.abspath(os.path.dirname(file_name)) # working directory
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
        self.run([command])


    # Open INP file in external text editor
    def editINP(self):
        if os.path.isfile(self.settings.path_editor):
            if os.path.isfile(self.inp):
                os.system(self.settings.path_editor + ' ' + self.inp)
            else:
                logging.error('File not found: ' + self.inp)
                logging.error('Write input first.')
        else:
            logging.error('Wrong path to text editor: ' + \
                self.settings.path_editor +\
                '. Configure it in File->Settings.')


    # Dialog window to filter fortran subroutines
    def openSubroutine(self):
        if os.path.isfile(self.settings.path_editor):
            file_name = QtWidgets.QFileDialog.getOpenFileName(None, \
                'Open a subroutine', \
                os.path.join(self.home_dir, 'ccx', 'ccx_free_form_fortran'), \
                'FORTRAN (*.f)')[0]
            if file_name:
                os.system(self.settings.path_editor + ' ' + file_name)
        else:
            logging.error('Wrong path to text editor: ' + \
                self.settings.path_editor +\
                '. Configure it in File->Settings.')


    # Recompile CalculiX sources with updated subroutines
    def rebuildCCX(self):

        # Windows
        if os.name=='nt':

            # App home path in cygwin
            home = '/cygdrive/' + \
                    self.home_dir[0].lower() + \
                    self.home_dir[2:].replace('\\', '/')

            # Path to ccx sources in cygwin
            ccx = home + '/ccx_windows/ccx_free_form_fortran'
           
            make = 'C:\\cygwin64\\bin\\make.exe' # Make command
            move = 'C:\\cygwin64\\bin\\mv.exe' # Move command

        # Linux
        else:

            # App home path
            home = self.home_dir

            # Path to ccx sources
            ccx = home + '/ccx_linux/ccx_free_form_fortran'

            make = 'make' # Make command
            move = 'mv' # Move command

        # Build CalculiX
        command1 = [make, '-f', 'Makefile_MT', '-C', ccx]

        # Move binary
        command2 = [move, ccx + '/ccx_2.15_MT', home + '/bin/ccx_2.15_MT' + self.extension]

        self.run([command1, command2], msg='Compiled!')


    # Submit INP to CalculiX
    def submit(self):
        if os.path.isfile(self.inp):
            os.environ['OMP_NUM_THREADS'] = str(os.cpu_count()) # enable multithreading
            command = [self.path_ccx, '-i', self.path]
            self.run([command])
        else:
            logging.error('File not found: ' + self.inp)
            logging.error('Write input first.')


    # Open log file in external text editor
    def viewLog(self):
        if os.path.isfile(self.settings.path_editor):
            if os.path.isfile(self.log):
                os.system(self.settings.path_editor + ' ' + self.log)
            else:
                logging.error('File not found: ' + self.log)
                logging.error('Submit analysis first.')
        else:
            logging.error('Wrong path to text editor: ' + \
                self.settings.path_editor +\
                '. Configure it in File->Settings.')


    # Open FRD in GraphiX
    def openCGX(self):
        if os.path.isfile(self.settings.path_cgx):
            if os.path.isfile(self.frd):
                    command = self.settings.path_cgx + ' -o ' + self.frd
                    os.system(command)
            else:
                logging.error('File not found: ' + self.frd)
                logging.error('Submit analysis first.')
        else:
            logging.error('Wrong path to CGX: ' + \
                self.settings.path_cgx +\
                '. Configure it in File->Settings.')


    # Convert FRD to VTU
    def exportVTU(self):
        if os.path.isfile(self.frd):
            converter_path = os.path.join(self.home_dir,
                'bin', 'ccx2paraview' + self.extension)
            command = [converter_path, self.frd, 'vtu']
            self.run([command])
        else:
            logging.error('File not found: ' + self.frd)
            logging.error('Submit analysis first.')


    # Open VTU in Paraview
    def openParaview(self):
        if os.path.isfile(self.settings.path_paraview):

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

            command = [self.settings.path_paraview, '--data=' + vtu_path]
            self.run([command])
        else:
            logging.error('Wrong path to Paraview: ' + \
                self.settings.path_paraview +\
                '. Configure it in File->Settings.')


    # Run commands and log stdout without blocking GUI
    def run(self, commands, msg=None):
        import subprocess, time, threading

        start = time.perf_counter() # start time
        os.chdir(self.dir)

        for command in commands:
            logging.info(' '.join(command))

            # Run command
            p = subprocess.Popen(command,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            # Make daemon to enqueue stdout
            q = queue.Queue()
            t = threading.Thread(target=enqueue_output, args=(p.stdout, q))
            t.daemon = True # thread dies with the program
            t.start()

            # Read and log stdout without blocking GUI
            with open(self.log, 'a') as lf:
                while True:
                    try:
                        line = q.get_nowait()
                        if line == 'END':
                            logging.info('Total {:.0f} seconds.'\
                                .format(time.perf_counter()-start)) # end time
                            break
                        logLine(line)
                        lf.write(line + '\n')
                    except queue.Empty:
                        QtWidgets.qApp.processEvents() # do not block GUI
                        time.sleep(0.1) # reduce CPU usage

        if msg:
            logging.info(msg)
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
    clean.cleanCache()
