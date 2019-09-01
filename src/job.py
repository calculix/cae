# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, September 2019
    Distributed under GNU General Public License v3.0

    Job submition and convertion. Run a detached process and
    send messages to queue. This scheme doesn't freeze CAE
    while analysis is running or files are converting.
"""


from path import Path
import os, sys, logging, subprocess, queue
from PyQt5 import QtWidgets
from log import logLine
from settings import Settings


class Job:


    # Create job object
    def __init__(self, settings, file_name):
        self.settings = settings
        self.p = Path()
        logging.info('Application\'s home directory is: ' + self.p.app_home_dir)

        # Windows
        if os.name=='nt':
            self.op_sys = 'windows' # OS name
            self.extension = '.exe' # file extension in OS

        # Linux
        else:
            self.op_sys = 'linux' # OS name
            self.extension = '' # file extension in OS

        if not len(file_name):
            file_name = settings.path_start_model
        self.rename(file_name)


    # Rename job
    def rename(self, file_name):
        self.dir = os.path.dirname(os.path.abspath(file_name)) # working directory
        logging.info('Work directory is: ' + self.dir)
        self.name = os.path.basename(file_name) # INP file name
        self.inp = os.path.abspath(file_name) # full path to INP file with extension
        self.path = self.inp[:-4] # full path to INP without extension
        self.frd = self.path + '.frd' # full path to job results file
        self.log = self.path + '.log' # full path to job log file
        self.unv = self.path + '.unv' # full path to imported UNV file

        # Log each job into file
        if len(logging.getLogger().handlers) > 1:
            logging.getLogger().handlers.pop()
        fh = logging.FileHandler(self.log, mode='a') # TODO with 'w' log is empty after saving INP
        fmt = logging.Formatter('%(module)s, %(levelname)s: %(message)s')
        fh.setFormatter(fmt)
        logging.getLogger().addHandler(fh)


    # Convert UNV to INP
    def importUNV(self):
        converter_path = os.path.join(self.p.bin, 'unv2ccx' + self.extension)
        # cmd1 = converter_path + ' ' + self.unv
        cmd1 = [converter_path, self.unv]
        self.run([(cmd1, ''), ])


    # Open INP file in external text editor
    def editINP(self):
        if os.path.isfile(self.settings.path_editor):
            if os.path.isfile(self.inp):
                # os.system(self.settings.path_editor + ' ' + self.inp)
                command = [self.settings.path_editor, self.inp]
                subprocess.Popen(command)
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
            file_name = QtWidgets.QFileDialog.getOpenFileName(None,
                'Open a subroutine', self.p.ccx, 'FORTRAN (*.f)')[0]
            if file_name:
                command = [self.settings.path_editor, file_name]
                subprocess.Popen(command)
        else:
            logging.error('Wrong path to text editor: ' + \
                self.settings.path_editor +\
                '. Configure it in File->Settings.')


    # Recompile CalculiX sources with updated subroutines
    def rebuildCCX(self):

        # Windows
        if os.name=='nt':

            # App home path in cygwin
            home = path2cygwin(self.home_dir)

            # Path to ccx sources
            ccx = home + '/ccx_' + self.op_sys + '/ccx_free_form_fortran'

            # Open bash
            cmd1 = 'C:\\cygwin64\\bin\\bash.exe --login'

            # Send command to build CalculiX
            send1 = '/bin/make -f Makefile_MT -C {}'.format(ccx)

            # Move built binary
            cmd2 = 'C:\\cygwin64\\bin\\mv.exe ' + ccx + '/ccx_2.15_MT ' + \
                    path2cygwin(self.p.ccx)

            self.run([(cmd1, send1), (cmd2, '')], msg='Compiled!')

        # Linux
        else:

            # Path to ccx sources
            ccx = self.home_dir + '/ccx_' + self.op_sys + '/ccx_free_form_fortran'

            # Build CalculiX
            cmd1 = ['make', '-f', 'Makefile_MT', '-C', ccx]

            # Move binary
            cmd2 = ['mv', ccx + '/ccx_2.15_MT', self.p.ccx]

            self.run([(cmd1, ''), (cmd2, '')], msg='Compiled!')


    # Submit INP to CalculiX
    def submit(self):
        if os.path.isfile(self.settings.path_ccx):
            if os.path.isfile(self.inp):
                os.environ['OMP_NUM_THREADS'] = str(os.cpu_count()) # enable multithreading
                cmd1 = [self.p.ccx, '-i', self.path]
                self.run([(cmd1, ''), ])
            else:
                logging.error('File not found: ' + self.inp)
                logging.error('Write input first.')
        else:
            logging.error('Wrong path to CCX: ' + \
                self.settings.path_ccx +\
                '. Configure it in File->Settings.')


    # Open log file in external text editor
    def viewLog(self):
        if os.path.isfile(self.settings.path_editor):
            if os.path.isfile(self.log):
                command = [self.settings.path_editor, self.log]
                subprocess.Popen(command)
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
                    command = [self.settings.path_cgx, '-o', self.frd]
                    subprocess.Popen(command)
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
            converter_path = os.path.join(self.p.bin,
                    'ccx2paraview' + self.extension)
            cmd1 = [converter_path, self.frd, 'vtu']
            self.run([(cmd1, ''), ], msg='Finished!')
        else:
            logging.error('File not found: ' + self.frd)
            logging.error('Submit analysis first.')


    # Open VTU in ParaView
    def openParaView(self):
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

            # command = self.settings.path_paraview + ' --data=' + vtu_path
            # os.system(command)
            command = [self.settings.path_paraview, '--data=' + vtu_path]
            subprocess.Popen(command)
        else:
            logging.error('Wrong path to ParaView: ' + \
                self.settings.path_paraview +\
                '. Configure it in File->Settings.')


    # Run multiple commands and log stdout without blocking GUI
    def run(self, commands, msg=None):
        import time, threading

        start = time.perf_counter() # start time
        os.chdir(self.dir)

        for cmd1, cmd2 in commands:
            logging.info(str(cmd1) + ' ' + cmd2)

            # Run command: works both in Linux and in Windows
            p = subprocess.Popen(cmd1, stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if len(cmd2):
                p.stdin.write(bytes(cmd2, 'utf8'))
                p.stdin.close()

            # Make daemon to enqueue stdout
            q = queue.Queue()
            t = threading.Thread(target=enqueue_output, args=(p.stdout, q))
            t.daemon = True # thread dies with the program
            t.start()

            # Read and log stdout without blocking GUI
            while True:
                try:
                    line = q.get_nowait()
                    if line == 'END':
                        break
                    logLine(line)
                except queue.Empty:
                    QtWidgets.qApp.processEvents() # do not block GUI
                    time.sleep(0.1) # reduce CPU usage

        os.chdir(self.home_dir)

        # Total time passed
        total_time = 'Total {:.0f} seconds.'.format(time.perf_counter() - start)
        logging.info(total_time)

        # Post final message
        if msg:
            logging.info(msg)


# Put stdout lines to queue
def enqueue_output(stdout, queue):
    for line in iter(stdout.readline, b''):
        line = line.decode().strip()
        queue.put(line)
    queue.put('END') # mark to break while loop
    stdout.close()


# Converts Windows path to Cygwin path
def path2cygwin(path):
    return '/cygdrive/' + \
            path[0].lower() + \
            path[2:].replace('\\', '/')


# Tests
if __name__ == '__main__':

    # Configure logging
    logging.log = print
    logging.debug = print
    logging.info = print
    logging.warning = print

    settings = Settings()

    # Create job
    j = Job(settings, '')

    # Rebuild CCX
    j.rebuildCCX()

    # Remove cached files
    import clean
    clean.cleanCache()
