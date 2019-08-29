# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Job submition and convertion. Run a detached process and
    send messages to queue. This scheme doesn't freeze CAE
    while analysis is running or files are converting.
"""


import os, sys, logging, subprocess, queue
from PyQt5 import QtWidgets
try:
    from .ccx_log import logLine
except:
    from ccx_log import logLine


class Job:


    # Create job object
    def __init__(self, settings, file_name):
        self.settings = settings
        self.home_dir = os.path.dirname(os.path.abspath(sys.argv[0])) # app. home directory

        # Windows
        if os.name=='nt':
            self.op_sys = 'windows' # OS name
            self.extension = '.exe' # file extension in OS
            if self.home_dir.endswith('\\src'):
                self.home_dir = self.home_dir[:-4]

        # Linux
        else:
            self.op_sys = 'linux' # OS name
            self.extension = '' # file extension in OS
            if self.home_dir.endswith('/src'):
                self.home_dir = self.home_dir[:-4]

        self.path_ccx = os.path.join(self.home_dir, 'bin', 'ccx_2.15_MT') + self.extension
        if not len(file_name):
            file_name = 'job.inp'
        self.rename(file_name)


    # Rename job
    def rename(self, file_name):
        self.dir = os.path.abspath(os.path.dirname(file_name)) # working directory
        self.name = os.path.basename(file_name) # INP file name
        self.inp = os.path.abspath(file_name) # full path to INP file with extension
        self.path = self.inp[:-4] # full path to INP without extension
        self.frd = self.path + '.frd' # full path to job results file
        self.log = self.path + '.log' # full path to job log file
        self.unv = self.path + '.unv' # full path to imported UNV file

        # Log each job into file
        if len(logging.getLogger().handlers) > 1:
            logging.getLogger().handlers.pop()
        fh = logging.FileHandler(self.log, mode='w')
        fmt = logging.Formatter('%(levelname)s: %(message)s')
        fh.setFormatter(fmt)
        logging.getLogger().addHandler(fh)


    # Convert UNV to INP
    def importUNV(self):
        converter_path = os.path.join(self.home_dir, 'bin', 'unv2ccx' + self.extension)
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
            file_name = QtWidgets.QFileDialog.getOpenFileName(None, 'Open a subroutine',
                os.path.join(self.home_dir, 'ccx_' + self.op_sys,
                    'ccx_free_form_fortran'), 'FORTRAN (*.f)')[0]
            if file_name:
                # os.system(self.settings.path_editor + ' ' + file_name)
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
            home = '/cygdrive/' + \
                    self.home_dir[0].lower() + \
                    self.home_dir[2:].replace('\\', '/')

            # Path to ccx sources
            ccx = home + '/ccx_' + self.op_sys + '/ccx_free_form_fortran'

            # Open bash
            cmd1 = 'C:\\cygwin64\\bin\\bash.exe --login'
            # p = subprocess.Popen(cmd1, stdin=subprocess.PIPE,
            #         stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Send command to build CalculiX
            send1 = '/bin/make -f Makefile_MT -C {}\n'.format(ccx)
            # logging.info(send1)
            # p.stdin.write(bytes(send1, 'utf8'))
            # p.stdin.close()

            # # Read response: without it doesn't work
            # out = p.stdout.read()
            # for line in out.decode().split('\n'):
            #     logging.info(line)

            # Move built binary to ./bin
            cmd2 = 'C:\\cygwin64\\bin\\mv.exe ' + ccx + '/ccx_2.15_MT ' + \
                    home + '/bin/ccx_2.15_MT' + self.extension
            # logging.info(cmd2)
            # p = subprocess.Popen(cmd2, stdin=subprocess.PIPE,
            #         stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # logging.info('Compiled!')

            self.run([(cmd1, send1), (cmd2, '')], msg='Compiled!')

        # Linux
        else:

            # App home path
            home = self.home_dir

            make = 'make' # Make command
            move = 'mv' # Move command

            # Path to ccx sources
            ccx = home + '/ccx_' + self.op_sys + '/ccx_free_form_fortran'

            # Build CalculiX
            cmd1 = [make, '-f', 'Makefile_MT', '-C', ccx]

            # Move binary
            cmd2 = [move, ccx + '/ccx_2.15_MT', home + '/bin/ccx_2.15_MT' + self.extension]

            self.run([(cmd1, ''), (cmd2, '')], msg='Compiled!')


    # Submit INP to CalculiX
    def submit(self):
        if os.path.isfile(self.inp):
            os.environ['OMP_NUM_THREADS'] = str(os.cpu_count()) # enable multithreading
            cmd1 = [self.path_ccx, '-i', self.path]
            self.run([(cmd1, ''), ])
        else:
            logging.error('File not found: ' + self.inp)
            logging.error('Write input first.')


    # Open log file in external text editor
    def viewLog(self):
        if os.path.isfile(self.settings.path_editor):
            if os.path.isfile(self.log):
                # os.system(self.settings.path_editor + ' ' + self.log)
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
                    # command = self.settings.path_cgx + ' -o ' + self.frd
                    # os.system(command)
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
            converter_path = os.path.join(self.home_dir,
                'bin', 'ccx2paraview' + self.extension)
            cmd1 = [converter_path, self.frd, 'vtu']
            self.run([(cmd1, ''), ], msg='Finished!')
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

            # command = self.settings.path_paraview + ' --data=' + vtu_path
            # os.system(command)
            command = [self.settings.path_paraview, '--data=' + vtu_path]
            subprocess.Popen(command)
        else:
            logging.error('Wrong path to Paraview: ' + \
                self.settings.path_paraview +\
                '. Configure it in File->Settings.')


    # Run multiple commands and log stdout without blocking GUI
    def run(self, commands, msg=None):
        import time, threading

        start = time.perf_counter() # start time
        os.chdir(self.dir)

        for cmd1, cmd2 in commands:
            logging.debug(str(cmd1) + ' ' + cmd2)

            # Run command: works both in Linux and in Windows
            p = subprocess.Popen(cmd1, stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
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


# Tests
if __name__ == '__main__':

    # Configure logging
    logging.log = print
    logging.debug = print
    logging.info = print
    logging.warning = print

    # Rebuild CCX
    Job(None, '').rebuildCCX()

    # Remove cached files
    import clean
    clean.cleanCache()
