#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Class Job represents a container for all jobs from the Job menu.

Job instance creation indicates start of the new session.
Old log file is erased and a new one is created.

If needed job method is called via a thread and method 'self.run'.
Terminal command is passed to the thread with self.run in this case.

# TODO All stuff with threads should be concentrated here.
"""

# Standard modules
import os
import sys
import time
import logging
import subprocess
import threading

# External modules
from PyQt5 import QtWidgets

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from path import p
from settings import s
from model import m
import log


def list_threads():
    """List currently running threads."""
    t_names = sorted([t.name for t in threading.enumerate() \
        if t.name != threading.main_thread().name])
    msg = '\nRunning threads:\n' + '\n'.join(t_names) + '\n'
    logging.debug(msg)


thread_counter = 0


class Job:

    def __init__(self):
        pass

    def generate(self, file_name=s.start_model):
        """Create job object. Is called from importer.py.
        NOTE file_name could be None
        """
        log.print_to_file(p.log, '\nCREATING JOB INSTANCE\n')

        self.dir = os.path.dirname(os.path.abspath(file_name)) # working directory
        self.name = os.path.basename(file_name) # INP file name
        self.inp = os.path.abspath(file_name) # full path to INP file with extension
        self.path = self.inp[:-4] # full path to INP without extension
        self.frd = self.path + '.frd' # full path to job results file
        self.log = self.path + '.log' # full path to job log file
        self.sta = self.path + '.sta' # full path to job status file

        # Remove old log_file
        if os.path.exists(self.log):
            os.remove(self.log)

        logging.info('Application home directory is: '\
            + p.app_home_dir)
        os.chdir(self.dir)
        logging.info('Work directory is: ' + self.dir)

    def convert_unv(self):
        """Convert UNV to INP."""
        import unv2ccx
        unv2ccx.Converter(self.path + '.unv').run()

    def monitor_status(self):
        """Open .sta file in external text editor."""
        # TODO Run as thread
        if os.path.isfile(s.path_editor):
            if os.path.isfile(self.sta):
                command = [s.path_editor, self.sta]
                subprocess.Popen(command)
            else:
                logging.error('File not found:\n' \
                    + self.sta \
                    + '\nSubmit analysis first.')
        else:
            logging.error('Wrong path to text editor:\n' \
                + s.path_editor \
                + '\nConfigure it in File->Settings.')

    """Menu Job."""

    def write_input(self, lines, file_name=''):
        """Write the whole model inp_code into the output .inp-file.
        Is called from menu 'Job -> Write input'.
        Reinitialize job because of possible file_name change.
        """
        if not len(file_name):
            file_name = QtWidgets.QFileDialog.getSaveFileName(None, \
                'Write INP file', self.dir, \
                'Input files (*.inp)')[0]
        if len(file_name):
            if not file_name.endswith('.inp'):
                file_name += '.inp'
            with open(file_name, 'w') as f:
                f.writelines(lines)
            logging.info('Input written to\n' + file_name)
            self.generate(file_name)

            # Reopen CGX
            has_nodes = len(m.Mesh.nodes)
            from gui import cgx
            cgx.open_inp(self.inp, has_nodes)

    def open_inp(self):
        """Open INP file in external text editor."""
        if os.path.isfile(s.path_editor):
            if os.path.isfile(self.inp):
                command = [s.path_editor, self.inp]
                p = subprocess.Popen(command)
                return p
            else:
                msg = 'File not found:\n' \
                    + self.inp \
                    + '\nWrite input first.'
                logging.error(msg)
                return None
        else:
            msg = 'Wrong path to text editor:\n' \
                + s.path_editor \
                + '\nConfigure it in File->Settings.'
            logging.error(msg)
            return None

    def open_subroutine(self):
        """Dialog window to filter fortran subroutines."""
        if os.path.isfile(s.path_editor):
            file_name = QtWidgets.QFileDialog.getOpenFileName(None,
                'Open a subroutine', p.ccx, 'FORTRAN (*.f)')[0]
            if file_name:
                command = [s.path_editor, file_name]
                subprocess.Popen(command)
        else:
            logging.error('Wrong path to text editor:\n' \
                + s.path_editor \
                + '\nConfigure it in File->Settings.')

    def rebuild_ccx(self):
        """Recompile CalculiX sources with updated subroutines."""
        global thread_counter

        # Windows
        if os.name == 'nt':

            # Path to ccx sources
            ccx = path2cygwin(p.ccx)

            # Open bash and send command to build CalculiX
            cmd1 = 'C:\\cygwin64\\bin\\bash.exe --login'
            send1 = '/bin/make -f ../../config/ccx_Makefile_MT_windows -C {}'.format(ccx)

            # Move binary
            cmd2 = 'C:\\cygwin64\\bin\\mv.exe -T ' \
                    + ccx + '/ccx ' \
                    + p.bin + '/ccx'

        # Linux
        else:

            # Build CalculiX
            cmd1 = ['make', '-f', '../../config/ccx_Makefile_MT_linux', '-C', p.ccx]
            send1 = ''

            # Move binary
            cmd2 = ['mv', '-T', p.ccx + '/ccx',
                    p.bin + '/ccx']

        # Build CalculiX
        if type(cmd1) == str:
            logging.info(cmd1 + ' ' + send1)
        else:
            logging.info(' '.join(cmd1) + ' ' + send1)
        thread_counter += 1
        t_name = 'thread_{}_rebuild_ccx'.format(thread_counter)
        t = threading.Thread(target=self.run,
            args=(cmd1, send1), name=t_name, daemon=True)
        t.start()

    def submit(self):
        """Submit INP to CalculiX. Calculation starts in self.run method,
        which is called via thread to avoid GUI freeze.
        """
        global thread_counter

        if not os.path.isfile(p.path_ccx):
            logging.error('CCX not found:\n' \
                + p.path_ccx)
            return

        if os.path.isfile(self.inp):
            os.environ['OMP_NUM_THREADS'] = str(os.cpu_count()) # enable multithreading
            cmd = [p.path_ccx, '-i', self.path]
            logging.info(' '.join(cmd))

            thread_counter += 1
            t_name = 'thread_{}_submit_ccx'.format(thread_counter)
            t = threading.Thread(target=self.run,
                args=(cmd, '', True), name=t_name, daemon=True)
            t.start()
        else:
            logging.error('File not found:\n' \
                + self.inp \
                + '\nWrite input first.')

    def view_log(self):
        """Open log file in external text editor."""
        if os.path.isfile(s.path_editor):
            if os.path.isfile(self.log):
                command = [s.path_editor, self.log]
                p = subprocess.Popen(command)
                return p
            else:
                msg = 'File not found:\n' \
                    + self.log \
                    + '\nSubmit analysis first.'
                logging.error(msg)
                return None
        else:
            msg = 'Wrong path to text editor:\n' \
                + s.path_editor \
                + '\nConfigure it in File->Settings.'
            logging.error(msg)
            return None

    def export_vtu(self):
        """Convert FRD to VTU."""
        from ccx2paraview import ccx2paraview
        if os.path.isfile(self.frd):
            ccx2paraview.Converter(self.frd, ['vtu']).run()
            os.system('mkdir paraview')
            os.system('mv *.vtu paraview')
            os.system('mv *.pvd paraview')
        else:
            logging.error('File not found:\n' \
                + self.frd \
                + '\nSubmit analysis first.')

    def open_paraview(self):
        """Open VTU in ParaView."""
        if os.path.isfile(s.path_paraview):

            # Count result VTU files
            file_list = []
            for f in os.listdir(self.dir + '/paraview'):
                f = os.path.basename(f)
                if f.lower() == self.name[:-4] + '.vtu':
                    file_list = [f]
                    break
                if f.lower().endswith('.vtu') and f.startswith(self.name[:-4]):
                    file_list.append(f)
            if len(file_list) > 1:
                #vtu_path = self.path + '...vtu'
                f_name, f_extension = os.path.splitext(f)
                f_name, f_extension = os.path.splitext(f_name)
                vtu_path = self.dir + '/paraview/' + f_name + '.pvd'
            elif len(file_list) == 1:
                vtu_path = self.path + '.vtu'
            else:
                logging.error('VTU file not found.\nExport VTU results first.')
                return

            command = [s.path_paraview, '--data=' + vtu_path]
            logging.info(' '.join(command))
            subprocess.Popen(command)
        else:
            logging.error('Wrong path to ParaView:\n' \
                + s.path_paraview \
                + '\nConfigure it in File->Settings.')

    def run(self, cmd, send='', read_output=True):
        """Run a single command, wait for its completion and log stdout.
        Doesn't block GUI, because is called in a separate thread.
        """
        while True:
            """Wait for previous thread to finish.
            t_names is a list of currently running threads names."""
            t_name = threading.current_thread().name
            t_names = sorted([t.name for t in threading.enumerate() \
                if '_ccx' in t.name])
            if t_name != t_names[0]:
                time.sleep(1)
            else:
                break

        # Run command
        process = subprocess.Popen(cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        if len(send):
            process.stdin.write(bytes(send, 'utf8'))
            process.stdin.close()

        # Start stdout reading and logging thread
        args = [process.stdout, 'read_stdout', read_output]
        from gui import stdout
        stdout.start_reader(*args)

        while process.poll() is None:
            time.sleep(1)


# Create job object with file logging handler
j = Job()


def path2cygwin(path):
    """Converts Windows path to Cygwin path."""
    return '/cygdrive/' + \
            path[0].lower() + \
            path[2:].replace('\\', '/')
