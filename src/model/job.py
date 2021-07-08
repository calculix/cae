#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Job submition and convertion. Run a detached process and
send messages to queue. This scheme doesn't freeze Window
while analysis is running or files are converting.

Job creation indicates start of the new session.
Old log file is erased and a new one is created.
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
import ccx2paraview
import unv2ccx

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
import path
import settings
import gui
import gui.window
import log
import tests


def copy_checks_log_contents_to(job_logfile):
    """Copy logs from startup checks into the job logfile.
    TODO Remove checks.log after.
    """
    checks_log = os.path.join(path.p.src, 'checks.log')
    lines = []
    with open(checks_log, 'r') as f:
        lines = f.readlines()
    with open(job_logfile, 'w') as f:
        f.writelines(lines)
        f.write('\nAPPLICATION START\n')
    # os.remove(checks_log)


class Job:

    def __init__(self, f, m, file_name=settings.s.start_model):
        """Create job object.
        Is called twice on startup:
        first call - from cae.py,
        second - from importer.py.
        TODO Fix it
        """
        self.f = f
        self.m = m
        self.dir = os.path.dirname(os.path.abspath(file_name)) # working directory
        self.name = os.path.basename(file_name) # INP file name
        self.inp = os.path.abspath(file_name) # full path to INP file with extension
        self.path = self.inp[:-4] # full path to INP without extension
        self.frd = self.path + '.frd' # full path to job results file
        self.log = self.path + '.log' # full path to job log file

        # Remove old log_file
        if os.path.exists(self.log):
            os.remove(self.log)

        # Handler to write the job log file
        log.remove_file_handler()
        copy_checks_log_contents_to(self.log)
        log.add_file_handler(self.log)

        if sys.argv[0].endswith('.py'):
            logging.debug('Running from sources.')
        else:
            logging.debug('Running from binaries.')
        logging.info('Application home directory is:\n'\
            + path.p.app_home_dir)
        logging.info('Work directory is:\n' + self.dir)

    def convert_unv(self):
        """Convert UNV to INP."""
        # converter_path = os.path.join(path.p.bin, 'unv2ccx' + path.p.extension)
        # cmd = [converter_path, self.path + '.unv']
        # logging.info(' '.join(cmd))
        # self.run(cmd)
        unv2ccx.Converter(self.path + '.unv').run()

    def write_input(self, lines):
        """Write the whole model inp_code into the output .inp-file.
        Is called from menu 'Job -> Write input'.
        Reinitialize job because of possible file_name change.
        """
        file_name = QtWidgets.QFileDialog.getSaveFileName(None, \
            'Write INP file', self.dir, \
            'Input files (*.inp)')[0]
        if file_name:
            with open(file_name, 'w') as f:
                f.writelines(lines)
            logging.info('Input written to\n' + file_name)
            self.__init__(self.f,\
                self.m, file_name[:-4] + '.inp')

            # Reopen CGX
            has_nodes = len(self.m.Mesh.nodes)
            gui.cgx.open_inp(self.f, self.inp, has_nodes)

    def edit_inp(self):
        """Open INP file in external text editor."""
        if os.path.isfile(settings.s.path_editor):
            if os.path.isfile(self.inp):
                command = [settings.s.path_editor, self.inp]
                subprocess.Popen(command)
            else:
                logging.error('File not found:\n' \
                    + self.inp \
                    + '\nWrite input first.')
        else:
            logging.error('Wrong path to text editor:\n' \
                + settings.s.path_editor \
                + '\nConfigure it in File->Settings.')

    def open_subroutine(self):
        """Dialog window to filter fortran subroutines."""
        if os.path.isfile(settings.s.path_editor):
            file_name = QtWidgets.QFileDialog.getOpenFileName(None,
                'Open a subroutine', path.p.ccx, 'FORTRAN (*.f)')[0]
            if file_name:
                command = [settings.s.path_editor, file_name]
                subprocess.Popen(command)
        else:
            logging.error('Wrong path to text editor:\n' \
                + settings.s.path_editor \
                + '\nConfigure it in File->Settings.')

    def rebuild_ccx(self):
        """Recompile CalculiX sources with updated subroutines."""

        # Windows
        if os.name == 'nt':

            # Path to ccx sources
            ccx = path2cygwin(path.p.ccx)

            # Open bash and send command to build CalculiX
            cmd1 = 'C:\\cygwin64\\bin\\bash.exe --login'
            send1 = '/bin/make -f Makefile_MT -C {}'.format(ccx)

            # Move binary
            cmd2 = 'C:\\cygwin64\\bin\\mv.exe -T ' \
                    + ccx + '/ccx ' \
                    + path.p.bin + '/ccx'

        # Linux
        else:

            # Build CalculiX
            cmd1 = ['make', '-f', 'Makefile_MT', '-C', path.p.ccx]
            send1 = ''

            # Move binary
            cmd2 = ['mv', '-T', path.p.ccx + '/ccx',
                    path.p.bin + '/ccx']

        # Build CalculiX
        if type(cmd1) == str:
            logging.info(cmd1 + ' ' + send1)
        else:
            logging.info(' '.join(cmd1) + ' ' + send1)
        t_name = 'thread_{}_rebuild_ccx'\
            .format(threading.active_count())
        t = threading.Thread(target=self.run,
            args=(cmd1, send1), name=t_name, daemon=True)
        t.start()

        # Move binary
        t_name = 'thread_{}_rebuild_ccx'\
            .format(threading.active_count())
        t = threading.Thread(target=self.run,
            args=(cmd2, '', False), name=t_name, daemon=True)
        t.start()

    def submit(self):
        """Submit INP to CalculiX."""
        if not os.path.isfile(path.p.path_ccx):
            logging.error('CCX not found:\n' \
                + path.p.path_ccx)
            return

        if os.path.isfile(self.inp):
            os.environ['OMP_NUM_THREADS'] = str(os.cpu_count()) # enable multithreading
            cmd = [path.p.path_ccx, '-i', self.path]
            logging.info(' '.join(cmd))

            t_name = 'thread_{}_submit_{}'\
                .format(threading.active_count(), int(time.time()))
            t = threading.Thread(target=self.run,
                args=(cmd, '', False), name=t_name, daemon=True)
            t.start()
        else:
            logging.error('File not found:\n' \
                + self.inp \
                + '\nWrite input first.')

    def view_log(self):
        """Open log file in external text editor."""
        if os.path.isfile(settings.s.path_editor):
            if os.path.isfile(self.log):
                command = [settings.s.path_editor, self.log]
                subprocess.Popen(command)
            else:
                logging.error('File not found:\n' \
                    + self.log \
                    + '\nSubmit analysis first.')
        else:
            logging.error('Wrong path to text editor:\n' \
                + settings.s.path_editor \
                + '\nConfigure it in File->Settings.')

    def export_vtu(self):
        """Convert FRD to VTU."""
        if os.path.isfile(self.frd):
            ccx2paraview.Converter(self.frd, ['vtu']).run()
        else:
            logging.error('File not found:\n' \
                + self.frd \
                + '\nSubmit analysis first.')

    def open_paraview(self):
        """Open VTU in ParaView."""
        if os.path.isfile(settings.s.path_paraview):

            # Count result VTU files
            file_list = []
            for f in os.listdir(self.dir):
                f = os.path.basename(f)
                if f.lower() == self.name[:-4] + '.vtu':
                    file_list = [f]
                    break
                if f.lower().endswith('.vtu') and f.startswith(self.name[:-4]):
                    file_list.append(f)
            if len(file_list) > 1:
                vtu_path = self.path + '...vtu'
            elif len(file_list) == 1:
                vtu_path = self.path + '.vtu'
            else:
                logging.error('VTU file not found.\nExport VTU results first.')
                return

            command = [settings.s.path_paraview, '--data=' + vtu_path]
            logging.info(' '.join(command))
            subprocess.Popen(command)
        else:
            logging.error('Wrong path to ParaView:\n' \
                + settings.s.path_paraview \
                + '\nConfigure it in File->Settings.')

    def run(self, cmd, send='', read_output=True):
        """Run a single command, wait for its completion and log stdout.
        Doesn't block GUI if called via thread.
        """
        while True:
            # Wait for previous thread to finish
            t_name = threading.current_thread().name
            t_names = sorted([t.name for t in threading.enumerate() \
                if '_rebuild_ccx' in t.name])
            if not len(t_names) or t_name == t_names[0]:
                break
            else:
                time.sleep(1)

        # Run command
        os.chdir(self.dir)
        process = subprocess.Popen(cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        if len(send):
            process.stdin.write(bytes(send, 'utf8'))
            process.stdin.close()
        os.chdir(path.p.app_home_dir)

        # Start stdout reading and logging thread
        sr = log.StdoutReader(process.stdout, 'read_stdout', read_output)
        self.f.stdout_readers.append(sr)
        sr.start()

        # Do not finish thread until the process end up
        while process.poll() is None:
            time.sleep(1)


def path2cygwin(path):
    """Converts Windows path to Cygwin path."""
    return '/cygdrive/' + \
            path[0].lower() + \
            path[2:].replace('\\', '/')


@tests.test_wrapper()
def test():
    j = Job(None, None)
    j.view_log()


if __name__ == '__main__':
    test() # run test
