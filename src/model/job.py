#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, July 2020
Distributed under GNU General Public License v3.0

Job submition and convertion. Run a detached process and
send messages to queue. This scheme doesn't freeze Window
while analysis is running or files are converting. """


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
sys_path = os.path.dirname(__file__)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys.path.append(sys_path)
import gui
from utils.ccx2paraview import ccx2paraview
from utils.unv2ccx import unv2ccx

class Job:

    # Create job object
    def __init__(self, p, s, w, m, file_name=''):
        self.p = p
        self.s = s
        self.w = w
        self.m = m
        if not len(file_name):
            file_name = self.s.start_model
        self.dir = os.path.dirname(os.path.abspath(file_name)) # working directory
        self.name = os.path.basename(file_name) # INP file name
        self.inp = os.path.abspath(file_name) # full path to INP file with extension
        self.path = self.inp[:-4] # full path to INP without extension
        self.frd = self.path + '.frd' # full path to job results file
        self.log = self.path + '.log' # full path to job log file

        # Handler to write the job's log file
        gui.log.add_file_handler(self.log)

        if sys.argv[0].endswith('.py'):
            logging.debug('Running from sources.')
        else:
            logging.debug('Running from binaries.')
        logging.info('Application\'s home directory is:\n'\
            + self.p.app_home_dir)
        logging.info('Work directory is:\n' + self.dir)

    # Convert UNV to INP
    def convert_unv(self):
        # converter_path = os.path.join(self.p.bin, 'unv2ccx' + self.p.extension)
        # cmd = [converter_path, self.path + '.unv']
        # logging.info(' '.join(cmd))
        # self.run(cmd)
        unv2ccx.Converter(self.path + '.unv').run()

    # Write the whole model's inp_code
    # into the output .inp-file.
    # Is called from menu 'Job -> Write input'
    # Reinitialize job because of possible file_name change
    def write_input(self, lines):
        file_name = QtWidgets.QFileDialog.getSaveFileName(None, \
            'Write INP file', self.dir, \
            'Input files (*.inp)')[0]
        if file_name:
            with open(file_name, 'w') as f:
                f.writelines(lines)
            logging.info('Input written to:\n' + file_name)
            self.__init__(self.p, self.s, self.w,\
                self.m, file_name[:-4] + '.inp')

            # Reopen CGX
            # gui.cgx.kill(self.w)
            has_nodes = len(self.m.Mesh.nodes)
            gui.cgx.open_inp(self.w, self.inp, has_nodes)

    # Open INP file in external text editor
    def edit_inp(self):
        if os.path.isfile(self.s.path_editor):
            if os.path.isfile(self.inp):
                command = [self.s.path_editor, self.inp]
                subprocess.Popen(command)
            else:
                logging.error('File not found:\n' \
                    + self.inp \
                    + '\nWrite input first.')
        else:
            logging.error('Wrong path to text editor:\n' \
                + self.s.path_editor \
                + '\nConfigure it in File->Settings.')

    # Dialog window to filter fortran subroutines
    def open_subroutine(self):
        if os.path.isfile(self.s.path_editor):
            file_name = QtWidgets.QFileDialog.getOpenFileName(None,
                'Open a subroutine', self.p.ccx, 'FORTRAN (*.f)')[0]
            if file_name:
                command = [self.s.path_editor, file_name]
                subprocess.Popen(command)
        else:
            logging.error('Wrong path to text editor:\n' \
                + self.s.path_editor \
                + '\nConfigure it in File->Settings.')

    # Recompile CalculiX sources with updated subroutines
    def rebuild_ccx(self):

        # Windows
        if os.name == 'nt':

            # Path to ccx sources
            ccx = path2cygwin(self.p.ccx)

            # Open bash and send command to build CalculiX
            cmd1 = 'C:\\cygwin64\\bin\\bash.exe --login'
            send1 = '/bin/make -f Makefile_MT -C {}'.format(ccx)

            # Move binary
            cmd2 = 'C:\\cygwin64\\bin\\mv.exe -T ' \
                    + ccx + '/ccx_' + self.p.ccx_version + '_MT ' \
                    + self.p.bin + '/ccx'

        # Linux
        else:

            # Build CalculiX
            cmd1 = ['make', '-f', 'Makefile_MT', '-C', self.p.ccx]
            send1 = ''

            # Move binary
            cmd2 = ['mv', '-T', self.p.ccx + '/ccx_' + self.p.ccx_version + '_MT',
                    self.p.bin + '/ccx']

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

    # Submit INP to CalculiX
    def submit(self):
        if not os.path.isfile(self.p.path_ccx):
            logging.error('CCX not found:\n' \
                + self.p.path_ccx)
            return

        if os.path.isfile(self.inp):
            os.environ['OMP_NUM_THREADS'] = str(os.cpu_count()) # enable multithreading
            cmd = [self.p.path_ccx, '-i', self.path]
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

    # Open log file in external text editor
    def view_log(self):
        if os.path.isfile(self.s.path_editor):
            if os.path.isfile(self.log):
                command = [self.s.path_editor, self.log]
                subprocess.Popen(command)
            else:
                logging.error('File not found:\n' \
                    + self.log \
                    + '\nSubmit analysis first.')
        else:
            logging.error('Wrong path to text editor:\n' \
                + self.s.path_editor \
                + '\nConfigure it in File->Settings.')

    # Convert FRD to VTU
    def export_vtu(self):
        if os.path.isfile(self.frd):
            ccx2paraview.Converter(self.frd, 'vtu').run()
        else:
            logging.error('File not found:\n' \
                + self.frd \
                + '\nSubmit analysis first.')

    # Open VTU in ParaView
    def open_paraview(self):
        if os.path.isfile(self.s.path_paraview):

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

            command = [self.s.path_paraview, '--data=' + vtu_path]
            logging.info(' '.join(command))
            subprocess.Popen(command)
        else:
            logging.error('Wrong path to ParaView:\n' \
                + self.s.path_paraview \
                + '\nConfigure it in File->Settings.')

    # Run a single command, wait for its completion and log stdout
    # Doesn't block GUI if called via thread
    def run(self, cmd, send='', read_output=True):

        # Wait for previous thread to finish
        while True:
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
        os.chdir(self.p.app_home_dir)

        # Start stdout reading and logging thread
        sr = gui.log.StdoutReader(process.stdout, 'read_stdout', read_output)
        self.w.stdout_readers.append(sr)
        sr.start()

        # Do not finish thread until the process end up
        while process.poll() is None:
            time.sleep(1)


# Converts Windows path to Cygwin path
def path2cygwin(path):
    return '/cygdrive/' + \
            path[0].lower() + \
            path[2:].replace('\\', '/')


# Run test
# if __name__ == '__main__':

#     # Configure logging
#     logging.log = print
#     logging.debug = print
#     logging.info = print
#     logging.warning = print

#     s = Settings()

#     # Create job
#     j = Job(...)

#     # Rebuild CCX
#     j.rebuildCCX()

#     # Remove cached files
#     import clean
#     clean.cache()
