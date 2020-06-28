#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, May 2020
Distributed under GNU General Public License v3.0

Job submition and convertion. Run a detached process and
send messages to queue. This scheme doesn't freeze Window
while analysis is running or files are converting. """


# Standard modules
import os
import logging
import subprocess

# External modules
from PyQt5 import QtWidgets

# My modules
import gui


class Job:

    # Create job object
    def __init__(self, p, s, w):
        self.p = p
        self.s = s
        self.w = w
        self.dir = p.examples
        msg = 'Application\'s home directory is:\n' + p.app_home_dir
        logging.info(msg)

    # Rename job
    def initialize(self, file_name):
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
        logging.info('Work directory is:\n' + self.dir)

    # Convert UNV to INP
    def convert_unv(self):
        converter_path = os.path.join(self.p.bin, 'unv2ccx' + self.p.extension)
        cmd1 = [converter_path, self.path + '.unv']
        self.run([(cmd1, ''), ])

    # Recursively write the whole model's INP_code
    # into the output .inp-file.
    # Is called from menu 'Job -> Write input'
    # TODO Doesn't write comments present in initial model
    def write_input(self, lines):
        file_name = QtWidgets.QFileDialog.getSaveFileName(None, \
            'Write INP file', self.dir, \
            'Input files (*.inp)')[0]
        if file_name:
            with open(file_name, 'w') as f:
                f.writelines(lines)
            self.initialize(file_name)
            logging.info('Input written to:\n' + file_name)

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

            # Open bash
            cmd1 = 'C:\\cygwin64\\bin\\bash.exe --login'

            # Send command to build CalculiX
            send1 = '/bin/make -f Makefile_MT -C {}'.format(ccx)

            # Copy binary
            cmd2 = 'C:\\cygwin64\\bin\\cp.exe ' \
                    + ccx + '/ccx_' + self.p.ccx_version + '_MT ' \
                    + self.s.path_ccx

            self.run([(cmd1, send1), (cmd2, '')])

        # Linux
        else:

            # Build CalculiX
            cmd1 = ['make', '-f', 'Makefile_MT', '-C', self.p.ccx]

            # Copy binary
            cmd2 = ['cp', self.p.ccx + '/ccx_' + self.p.ccx_version + '_MT',
                    self.s.path_ccx]

            self.run([(cmd1, ''), (cmd2, '')])

    # Submit INP to CalculiX
    def submit(self):
        if os.path.isfile(self.s.path_ccx):
            if os.path.isfile(self.inp):
                os.environ['OMP_NUM_THREADS'] = str(os.cpu_count()) # enable multithreading
                cmd1 = [self.s.path_ccx, '-i', self.path]
                self.run([(cmd1, ''), ])
            else:
                logging.error('File not found:\n' \
                    + self.inp \
                    + '\nWrite input first.')
        else:
            logging.error('Wrong path to CCX:\n' \
                + self.s.path_ccx \
                + '\nConfigure it in File->Settings.')

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

    # Open INP in GraphiX
    def cgx_inp(self, m):
        if os.path.isfile(self.inp):
            if not len(m.Mesh.nodes):
                gui.cgx.kill(self.w) # close old CGX
                logging.warning('Empty mesh, CGX will not start!')
                return
            self.w.run_cgx(self.s.path_cgx + ' -c ' + self.inp)
        else:
            logging.error('File not found:\n' + self.inp)

    # Open FRD in GraphiX
    def cgx_frd(self):
        if os.path.isfile(self.frd):
            self.w.run_cgx(self.s.path_cgx + ' -o ' + self.frd)
        else:
            logging.error('File not found:\n' \
                + self.frd \
                + '\nSubmit analysis first.')

    # Convert FRD to VTU
    def export_vtu(self):
        if os.path.isfile(self.frd):
            converter_path = os.path.join(self.p.bin,
                    'ccx2paraview' + self.p.extension)
            cmd1 = [converter_path, self.frd, 'vtu']
            self.run([(cmd1, ''), ])
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
                    file_list = []
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
            subprocess.Popen(command)
        else:
            logging.error('Wrong path to ParaView:\n' \
                + self.s.path_paraview \
                + '\nConfigure it in File->Settings.')

    # Run multiple commands and log stdout without blocking GUI
    def run(self, commands):
        os.chdir(self.dir)
        for cmd1, cmd2 in commands:
            logging.info(' '.join(cmd1) + ' ' + cmd2)
            process = subprocess.Popen(cmd1,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
            if len(cmd2):
                process.stdin.write(bytes(cmd2, 'utf8'))
                process.stdin.close()
            gui.log.read_output(process.stdout, self.w.cgx_process)
        os.chdir(self.p.app_home_dir)


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
