© Ihor Mirzov, 2019-2021  
Distributed under GNU General Public License v3.0

[![Github All Releases](https://img.shields.io/github/downloads/calculix/cae/total.svg)](https://github.com/calculix/cae/releases)

<br/><br/>



---

[Downloads](https://github.com/calculix/cae/releases) |
[Features](#features) |
[Screenshots](#screenshots) |
[How to use](#how-to-use) |
[Your help](#your-help) |
[For developers](#for-developers) |
[What's new](#whats-new) |
[TODO](#todo)

---

<br/><br/>



# CalculiX Adanced Environment (CAE)

CAE is a software package mainly consisting of CalculiX [GraphiX](http://calculix.de/), [CrunchiX](http://dhondt.de/) and keyword editor. The last one - is a very simple, free and open source GUI/pre-processor for CalculiX. Program is based on the keywords hierarchy, is designed to guide you through the keywords creation process and is aimed to help you reach the correct input file with no mistakes. Keyword editor is written in Python3 and utilizes PyQt5.

It is implied that you have already created geometry and generated mesh in some other software like [FreeCAD](https://www.freecadweb.org/) or [Salome-platform](https://www.salome-platform.org/).

<br/><br/>



# Features

- CAE comes with the latest CGX and CCX compiled for both Linux and Windows. Analysis [runs on all cores](https://github.com/calculix/cae/issues/22) automatically.

- INP format for all needs: program parses .inp-file and generates model on the fly, so separate format for the model is not needed. Final model is saved also with .inp format ready to be calculated with CCX. Robust INP importer algorithm is tested on over 20 000 INP files, including Abaqus models ([see log](src/importer.log)).


- Official [HTML documentation](doc) is natively integrated into the keyword edit dialogs (Help button).

- Calculix keywords with all attributes is maintaned in [editable XML file](config/kw_list.xml).

- [Solid mesh parser](src/model/parsers/mesh.py) supports includes in the input file. Tested on the all official CalculiX examples. See [mesh.log](src/model/parsers/mesh.log).

- Application's global settings could be set up in the File->Settings menu. Settings are maintained in editable Python file - it is automatically overwritten during the workflow.

- [Optimal job management](src/model/job.py):

    - if you use subroutines, CalculiX sources could be automatically recompiled from GUI;
    - run analysis directly from GUI;
    - open results in GraphiX or convert to VTU format and open them in [Paraview](https://www.paraview.org).

- CAE uses converters:

    - [ccx2paraview](https://github.com/calculix/ccx2paraview) - CalculiX to Paraview converter (frd to vtk/vtu)
    - [unv2ccx](https://github.com/calculix/unv2ccx) - Salome universal to CalculiX converter (unv to inp)

- Cute modern design with [nice icons](img).


<br/><br/>



# Screenshots

Keyword editor and CGX windows with imported mesh and painted surfaces:
![img1](img/img_surfaces.png "Painted surfaces")

"New keyword" dialog with corresponding chapter from HTML manual:
![img2](img/img_dialog.png "Keyword dialog")

Keyword editor and CGX windows with calculated model:
![img3](img/img_social.png "Calculation results")

New CGX colormaps - viridis, inferno, turbo:
![img4](img/img_cmaps.png "New CGX colormaps")

Calculation result exported to Paraview:
![img5](img/img_paraview.png "Results in Paraview")

<br/><br/>



# How to use

Get [Python3](https://www.python.org/downloads/) (>3.8) and install it.

[Download latest release](https://github.com/calculix/cae/releases) (source code) and unpack it.

Open terminal in a directory where you've unpacked CAE. Allow all files to be executed (Linux only):

    sudo chmod -R 777 ./*

Install dependencies with command:

    pip3 install -r requirements.txt

It's OK if Xlib installation will fail in Windows. PyQt5 in Ubuntu also could be installed in this way:
    
    sudo apt install python3-pyqt5 python3-pyqt5.qtwebengine

Additional fonts for CGX in Ubuntu:

    sudo apt get install fonts-terminus-otb xfonts-terminus xfonts-terminus-oblique


That's all, enjoy!

Run the software with command (or double click those files):

    in Linux:       ./cae.sh
    in Windows:     cae.bat

You can edit default startup model name in File->Settings or leave it empty or even pass it as an argument to open on startup:

    in Linux:       ./cae.sh -inp yourmodel.inp
                    ./cae.sh -inp yourmodel.unv
    in Windows:     cae.bat -inp yourmodel.inp
                    cae.bat -inp yourmodel.unv

Another way to run CAE - is directly from 'src' directory:

    python3 ./src/cae.py
    python3 ./src/cae.py -inp yourmodel.inp
    python3 ./src/cae.py -inp yourmodel.unv

The intended workflow is:

- create geometry and mesh in [FreeCAD](https://www.freecadweb.org/) or [Salome-platform](https://www.salome-platform.org/),
- save mesh as UNV or export it to INP with [Salome to CalculiX mesh exporter](https://github.com/psicofil/SalomeToCalculix),
- import INP or UNV mesh into CAE and continue creating model,
- if needed, edit Fortran subroutines and rebuild CCX (menu Job),
- submit job from CAE,
- view job result in GraphiX or export it to post-process in Paraview.

In Windows to work with subroutines and to recompile CalculiX sources from CAE you'll need [cygwin](https://www.cygwin.com/). Install it to 'C:\\cygwin64' with:
- gcc-g++ 11.3.0-1
- gcc-fortran 11.3.0-1
- make 4.3-1
- perl 5.32.1-2
- libglut3 3.2.1-1

**Attention!** Do not move CCX or CGX from 'bin' directory!

<br/><br/>



# Your help

Please, you may:

- Star this project.
- Simply use this software and ask questions.
- Share your models and screenshots.
- Report problems by [posting issues](https://github.com/calculix/cae/issues). Please, attach ./cae.log (or other logs) to the problem report.
- Do something from the [TODO-list](#TODO).

<br/><br/>



# For developers

You may also need libraries:

    qttools5-dev-tools
    pycallgraph

<br/><br/>



# What's new

What's new in future v0.9.0:

+ CalculiX 2.20.
+ Significantly improved INP importer algorithm. Now parser supports keyword line continuation. Tested on over 20 000 INP files, including Abaqus models.
+ New checker module. Checks are called on the application start. OS name, Python version, CAE version and default web browser are logged. Requirements are installed automatically via pip. Statrup checks could be switched off from the settings.
+ Refactored logging system. Now it is a new dedicated module.
+ Significantly improved window connectivity (master/slave).
+ Improved robustness - now almost every module has a test. Test system uses 'unittest' standard library.
+ Simplified code to omit redundant arguments passing between objects. Class instances are created directly in modules and imported to any other module.
+ Using external python packages *ccx2paraview* and *unv2ccx* as required dependacies.
+ Added buttons to select paths in Settings dialog.
+ cae.bat and cae.sh - run source code
+ Unittest: each keyword from kw_list.xml has corresponding html file in doc.
+ Keywords hierarchy and keywords list with arguments now are driven by two separate xml-configs.
+ Added validation to the KeywordDialog widgets.
+ Bugfix: import menu did nothing in empty model.
+ Keyword Dialog: support number input (int/float).

<br/><br/>



# TODO

- Table widget: redundant comma in the end of line.
- Text widget has to be resizable.
- Argument with "use": link to create a new IMPLEMENTATION.
- Insert Combo widget into Table cell.
- Draw a separator (black line) in GUI for group with newline=1.
- Argument form='Combo' with sub-arguments: regenerate GUI on Combo index change. Current form='VOr' is not user friendly.
- Unittests for KeywordDialog: compare generated INP code with the one which should be.
- Code, docs and KOM for v 2.20.

https://www.x.org/releases/X11R7.5/doc/man/man3/XSetFont.3.html

<pre>
INFO: /media/ihor/WORK/Programming/Calculix/cae/bin/cgx -o "/media/ihor/WORK/Programming/MIPT/8_term_FEM/model.frd"
DEBUG: Starting new window: 0x04200001, 50087, CalculiX GraphiX
X protocol error:
<class 'Xlib.error.BadWindow'>: code = 3, resource_id = <class 'Xlib.xobject.resource.Resource'>(0x04200001), sequence_number = 9, major_opcode = 42, minor_opcode = 0
X protocol error:
<class 'Xlib.error.BadWindow'>: code = 3, resource_id = <class 'Xlib.xobject.resource.Resource'>(0x04200001), sequence_number = 15, major_opcode = 12, minor_opcode = 0
WARNING: No slave window.
WARNING: No slave window.
</pre>

<pre>
Traceback (most recent call last):
  File "/media/ihor/WORK/Programming/Calculix/cae/src/gui/tree.py", line 234, in clicked
    _set.append(int(n))
ValueError: invalid literal for int() with base 10: 'Nx0'

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/media/ihor/WORK/Programming/Calculix/cae/src/gui/tree.py", line 237, in clicked
    _set.extend([n.num for n in m.Mesh.nsets[n].items])
KeyError: 'Nx0'
</pre>


Examples and training materials:
- Prool's INP-templates and snippets.
- Interface for [Martin Kraska examples](https://github.com/mkraska/CalculiX-Examples).
- Menu Help: wiki, YouTube videos for beginners. 

Python API:
- Invent some kind of Python API. Possibly use [pycalculix](https://github.com/spacether/pycalculix), [pyccx](https://github.com/drlukeparry/pyccx) and [pycgx](https://github.com/zagl/pycgx). Here are [pycgx examples](https://github.com/zagl/ccx-user).
- Pay attention to [pyvista](https://www.pyvista.org/) and [Python tetgen interface](https://tetgen.pyvista.org/).
- [Meshpresso Mesh Converter](https://cofea.readthedocs.io/en/latest/meshpresso/index.html)
- Python code (*PYTHON keyword) in the INP_code for step repetition and other kind of model generation.
- Martin Kraska: support the param.py parametrization mechanism.

FreeCAD:
- Embeded CAE into FreeCAD.
- Connector for FreeCAD: listen to port and obtain meshed geometry.

CGX:
- Master window / toolset with list of CGX commands.
- Highlight loads and boundary conditions.

Parsers:
- Reparse mesh/model after tree.actionDeleteImplementation.
- Reparse mesh/model after new item creation. New elsets/nodesets are not highlighted in CGX.
- Parsers for loads and boundary conditions.

Importer:
- treeView: import certain keyword with descendants (tree branch) from INP-file.
- Interface for materials import: right click on *Material -> Import.
- [Meshio converter](https://github.com/nschloe/meshio).
- Import mesh from FRD, [VTK](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadLegacyUnstructuredGrid/), [VTU](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadUnstructuredGrid/), [Gmsh](http://gmsh.info/).
- Open .fbd/.fbl and forward to CGX. Then import generated model.

Other:
- Main Window: filter tree items - search item. 
- Main Window: treeView - highlight background of the collections.