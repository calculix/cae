© Ihor Mirzov, July 2020  
Distributed under GNU General Public License v3.0

<br/><br/>



---

[Features](#features) |
[Screenshots](#screenshots) |
[Downloads](https://github.com/calculix/cae/releases) |
[How to use](#how-to-use) |
[Your help](#your-help) |
[For developers](#for-developers) |
[What's new](#whats-new) |
[TODO](#todo)

---

<br/><br/>



# CalculiX Adanced Environment (CAE)

CAE is a software package mainly consisting of CalculiX [GraphiX](http://calculix.de/), [CrunchiX](http://dhondt.de/) and keyword editor. The last one - is a very simple, free and open source GUI/pre-processor for CalculiX. Program is based on the keywords hierarchy, is designed to guide you through the keywords creation process and is aimed to help you reach correct input file with no mistakes. Keyword editor is written in Python3 and utilizes PyQt5.

It is implied that you have already created geometry and generated mesh in some other software like [FreeCAD](https://www.freecadweb.org/) or [Salome-platform](https://www.salome-platform.org/).

<br/><br/>



# Features

- CAE comes with the latest CGX and CCX compiled for both Linux and Windows. Windows version of CGX is taken from [CalculiX Launcher v3.4](http://www.calculixforwin.com/) 

- INP format for all needs: program parses .inp-file and generates model on the fly, so separate format for the model is not needed. Final model is saved also with .inp format ready to be calculated with CCX.

- Official [HTML documentation](doc) is natively integrated into the keyword edit dialogs (Help button).

- Calculix keywords hierarchy with all attributes is maintaned in [editable XML file](config/kom.xml).

- [Solid mesh parser](src/model/parsers/mesh.py) supports includes in the input file. Tested on the all official CalculiX examples. See [mesh.log](src/model/parsers/mesh.log).

- Application's global settings could be set up in the File->Settings menu. Settings are maintained in editable Python file - it is automatically overwritten during the workflow.

- [Optimal job management](src/model/job.py):

    - if you use subroutines, CalculiX sources could be automatically recompiled from GUI;
    - run analysis directly from GUI;
    - open results in GraphiX or convert to VTU format and open them in [Paraview](https://www.paraview.org).

- Supports UNV mesh import.

- Cute modern design with [nice icons](img).


<br/><br/>



# Screenshots

Keyword editor and CGX windows with imported mesh and painted surfaces:
![img1](img/img_surfaces.png "Painted surfaces")

"New keyword" dialog with corresponding chapter from HTML manual:
![img2](img/img_dialog.png "Keyword dialog")

Keyword editor and CGX windows with calculated model:
![img3](img/img_social.png "Calculation results")

Calculation result exported to Paraview:
![img4](img/img_paraview.png "Results in Paraview")

<br/><br/>



# How to use

[Download latest release](https://github.com/calculix/cae/releases) (binaries or source code), unpack it and allow to be executed (give permissions).

The intended workflow is:

- create geometry and mesh in [FreeCAD](https://www.freecadweb.org/) or [Salome-platform](https://www.salome-platform.org/),
- save mesh as UNV or export it to INP with [Salome to CalculiX mesh exporter](https://github.com/psicofil/SalomeToCalculix),
- import INP or UNV mesh into CAE and continue creating model,
- if needed, edit Fortran subroutines and rebuild CCX (menu Job->Rebuild CalculiX),
- submit job from CAE,
- view job result in GraphiX or export it to the Paraview post-processor.

Run the software with command (works both for binaries and source code):

    in Linux:       ./cae.sh
    in Windows:     cae.bat

**Attention!** To run sources you'll need Python3 with PyQt5. Linux also needs Xlib. See [requirements.txt](./requirements.txt).

You can edit default startup model name in File->Settings or leave it empty or even pass it as an argument to open on startup:

    in Linux:       ./cae.sh -inp yourmodel.inp
                    ./cae.sh -inp yourmodel.unv
    in Windows:     cae.bat -inp yourmodel.inp
                    cae.bat -inp yourmodel.unv

Another way to run CAE is (from 'src' directory):

    python3 cae.py
    python3 cae.py -inp yourmodel.inp
    python3 cae.py -inp yourmodel.unv

<br/><br/>



# Your help

Please, you may:

- Simply use this software and ask questions.
- Share your models and screenshots.
- Do something from the [TODO-list](#TODO).
- Report problems by [posting issues](https://github.com/calculix/cae/issues). Please, attach ./examples/default.log to the problem report.

<br/><br/>



# For developers

Suppose you already have [Python 3](https://www.python.org/downloads/). Open terminal in CAE directory and install dependencies with command:

    pip3 install -r requirements.txt

You may also need libraries:

    qttools5-dev-tools
    pyinstaller
    pycallgraph

Use [make_release.py](make_release.py) to create releases (pyinstaller needed):

    python3 make_release.py

CAE uses converters:

- [ccx2paraview](https://github.com/calculix/ccx2paraview) - CalculiX to Paraview converter (frd to vtk/vtu)
- [unv2ccx](https://github.com/calculix/unv2ccx) - Salome universal to CalculiX converter (unv to inp)

The source code for the CalculiX solver is taken from the [free_form_fortran project](https://github.com/calculix/free_form_fortran) - this version has much better code folding and readability.

In Windows to work with subroutines and to recompile CalculiX sources from CAE you'll need [cygwin](https://www.cygwin.com/). Install it to 'C:\\cygwin64' with:
- gcc-g++ 7.3.0-3
- gcc-fortran 7.3.0-3
- make 4.2.1-2
- perl 5.26.3-2

**Attention!** Do not move CCX or CGX from 'bin' directory!

<br/><br/>



# What's new

[Materials library](https://github.com/calculix/examples/tree/master/materials)!

Do not skip comments while parsing the INP file. Keep comments on file save.

Replaced converters binaries with source code version. It reduced distro weight from 15 to 30 MB. 

Both source code and binaries could be started via cae.sh and cae.bat.

Interface and usability:
+ Notify on job completion
+ treeView: collapse/expand items with left/right buttons.
+ Run CGX without cmd in Windows
+ Do not start CGX with no mesh.
+ Allow to start CAE with empty model.
+ Paint CGX logs with #darkslategray.
+ Updated a few images, brighter logo icon.
+ Menu Help: README.pdf.
+ New menu - "CGX" - devoted to GraphiX. Moved CGX actions from Job to CGX menu.

Settings:
+ New setting "Start CGX by default".
+ Made Gedit default text editor for Linux.
+ Now a user can't choose path to CGX and CCX - only embed version should be used.

Default model - default.inp:
+ Updated surfaces.
+ Added some comments.

Bugfixes:
+ Allow app to run even without WIDs. Now CAE should always start.
+ Stop all stdout reading threads on file import - app doesn't die on logs scrolling.
+ Fixed hiding Import dialog (QFileDialog visibility).
+ Fixed Jobs. Now they run via thread and doesn't block GUI.
+ Fixed small bug in surface parser.
+ Fixed Job -> open_paraview.
+ Fixed command CGX -> Paint sets.

Logging system:
+ Stdout reading class for better logging of child processes.
+ Logging of how the app is running: from source code or as binary.
+ Improved logging for mesh and keywords parsers.
+ Switch off redundant logging in keyword dialog.
+ Flush CGX buffer on startup - it enforces logging to start without unwanted delay.
+ Log window list on error getting WID.
+ Allow to log empty lines - for some space between messages.

CGX:
+ Correctly kill CGX in Windows.
+ Filter accidental key presses on CGX.
+ Reduced methods dependencies in cgx.py.
+ Register additional colors in CGX on startup via .fbd commands.
+ Workaround for iso view in Windows (with .fbd commands).
+ Slightly more precise colors for surfaces and sets.

General code improvements:
+ Simplified stop_stdout_readers().
+ Job class without initialize() method.
+ Better window match by title.
+ Submit job via thread.
+ No psutil dependency.
+ Make KOM independent on Path and Settings.
+ Renamed KOM items to follow PEP8.
+ Correctly initialize window's parent class.
+ New save method for Settings.
+ New tests for KOM class and Importer. Updated test for mesh parser.

<br/><br/>



# TODO

- Implement keyword dialog as tab on main window.
- Settings Dialog: add buttons to Pathes to open files.

- Invent some kind of Python API for CAE. Possibly use [pycalculix](https://github.com/spacether/pycalculix) and [pyccx](https://github.com/drlukeparry/pyccx).
- Python code (*PYTHON keyword) in the INP_code for step repetition and other kind of model generation.

- treeView: import certain keyword with descendants (tree branch) from INP-file.

- Interface for materials import: right click on *Material -> Import.
- Prool's INP-templates and snippets.

- [meshio converter](https://github.com/nschloe/meshio)
- Import mesh from FRD, [VTK](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadLegacyUnstructuredGrid/), [VTU](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadUnstructuredGrid/), [Gmsh](http://gmsh.info/), Abaqus INP and ODB.


- Forward .fbd/.fbl to CGX.
- [Martin Kraska examples](https://github.com/mkraska/CalculiX-Examples)
- M.Kraska: support the param.py parametrization mechanism.

- Menu Help: wiki, YouTube videos for beginners. 

- CAE should be independent and be able to be embeded into FreeCAD.
- Connector for FreeCAD: listen to port and obtain meshed geometry.

- Change CGX colormaps via Xlib.

- Parser: reparse mesh/model after tree.actionDeleteImplementation.

- Is it possible to make release binaries without pyinstaller?

- Thread for periodic checking of CGX WID. Connect/disconnect CGX automatically in the background.

- https://mechanicalhacks.wordpress.com/2011/03/29/building-calculix-to-run-native-on-64-bit-windows/
