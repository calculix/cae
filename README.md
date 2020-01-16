© Ihor Mirzov, January 2020  
Distributed under GNU General Public License v3.0

<br/><br/>



---

[Features](#features) |
[Screenshots](#screenshots) |
[Downloads](#downloads) |
[How to use](#how-to-use) |
[Your help](#your-help) |
[For developers](#for-developers) |
[What's new](#whats-new) |
[TODO](#todo)

---

<br/><br/>



# CalculiX CAE

GUI/pre-processor for [CalculiX CrunchiX](http://dhondt.de/). Very simple, free and open source. Program is based on CalculiX keywords hierarchy. Written in Python3, utilizes PyQt5 and VTK.

It is implied that you have already created geometry and generated mesh in some other software like [FreeCAD](https://www.freecadweb.org/) or [Salome-platform](https://www.salome-platform.org/). CacluliX CAE is designed to guide you through the keywords creation sequence and is aimed to help you reach correct input file with no mistakes.

<br/><br/>



# Features

- Calculix keywords hierarchy with all attributes is maintaned in [editable XML file](config/kom.xml).

- "New keyword" dialog shows apropriate chapter of the official [HTML documentation](doc).

- INP format for all needs: program parses .inp-file and generates model on the fly, so separate format for the model is not needed. Final model is saved also with .inp format ready to be calculated by CalculiX.

- [Solid mesh parser](src/model/parsers/mesh.py) supports includes in the input file. Tested on the all official CacluliX examples. See [tests.log](src/tests.log).

- Surfaces and sets of the imported mesh could be highlighted in the [VTK module](src/gui/vtk.py).

- Application's global settings could be set up in the File->Settings menu. Settings are maintained in editable env-file with Python syntax. The file is automatically overwritten during the workflow.

- [Optimal job management](src/model/job.py):

    - if you use subroutines, CalculiX sources could be automatically recompiled from GUI;
    - run analysis directly from GUI - you'll be notified on job completion;
    - open results in GraphiX or convert to VTU format and open them in [Paraview](https://www.paraview.org).

- Supports UNV mesh import.

- Cute modern design with [nice icons](img).


<br/><br/>



# Screenshots

Main window with imported mesh and highlighted node set:
![Main window](img/img_social.png "Main window")

"New keyword" dialog with corresponding chapter from HTML manual:
![Create keyword dialog](img/img_dialog.png "Create keyword dialog")

Main window and "new keyword" dialog in a simple view mode with hidden VTK widget and HTML help:
![Simple view mode](img/img_simple.png "Simple view mode")

Calculation result exported to Paraview:
![Results in Paraview](img/img_paraview.png "Results in Paraview")

<br/><br/>



# Downloads

Both sources and binaries could be found on [the releases page](https://github.com/imirzov/ccx_cae/releases).

<br/><br/>



# How to use

CacluliX CAE is portable software and doesn't need to be installed. Just extract archive, make binary executable and run it with double click. **To view results in GraphiX and Paraview configure pathes in File->Settings**.

The intended workflow is:

- create geometry and mesh in [Salome-platform](https://www.salome-platform.org/),
- save mesh as UNV or export it to INP with [salome2ccx GUI tool](https://github.com/imirzov/salome2ccx),
- import INP or UNV mesh to CAE and continue creating model,
- if needed, edit Fortran subroutines and rebuild ccx (Job->Rebuild CalculiX),
- submit job from CAE,
- export job result to the Paraview post-processor or view it in GraphiX.

Edit default startup model name in in File->Settings or you can pass it as an argument to open at start:

    in Linux:       ./cae.sh -inp model.inp
                    ./cae.sh -inp model.unv
    in Windows:     cae.bat -inp model.inp
                    cae.bat -inp model.unv

<br/><br/>



# Your help

Please, you may:

- Simply use this software and ask questions.
- Share your models and screenshots.
- Report problems by [posting issues](https://github.com/imirzov/ccx_cae/issues).
- Do something from the [TODO-list](#TODO).

<br/><br/>



# For developers

Suppose you already have [Python 3](https://www.python.org/downloads/). Open terminal in CalculiX CAE directory and install dependancies with command:

    pip3 install -r requirements.txt

From 'src' directory run source code with one of the commands:

    python3 cae.py
    python3 cae.py -inp model.inp
    python3 cae.py -inp model.unv

Use [make_release.py](make_release.py) to create releases:

    python3 make_release.py

CalculiX CAE uses external converters:

- [ccx2paraview](https://github.com/imirzov/ccx2paraview) - CalculiX to Paraview converter (frd to vtk/vtu)
- [unv2ccx](https://github.com/imirzov/unv2ccx) - Salome universal to CalculiX converter (unv to inp)

The source code for the CalculiX solver is taken from the [free_form_fortran project](https://github.com/imirzov/ccx_free_form_fortran) - this version has much better code folding and readability.

In Windows to work with subroutines and to recompile CalculiX sources from CAE you'll need *cygwin*. Install it to 'C:\\cygwin64' with:
- gcc-g++ 7.3.0-3
- gcc-fortran 7.3.0-3
- make 4.2.1-2
- perl 5.26.3-2

**Attention!** In Windows CCX may not work if placed outside 'bin' directory. It needs Cygwin DLLs! So, if you change path to CCX, place there Cygwin DLLs from 'bin' directory.

<br/><br/>



# What's new

- CalculiX 2.16
- Mesh object is reparsed after edit in the treeView.

<br/><br/>



# TODO

## Improvements and bugfixes

- Better ViewSurfaceWithEdges - implement it manually without EdgeVisibilityOn().

- Better surface highlight: remove original face for pure color.

- During keyword's edit parse arguments and pass them to Dialog.

## New features

- It would be nice if different ELSETS were assigned different colors for visualization.

- treeView: show implementations only.

- Visualize BC's, constraints and loads in VTK.

- Python code (*PYTHON keyword) in INP_code for step repetition and other kind of model generation.

- Import certain keyword with descendants (tree branch) from INP-file.

- Interface for materials import. Enrich Materials library.

- Import mesh from FRD, [VTK](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadLegacyUnstructuredGrid/), [VTU](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadUnstructuredGrid/), [Gmsh](http://gmsh.info/), Abaqus INP and ODB.

- Prool's INP-templates and snippets.

- Include Martin Kraska examples.

- Menu Help:
    - PDF manual / YouTube videos for beginners.
