© Ihor Mirzov, July 2019.

Distributed under GNU General Public License, version 2.

<br/><br/>



# CalculiX CAE

GUI/pre-processor for [CalculiX CrunchiX](http://dhondt.de/). Very simple, free and open source. Program is based on CalculiX 2.15 keywords hierarchy. Written in Python3, utilizes PyQt5 and VTK.

It is implied that you have already created geometry and generated mesh in some other software like [Salome-platform](https://www.salome-platform.org/). CacluliX CAE is designed to guide you through the keywords creation sequence and is aimed to help you reach correct input file with no mistakes.

<br/><br/>



# Features

- Calculix keywords hierarchy with all attributes is maintaned in simple editable text file - [ccx_dom.inp](ccx_dom.inp). Keyword's edit dialog is based on it.

- Import/export of .inp-files - one format for all needs. Program parses .inp-file and generates model on the fly, so separate format for the model is not needed. Final model is saved also with .inp format ready to be calculated by ccx.

- Solid mesh parser - [ccx_mesh.py](ccx_mesh.py) - tested on the all example CacluliX models - with ability to highlight surfaces, node and element sets after mesh import. See [test.log](test.log).

- Keyword's edit dialog shows apropriate chapter of the official HTML documentation. Edit dialog is generated on the fly reading keyword's attributes from [ccx_dom.inp](ccx_dom.inp).

- Nice icons for all keywords/tree objects.

<br/><br/>



# Screenshots

Main window with imported mesh and highlighted node set:
![Main window](img_social.png "Main window")

"New keyword" dialog:
![Create keyword dialog](img_dialog.png "Create keyword dialog")

<br/><br/>



# Prerequisites

Suppose you already have [Python 3](https://www.python.org/downloads/). So first of all install VTK, PyQt5 and PyQtWebEngine. Then run [ccx_cae.py](ccx_cae.py):

    pip3 install vtk PyQt5 PyQtWebEngine
    python3 ccx_cae.py

Here are [some examples](https://lorensen.github.io/VTKExamples/site/Python/) for getting started with VTK in Python.

To edit GUI-forms [ccx_dialog.ui](ccx_dialog.ui) and [ccx_cae.ui](ccx_cae.ui) you'll need QT designer:

    sudo apt install qttools5-dev-tools

Also you'll need [pyinstaller](https://www.pyinstaller.org/) to build release binaries:

    pip3 install pyinstaller

<br/><br/>



# How to use

Despite of your operating system it is recommended to use the source code. Simply run [ccx_cae.py](ccx_cae.py) from terminal:

    python3 ccx_cae.py

Otherwise you may download binaries (for now only Ubuntu is supported), extract archive and run binary *ccx_cae*:

    ./ccx_cae

<br/><br/>



# Downloads

Both sources and binaries could be found on [the release page](https://github.com/imirzov/ccx_cae/releases).

<br/><br/>



# TODO

**Improvements and bugfixes:**

+ Tested VTK camera for all example models.

- Windows binary

- Better ViewSurfaceWithEdges - implement it manually without EdgeVisibilityOn().

- Better surface highlight: remove original face for pure color.

- Get all mutually exclusive keyword arguments from the manual.

- During keyword's edit parse arguments and pass them to Dialog. Reparse *NODE, *ELEMENT, *NSET or *ELSET after edit.

- Do not remove comments. Include them to apropriate tree item's INP_code.

**New features:**

- Visualize BC's, constraints and loads in VTK.

- Python code (*PYTHON keyword) in INP_code for step repetition and other kind of model generation.

- Save display options.

- Import certain keyword with descendants from .inp-file.

- Import mesh from FRD, [VTK](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadLegacyUnstructuredGrid/), [VTU](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadUnstructuredGrid/).

<!--
    DistanceBetweenPoints:  
    https://lorensen.github.io/VTKExamples/site/Python/SimpleOperations/DistanceBetweenPoints/

    Screenshot:  
    https://lorensen.github.io/VTKExamples/site/Python/Utilities/Screenshot/

    Text actor for displaying model info:  
    https://lorensen.github.io/VTKExamples/site/Python/GeometricObjects/TextActor/
-->