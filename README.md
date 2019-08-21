© Ihor Mirzov, August 2019  
Distributed under GNU General Public License v3.0

<br/><br/>



# CalculiX CAE

GUI/pre-processor for [CalculiX CrunchiX](http://dhondt.de/). Very simple, free and open source. Program is based on CalculiX keywords hierarchy. Written in Python3, utilizes PyQt5 and VTK.

It is implied that you have already created geometry and generated mesh in some other software like [Salome-platform](https://www.salome-platform.org/). CacluliX CAE is designed to guide you through the keywords creation sequence and is aimed to help you reach correct input file with no mistakes.

<br/><br/>



# Features

- Calculix keywords hierarchy with all attributes is maintaned in simple editable text file - [ccx_dom.inp](ccx_dom.inp). Keyword's edit dialog is based on it.

- Import/export of .inp-files - one format for all needs. Program parses .inp-file and generates model on the fly, so separate format for the model is not needed. Final model is saved also with .inp format ready to be calculated by ccx.

- Solid mesh parser - [ccx_mesh.py](ccx_mesh.py) - tested on the all example CacluliX models - with ability to highlight surfaces, node and element sets after mesh import. See [tests.log](tests.log).

- "New keyword" dialog shows apropriate chapter of the official HTML documentation. The dialog is generated on the fly reading keyword's attributes from [ccx_dom.inp](ccx_dom.inp).

- Nice icons for all keywords/tree objects.

- Run job analysis directly from CAE.

- Application's global settings are maintained in text file [Settings.env](Settings.env) with Python syntax. Settings are automatically saved during the workflow.


<br/><br/>



# Screenshots

Main window with imported mesh and highlighted node set:
![Main window](img_social.png "Main window")

"New keyword" dialog with corresponding chapter from HTML manual:
![Create keyword dialog](img_dialog.png "Create keyword dialog")

<br/><br/>



# Downloads

Both sources and binaries could be found on [the releases page](./releases).

<br/><br/>



# How to use

CacluliX CAE is portable software and doesn't need to be installed. Just extract archive and run *ccx_cae* or *ccx_cae.exe*. Don't forget to make binaries executable.

To run analysis from CAE configure path to CalculiX in Settings.env.

It is recommended to use the source code - it's crossplatform.

    crossplatform: python3 ccx_cae.py
    in Linux:      ./ccx_cae
    in Windows:    ccx_cae.exe

Pass name of your INP file as an argument to open it at start:

    crossplatform: python3 ccx_cae.py -inp model.inp
    in Linux:      ./ccx_cae -inp model.inp
    in Windows:    ccx_cae.exe -inp model.inp

So, the intended workflow is:
- create geometry and mesh in [Salome-platform](https://www.salome-platform.org/)
- export Salome mesh to INP with [SalomeToCalculix GUI tool](https://github.com/psicofil/SalomeToCalculix),
- import INP mesh to CAE and continue creating model, submit job from CAE,
- convert FRD results to VTK/VTU to post-process in Paraview with the [ccx2paraview converter](https://github.com/imirzov/ccx2paraview).

<br/><br/>



# Your help

Please, you may:

- Simply use this software and ask questions.
- Share your models and screenshots.
- Report problems by [posting issues](./issues).
- Do something from the [TODO-list](#TODO).
- Follow discussion in the [Yahoo CalculiX Group](https://groups.yahoo.com/neo/groups/CALCULIX/conversations/topics/15616)

<br/><br/>



# For developers

Suppose you already have [Python 3](https://www.python.org/downloads/). Install dependancies with command:

    pip3 install -r requirements.txt

Use [pyinstaller](https://www.pyinstaller.org/) or [release_binaries.py](release_binaries.py) to create releases:

    pyinstaller ccx_cae.py
    python3 release_binaries.py

<br/><br/>



# TODO

**Done for the next release**

- Import mesh from [UNV](https://github.com/imirzov/unv2ccx).

- Notification on job completion.

- Convert calculation results: FRD -> VTU.

- Jobs run without GUI blocking.

**Improvements and bugfixes**

- Better ViewSurfaceWithEdges - implement it manually without EdgeVisibilityOn().

- Better surface highlight: remove original face for pure color.

- During keyword's edit parse arguments and pass them to Dialog. Reparse mesh objects after edit.  

**New features**

- PDF manual or YouTube videos for beginners.

- Visualize BC's, constraints and loads in VTK.

- Python code (*PYTHON keyword) in INP_code for step repetition and other kind of model generation.

- Import certain keyword with descendants (tree branch) from INP-file.

- Interface for materials import. Enrich Materials library.

- Import mesh from FRD, [VTK](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadLegacyUnstructuredGrid/), [VTU](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadUnstructuredGrid/), [Gmsh](http://gmsh.info/).

- Prool's INP-templates and snippets.

- Include Martin Kraska examples.

- Menu Help: VTK Keyboard Shortcuts and Controls.
