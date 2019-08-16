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

- Application's global settings are maintained in [text file](Settings.env) with Python syntax. Settings are automatically saved during the workflow.


<br/><br/>



# Screenshots

Main window with imported mesh and highlighted node set:
![Main window](img_social.png "Main window")

"New keyword" dialog with corresponding chapter from HTML manual:
![Create keyword dialog](img_dialog.png "Create keyword dialog")

<br/><br/>



# How to use

Despite of your operating system it is recommended to use the source code - simply call [ccx_cae.py](ccx_cae.py) from terminal. Otherwise you may [download binaries](https://github.com/imirzov/ccx_cae/releases). CacluliX CAE is portable software and doesn't need to be installed. Just extract archive and run *ccx_cae* or *ccx_cae.exe*. Don't forget to make all binaries executable.

    crossplatform: python3 ccx_cae.py
    in Linux:      ./ccx_cae
    in Windows:    ccx_cae.exe

Pass name of your INP file as an argument to open it at start:

    crossplatform: python3 ccx_cae.py -inp=your_model.inp
    in Linux: ./ccx_cae -inp=your_model.inp
    in Windows: ccx_cae.exe -inp=your_model.inp

<br/><br/>



# Downloads

Both sources and binaries could be found on [the release page](https://github.com/imirzov/ccx_cae/releases).

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

Suppose you already have [Python 3](https://www.python.org/downloads/). So first of all install VTK, PyQt5 and PyQtWebEngine:

    pip3 install vtk PyQt5 PyQtWebEngine

That's all if you're going only to use CalculiX CAE.

For developers to edit GUI-forms [ccx_dialog.ui](ccx_dialog.ui) and [ccx_cae.ui](ccx_cae.ui) you'll need QT designer:

    sudo apt install qttools5-dev-tools

To create release binaries install [pyinstaller](https://www.pyinstaller.org/):

    pip3 install pyinstaller

Then create release binaries with command:

    pyinstaller ccx_cae.py

<br/><br/>



# TODO

**Done for the new release**

+ Fixed *boundary line in ccx_dom.inp: added missig NEW|MOD values of parameter OP.

+ In the HTML help window line breaks are replaced by spaces, avoiding words to be glued together.

+ Fixed camera's clipping range bug due to nonconsequent point numbering in the VTK unstructured grid.

+ Fixed VTK cells numbering.


+ Now exporter writes INP code with paddings - for code folding.

+ Added debug level messages in logging.

+ Added support for mutually exclusive arguments in ccx_dom.inp. Example - *MODAL DAMPING.

+ Got all keyword arguments from the manual (still not accurately).

+ In 'New keyword' dialog propose implementation variants for corresponding keyword arguments.

+ Settings with autosave + path to CCX and CGX.

+ Run job from CAE.

+ Improved treeView:
    - non-clickable keywords are grayed out - 'active' parameter for Class item;
    - bold keywords implementations.

**TODO: improvements and bugfixes**

- Better ViewSurfaceWithEdges - implement it manually without EdgeVisibilityOn().

- Better surface highlight: remove original face for pure color.

- During keyword's edit parse arguments and pass them to Dialog. Reparse mesh objects after edit.

**TODO: new features**

- PDF manual for beginners: "I am not quite sure, what the intended workflow is. I'd expect some representation (object tree) of what has been imported and/or has been already defined. The tree in the CAE window seems to be just sort of a menu but not an actual object tree."

- Visualize BC's, constraints and loads in VTK.

- Python code (*PYTHON keyword) in INP_code for step repetition and other kind of model generation.

- Import certain keyword with descendants (tree branch) from INP-file. Interface for materials import. Enrich Materials library.

- Import mesh from FRD, [VTK](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadLegacyUnstructuredGrid/), [VTU](https://lorensen.github.io/VTKExamples/site/Python/IO/ReadUnstructuredGrid/), [Gmsh](http://gmsh.info/), UNV.

- Prool's INP-templates and snippets.

- Include 'Examples' into releases. Martin Kraska examples. HTML page with tags/groups for all INP files.

- Interactor to select nodes and elements for model construction.

- Open 'New keyword' dialog in QTabWidget, not in separate window. Ctrl+Tab to switch between tabs.

- Embed CGX into additional QTabWidget.

- Menu Help: VTK Keyboard Shortcuts and Controls.
