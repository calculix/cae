© Ihor Mirzov, 2019.  
Distributed under GNU General Public License, version 2.

<br/><br/>



# CalculiX CAE

GUI/pre-processor for [CalculiX CrunchiX](http://dhondt.de/). Very simple, free and open source. It is implied that you have already created geometry and generated mesh in some other software like [Salome-platform](https://www.salome-platform.org/). Program is based on CalculiX 2.15 keywords hierarchy. Written in Python3. Utilizes PyQt5 and VTK.

<br/><br/>



# Features

- Calculix keywords hierarchy with all attributes is maintaned in simple editable text file - ccx_dom.inp

- Import/export of .inp-files - one format for all needs. Program parses .inp-file and generates model on the fly, so separate format for the model is not needed.

- Highlight of surfaces, node and element sets after mesh import.

- Keyword's edit dialog shows apropriate chapter of the official HTML documentation.

- Nice icons for all keywords/tree objects.

<br/><br/>



# Screenshots

Main window:
![Main window](img_main.png "Main window")

Create keyword dialog:
![Create keyword dialog](img_dialog.png "Create keyword dialog")

<br/><br/>



# Downloads

First of all you can use the source code and run *ccx_cae.py* from terminal. Pay attention to [prerequisites](#prerequisites). Or download Ubuntu binaries, extract archive and run *ccx_cae*. Both sources and binaries are on [Google Drive](https://drive.google.com/drive/folders/1rrnOfIcu9gZ42I-r5-pVQMWRTp9Pp3IW).

<br/><br/>



# Prerequisites

To run ccx_cae.py or contribute first of all install VTK, PyQt5 and QT designer:

    pip3 install vtk
    pip3 install PyQt5
    pip3 install PyQtWebEngine --user
    sudo apt install qttools5-dev-tools

[Some examples](https://lorensen.github.io/VTKExamples/site/Python/) for getting started with VTK in Python. Also you'll need [pyinstaller](https://www.pyinstaller.org/) to create executable:

    pip3 install pyinstaller

<br/><br/>



# TODO

More tests for all features.

Visualize BC's, constraints and loads in VTK.

Reparse *NODE, *ELEMENT, *NSET or *ELSET after edit 

Do not remove comments. Include them to apropriate tree item's INP_code.

Get all mutually exclusive arguments from the manual.

Python code (*PYTHON keyword) in INP_code for step repetition and other kind of model generation.

During keyword's edit parse arguments and pass them to Dialog.

Save display options.



<!--
Import mesh from FRD, VTK, VTU:  
https://lorensen.github.io/VTKExamples/site/Python/IO/ReadUnstructuredGrid/  
https://lorensen.github.io/VTKExamples/site/Python/IO/ReadLegacyUnstructuredGrid/

DistanceBetweenPoints:  
https://lorensen.github.io/VTKExamples/site/Python/SimpleOperations/DistanceBetweenPoints/

Screenshot:  
https://lorensen.github.io/VTKExamples/site/Python/Utilities/Screenshot/

Text actor for displaying model info:  
https://lorensen.github.io/VTKExamples/site/Python/GeometricObjects/TextActor/
-->