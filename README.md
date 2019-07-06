© Ihor Mirzov, July 2019.  
Distributed under GNU General Public License, version 2.

<br/><br/>



# CalculiX CAE

GUI/pre-processor for CalculiX CrunchiX. Very simple, lightweight, free and open source. It is implied that you have already created geometry and generated mesh in some other software like Salome-platform.

Program is based on CalculiX 2.15 keywords hierarchy. Written in Python3. Implements PyQt5 and VTK.

<br/><br/>



# Info for developers

To contribute first install vtk, PyQt5 and QT designer:

    pip3 install vtk
    pip3 install PyQt5
    pip3 install PyQtWebEngine --user
    sudo apt install qttools5-dev-tools
    <!-- sudo apt install texlive-latex-base -->

Some examples for getting started with VTK in Python:

    https://lorensen.github.io/VTKExamples/site/Python/

<br/><br/>



# TODO

Visualize all sets, surfaces, BC's and loads in VTK
Update VTK after importINP and implementation edit 

Comments code

Show PDF documentation for each keyword's dialog.

Get all mutually exclusive arguments from the manual - now they are supported.

Styles for treeView elements.

Save display options.

Python code (*PYTHON keyword) in INP_code for step repetition and other kind of model generation.

Parse keyword's arguments and pass them to Dialog




Import mesh from FRD, VTK, VTU:  
https://lorensen.github.io/VTKExamples/site/Python/IO/ReadUnstructuredGrid/  
https://lorensen.github.io/VTKExamples/site/Python/IO/ReadLegacyUnstructuredGrid/

DistanceBetweenPoints:  
https://lorensen.github.io/VTKExamples/site/Python/SimpleOperations/DistanceBetweenPoints/

Screenshot:  
https://lorensen.github.io/VTKExamples/site/Python/Utilities/Screenshot/

Text actor for displaying model info:  
https://lorensen.github.io/VTKExamples/site/Python/GeometricObjects/TextActor/
