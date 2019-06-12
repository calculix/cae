© Ihor Mirzov, UJV Rez, June 2019.  
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
    sudo apt install qttools5-dev-tools

Some examples for getting started with VTK in Python:

    https://lorensen.github.io/VTKExamples/site/Python/

<br/><br/>



# TODO

treeView context menu: Show/Hide empty containers, Expand/Collapse all.

Delete implementation from treeView. Edit implementation.

Get all mutually exclusive arguments from the manual - now they are supported.

Styles with icons for treeView elements.

Parse INP-file into DOM during mesh import. Model import with all keywords.


Import mesh from FRD, VTK, VTU:  
https://lorensen.github.io/VTKExamples/site/Python/IO/ReadUnstructuredGrid/  
https://lorensen.github.io/VTKExamples/site/Python/IO/ReadLegacyUnstructuredGrid/

DistanceBetweenPoints:  
https://lorensen.github.io/VTKExamples/site/Python/SimpleOperations/DistanceBetweenPoints/

Screenshot:  
https://lorensen.github.io/VTKExamples/site/Python/Utilities/Screenshot/

Text actor for displaying model info:  
https://lorensen.github.io/VTKExamples/site/Python/GeometricObjects/TextActor/

Camera: perspective mode, parallel mode, projections, reset (rotate to initial position + fitView)
