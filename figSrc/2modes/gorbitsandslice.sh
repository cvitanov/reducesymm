#!/bin/bash

python gtans2vtk.py > paraview/gorbitsandslice/tangents.vtk
python gorbitone2vtk.py > paraview/gorbitsandslice/gorbit1.vtk
python gorbittwo2vtk.py > paraview/gorbitsandslice/gorbit2.vtk
python gorbitthree2vtk.py > paraview/gorbitsandslice/gorbit3.vtk
python slicehplane2vtk.py > paraview/gorbitsandslice/slicehplane.vtk

rm *.pyc
