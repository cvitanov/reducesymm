#!/bin/bash

python vfield2vtk.py > tangents.vtk
python gorbitone2vtk.py > gorbit1.vtk
python gorbittwo2vtk.py > gorbit2.vtk
python gorbitthree2vtk.py > gorbit3.vtk
python slicehplane2vtk.py > slicehplane.vtk

rm *.pyc
