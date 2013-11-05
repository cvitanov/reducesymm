#!/bin/bash

python setparametersO2.py
python comovsolver.py > data/solutioncomov.dat
python comov3dplotter.py
