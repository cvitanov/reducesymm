#!/bin/bash

python setparametersSO2.py
python invpolsolver.py > data/invpolsolution.dat
python invpol3dplotter.py
