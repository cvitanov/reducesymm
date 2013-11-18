#!/bin/bash

python setparametersSO2.py
python onslicesolver.py > data/solutiononslice.dat
python onslice3dplotter.py
