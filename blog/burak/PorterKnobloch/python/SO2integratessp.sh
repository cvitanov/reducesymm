#!/bin/bash

python setparametersSO2.py
python sspsolver.py > data/sspsolution.dat
python ssp3dplotter.py
