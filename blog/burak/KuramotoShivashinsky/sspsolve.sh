#!/bin/bash

octave solveKS.m
sed -i.bak -e '1,5d' data/sspsolution.dat #Delete first 5 lines to make it readable from python
#python sspsolver.py > data/sspsolution.dat
#python ssp3dplotter.py
