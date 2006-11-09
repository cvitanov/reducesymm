#!/bin/bash
if [ $# -ne 1 ]; then
	echo 1>&2 Usage: $0 directory
	exit 127
fi
d0=`pwd`
shopt -s expand_aliases
source ~/.Mathematica/bash_init
cd $1; nohup math -noprompt -run "<<$d0/invokePlotEquil.m" > mathOut & 
#cd "$d0"