#!/bin/bash
jmin=0
jmax=28
j=0
ls ../ | grep ks22rpo | while read i
do
	echo "$i"
	let j=j+1
	echo "$j"
	if [ "$j" -lt "$jmax" ]; then
        	nice ./ks "../$i/" 
	fi
done

