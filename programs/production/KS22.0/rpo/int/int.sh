#!/bin/bash
jmin=100
jmax=156
j=0
ls ../ | grep ks22rpo | while read i
do
	echo "$i"
	let j=j+1
	echo "$j"
	if [ "$j" -le "$jmin" ]; then
        	nice ./ks "../$i/" 
	fi
done

