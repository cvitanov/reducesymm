#!/bin/bash
jmin=100
#jmax=10
j=0
ls | grep ks22rpo | while read i
do
	echo "$i"
	let j=j+1
	echo "$j"
	if [ "$j" -gt "$jmin" ]; then
        	nice ./ks $i 
	fi
done

