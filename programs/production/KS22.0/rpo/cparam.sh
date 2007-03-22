#!/bin/bash
jmin=9
jmax=10
j=0
ls | grep ks22rpo | while read i
do
	echo "$i"
	let j=j+1
	echo "$j"
	if [ "$j" -gt "$jmin" ]; then
        	cp  'ks22rpo040.14_03.700/parameters.dat' $i/ 
	fi
done

