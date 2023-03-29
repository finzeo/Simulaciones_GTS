#!/bin/bash

# Define values
dt=0.00075
nt=63000
filename=forcesnew2.txt

dtnum=$(echo "scale=5; ${dt}" | bc -l )
#echo $dtnum
#let tf=47.25
cp ${filename} "${filename}.bk"
for (( k=1; k<=${nt}; k++ ))
do
	#echo $k
	#((t="k/8"))
	p=$(echo "scale=5; ${dtnum}*${k}" | bc -l )
	sed -i "${k}s/^/time ${p}s: /" ${filename}
done
