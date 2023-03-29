#!/bin/bash
# define variables
filename=forces_CBRao.txt
outputfile_plain=Clat_CBRao_plain.txt
outputfile_t=Clat_CBRao_t.txt
dt=0.00075
factor=7.435356812

dtnum=$(echo "scale=5; ${dt}" | bc -l )
((it=1))

# backup
cp ${filename} "${filename}.bk"

while read aux1 fx fy fz; do
	clat=$(echo "scale=7; ${factor}*${fy}" | bc -l )
	echo "$fy"
done < $filename >>$outputfile_plain

while read aux1 fx fy fz; do
	p=$(echo "scale=5; ${dtnum}*${it}" | bc -l )
	clat=$(echo "scale=7; ${factor}*${fy}" | bc -l )
	echo "time ${p}: $clat"
	((it=it+1))
done < $filename >>$outputfile_t
