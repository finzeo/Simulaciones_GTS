#!/bin/bash
# define variables
filename=forces_CBRao.txt
outputfile_plain=forces_CBRao_plain.txt
outputfile_t=forces_CBRao_t.txt
dt=0.00075

dtnum=$(echo "scale=5; ${dt}" | bc -l )
((it=1))

# backup
cp ${filename} "${filename}.bk"

while read aux1 fx fy fz; do
    echo "$fx $fy $fz"
done < $filename >>$outputfile_plain

while read aux1 fx fy fz; do
	p=$(echo "scale=5; ${dtnum}*${it}" | bc -l )
	echo "time ${p}: $fx $fy $fz"
	((it=it+1))
done < $filename >>$outputfile_t

