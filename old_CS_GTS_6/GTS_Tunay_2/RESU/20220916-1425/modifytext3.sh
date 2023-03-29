#!/bin/bash
# define variables
keyword_name="blank"
dt=0.00075

dtnum=$(echo "scale=5; ${dt}" | bc -l )
((it=1))

filename=forces_${keyword_name}.txt
outputfile_plain=forces_${keyword_name}_plain.txt
outputfile_t=forces_${keyword_name}_t.txt

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

