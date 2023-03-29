#!/bin/bash
# define variables
# arg1: slurm-name="slurm-140926.out"
# arg2: keyword_name="blank"
# arg3: dt=0.00075

dtnum=$(echo "scale=5; $3" | bc -l )
((it=1))

filename=forces_$2.txt
outputfile_plain=forces_$2_plain.txt
outputfile_t=forces_$2_t.txt

# backup
#cp ${filename} "${filename}.bk"

while read aux1 fx fy fz; do
    echo "$fx $fy $fz"
done < $filename >>$outputfile_plain

while read aux1 fx fy fz; do
	p=$(echo "scale=5; ${dtnum}*${it}" | bc -l )
	echo "time ${p}: $fx $fy $fz"
	((it=it+1))
done < $filename >>$outputfile_t

