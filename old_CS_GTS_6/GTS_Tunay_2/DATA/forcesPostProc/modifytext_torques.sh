#!/bin/bash
# define variables
# arg1: slurm-name="slurm-140926.out"
# arg2: keyword_name="blank"
# arg3: dt=0.00075

dtnum=$(echo "scale=5; $3" | bc -l )
((it=1))

filename=torques_$2.txt
outputfile_plain=torques_$2_plain.txt
outputfile_t=torques_$2_t.txt

# backup
#cp ${filename} "${filename}.bk"

while read aux1 aux2 aux3 mx aux5 my aux7 mz aux9; do
    echo "$mx $my $mz"
done < $filename >>$outputfile_plain

while read aux1 aux2 aux3 mx aux5 my aux7 mz aux9; do
	p=$(echo "scale=5; ${dtnum}*${it}" | bc -l )
	echo "time ${p}: $mx $my $mz"
	((it=it+1))
done < $filename >>$outputfile_t
