#!/bin/bash
# define variables
filename=torques_CBRao.txt
outputfile_plain=torques_CBRao_plain.txt
outputfile_t=torques_CBRao_t.txt
dt=0.00075

dtnum=$(echo "scale=5; ${dt}" | bc -l )
((it=1))

# backup
cp ${filename} "${filename}.bk"

while read aux1 aux2 aux3 mx aux5 my aux7 mz aux9; do
    echo "$mx $my $mz"
done < $filename >>$outputfile_plain

while read aux1 aux2 aux3 mx aux5 my aux7 mz aux9; do
	p=$(echo "scale=5; ${dtnum}*${it}" | bc -l )
	echo "time ${p}: $mx $my $mz"
	((it=it+1))
done < $filename >>$outputfile_t
