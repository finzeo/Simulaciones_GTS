#!/bin/bash
# define variables
# arg1: slurm-name="slurm-140926.out"
# arg2: keyword_name="blank"
# arg3: dt=0.00075

dtnum=$(echo "scale=5; $3" | bc -l )
((it=1))

filename=Cd_$2.txt
outputfile_plain=Cd_$2_plain.txt
outputfile_t=Cd_$2_t.txt

# backup
#if [ -f ${filename} ]; then
#	cp ${filename} "${filename}.bk";
#fi


while read aux1 aux2 Cd; do
	echo "$Cd"
done < $filename >>$outputfile_plain

while read aux1 aux2 Cd; do
	p=$(echo "scale=5; ${dtnum}*${it}" | bc -l )
	echo "time ${p}: ${Cd}"
	((it=it+1))
done < $filename >>$outputfile_t
