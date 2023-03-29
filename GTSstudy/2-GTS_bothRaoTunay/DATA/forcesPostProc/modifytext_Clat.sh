#!/bin/bash
# define variables
# arg1: slurm-name="slurm-140926.out"
# arg2: keyword_name="blank"
# arg3: dt=0.00075

#factor=7.435356812

filename=Clat_$2.txt
outputfile_plain=Clat_$2_plain.txt
outputfile_t=Clat_$2_t.txt

dtnum=$(echo "scale=5; $3" | bc -l )
((it=1))

# backup
#cp ${filename} "${filename}.bk"

while read aux1 aux2 Clat; do
	echo "$Clat"
done < $filename >>$outputfile_plain

while read aux1 aux2 Clat; do
	p=$(echo "scale=5; ${dtnum}*${it}" | bc -l )
	echo "time ${p}: ${Clat}"
	((it=it+1))
done < $filename >>$outputfile_t

