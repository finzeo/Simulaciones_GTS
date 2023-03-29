#!/bin/bash
# define variables
# arg1: dt=0.00075
# arg2: 
# arg3: 

dtnum=$(echo "scale=5; $1" | bc -l )
((it=1))

filename=timer_stats.csv
outputfile_plain_aux=timestats_plain_aux.txt
outputfile_plain=timestats_plain.txt
outputfile_t=timestats_t.txt

# backup
#if [ -f ${filename} ]; then
#	cp ${filename} "${filename}.bk";
#fi


while read step time aux1 aux2 aux3 aux4 aux5; do
	echo "$time"
done < $filename >>$outputfile_plain_aux

# Quitar ultimo caracter de cada linea (es una coma)
sed 's/.$//' $outputfile_plain_aux > $outputfile_plain

#while read aux1 aux2 Cd; do
#	p=$(echo "scale=5; ${dtnum}*${it}" | bc -l )
#	echo "time ${p}: ${Cd}"
#	((it=it+1))
#done < $filename >>$outputfile_t
