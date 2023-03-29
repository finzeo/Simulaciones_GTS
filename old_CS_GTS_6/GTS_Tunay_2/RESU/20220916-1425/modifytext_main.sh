#!/bin/bash
# define variables
# arg1: slurmname="slurm-140926.out"
# arg2: keyword_name="blank"
# arg3: dt=0.00075

cat $1 | grep force > "forces_$2.txt"
cat $1 | grep torque > "torques_$2.txt"
cat $1 | grep Drag > "Cd_$2.txt"
cat $1 | grep "Lateral coefficient" > "Clat_$2.txt"
cat $1 | grep Lift > "Cz_$2.txt"

./modifytext_forces.sh $1 $2 $3
./modifytext_torques.sh $1 $2 $3
./modifytext_Cd.sh $1 $2 $3
./modifytext_Clat $1 $2 $3
./modifytext_Cz $1 $2 $3



