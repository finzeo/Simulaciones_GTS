#!/bin/bash
#SBATCH --job-name=S01test # nombre para identificar el trabajo. Por defecto es el nombre del script # Autor-NumCaso-Param
#SBATCH --ntasks=80 # cantidad de cores pedidos
#SBATCH --ntasks-per-node=20 # cantidad de cores por nodo, para que agrupe o distribuya procesos
# la linea siguiente es ignorada por Slurm porque empieza con ##
##SBATCH --mem-per-cpu=4G # cantidad de memoria por core
##SBATCH --output=trabajo-%j-salida.txt # la salida y error estandar van a este archivo. Si no es especifca es slurm-%j.out (donde %j es el Job ID)
##SBATCH --error=trabajo-%j-error.txt # si se especifica, la salida de error va por separado a este archivo
#SBATCH --time=14-0 # tiempo máximo de ejecución, el formato es: dias-horas / dias-horas:minutos / horas:minutos:segundos

# aqui comienzan los comandos
# code_saturne run --param setup.xml
make
