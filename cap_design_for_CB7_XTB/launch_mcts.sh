#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=48:0:0
#$ -l mem=3G
#$ -l tmpfs=10G
#$ -N mcxtb
#$ -pe smp 4
#$ -t 1-20
number=$SGE_TASK_ID
# paramfile=$(pwd)/paramfile
# index=$(sed -n ${number}p $paramfile | awk '{print $1}')
# infile=$(sed -n ${number}p $paramfile | awk '{print $2}')
# export GAUSS_EXEDIR="/home/uccahcl/g16"
# export GAUSS_SCRDIR="/home/uccahcl/Scratch/FIND_CAP/GAUSS_SCRATCH"

#module load gaussian/g16-a03/pgi-2016.5
#mkdir -p $GAUSS_SCRDIR

export KMP_INIT_AT_FORK=FALSE
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OMP_STACKSIZE=4G


~/anaconda2/envs/chemts/bin/python mcts_ligand.py $number
