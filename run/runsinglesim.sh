#!/bin/sh
#1: namehead
#2: model
#3: modelfile
module load Python/3.6.1-foss-2016b
newdir="$1$2"s"$SLURM_ARRAY_TASK_ID"
mkdir $newdir
cp ../*.py $newdir
cd $newdir
python cellEvolver.py $newdir $3

