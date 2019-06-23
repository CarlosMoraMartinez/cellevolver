#!/bin/sh

modelformat=".cemod"
namehead="181101set"
for model in 4cell_mce4fix 4cell_mce5fix 4cell_mce8fix 4cell_mce9fix 6cell_mce6fix 6cell_mce7fix 4cell_mce0fixIntb
do
    modelfile=../../models_to_run/$model$modelformat
    mkdir $model
    cd $model
    echo $model " "
    sbatch --array=0-200 --cpus-per-task 24 --mincpus 4 ../runsinglesim.sh $namehead $model $modelfile
    cd ../
done

