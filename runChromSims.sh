modelformat=".cemod"
namehead="180601_"
for model in 4cell_mce1inh5Xss 4cell_mce2inh5Xss
do
    modelfile=~/cellevolver/chromatin_trials/models_to_run/$model$modelformat
    mkdir $model
    cd $model
    for i in {1..10};
    do
        newdir="$namehead$model"s"$i"
        echo $newdir
        mkdir $newdir
        cp ../*.py $newdir
        cd $newdir
        nohup python cellEvolver.py $newdir $modelfile
        echo "***********************"$model"____"$i
        cd ../
    done
    cd ../
done
#nohup python cellEvolver.py 180530_4cell_mce0Xaa_try1 4cell_mce0Xaa.cemod
#nohup python cellEvolver.py 180530_4cell_mce0Xs2_try1 4cell_mce0Xs.cemod
#nohup python cellEvolver.py 180530_4cell_mce0Xv2_try1 4cell_mce0Xv.cemod
