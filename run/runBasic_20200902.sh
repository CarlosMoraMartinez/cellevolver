modelformat=".cemod"
namehead=$1
numsims=$2
maxCPU=10
free=0

#4cell_mce1inh5Xss 4cell_mce2inh5Xss 4cell_mce0
for model in 4cell_mce0fixb 4cell_mce2fixb 4cell_mce1fixb 4cell_mce0 4cell_mce2inh5 4cell_mce2inh5Xss 4cell_mce0Xss
do
    modelfile=./models_to_run/$model$modelformat
    mkdir $model
    cp $modelfile $model
    cd $model
    for i in {1..10};
    do
        while [ $free -lt 1 ]
        do
            modelrunning=$(ps -aef | pgrep -f cellEvolver.py | wc -l)  #how many CPUs are running?
            vertex=$(ps -aef | pgrep vertex | wc -l)
            cpuRUNNING=$(($modelrunning + $vertex))
            echo '  Counting CPUs...: '$cpuRUNNING
            if [ $cpuRUNNING -lt $maxCPU ]
            then
                newdir="$namehead$model"s"$i"
                echo '      CPUs lower than '$maxCPU
                echo '      sending job: '$newdir
                mkdir $newdir
                cp ../*.py $newdir
                cd $newdir
                nohup python3 cellEvolver.py $newdir ../$model$modelformat >$newdir'.out' 2>newdir'.err' &
                cd ../
                free=$(($free + 1))
                sleep 3
            else
                sleep 30
            fi
        done
        free=0
    done
    cd ../
done
#nohup python cellEvolver.py 180530_4cell_mce0Xaa_try1 4cell_mce0Xaa.cemod
#nohup python cellEvolver.py 180530_4cell_mce0Xs2_try1 4cell_mce0Xs.cemod
#nohup python cellEvolver.py 180530_4cell_mce0Xv2_try1 4cell_mce0Xv.cemod
