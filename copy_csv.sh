
#4cell_mce4fix 4cell_mce5fix 4cell_mce8fix 4cell_mce9fix 
#6cell_mce6fix 6cell_mce7fix 4cell_mce0fixIntc 4cell_mce0fixIntb
#4cell_mce2fixIntc 4cell_mce2fixIntb 4cell_mce10fixIntb
mkdir sim_tables
for model in 4cell_mce4fix 4cell_mce5fix 4cell_mce8fix 4cell_mce9fix 6cell_mce6fix 6cell_mce7fix
do
    modelfile=../../models_to_run/$model$modelformat
    mkdir sim_tables/$model
    cp $model/*/*.csv sim_tables/$model
    echo $model " "
done

