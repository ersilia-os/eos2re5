source $CONDA_PREFIX_1/etc/profile.d/conda.sh
conda activate eos2re5-py27
python $1/code/python_2_main.py $2 $3
conda deactivate
