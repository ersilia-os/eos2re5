source $CONDA_PREFIX_1/etc/profile.d/conda.sh
CONDA_PATH=$(dirname $(dirname $(which conda)))
PYTHON_ENV_PATH="${CONDA_PATH}/envs/eos2re5-py27/bin/python"
$PYTHON_ENV_PATH $1/code/main.py $2 $3
