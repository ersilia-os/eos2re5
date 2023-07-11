FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN sudo apt update -y
RUN sudo apt install python2.7 -y
RUN sudo apt install gfortran-7 -y
RUN conda create -c conda-forge -n eos2re5-py27 python=2.7 -y
RUN conda install -n eos2re5-py27 -c rdkit rdkit=2018.09.3 -y
RUN wget https://anaconda.org/conda-forge/openbabel/3.0.0/download/linux-64/openbabel-3.0.0-py27hdef5451_1.tar.bz2
RUN conda install -n eos2re5-py27 openbabel-3.0.0-py27hdef5451_1.tar.bz2 -y
RUN conda install -n eos2re5-py27 -c conda-forge mopac -y
RUN conda install -n eos2re5-py27 pip -y
RUN conda init bash && source ~/.bashrc && conda activate eos2re5-py27 && pip install scikit-learn==0.17.1 && pip install scipy==1.2.3 && conda deactivate

WORKDIR /repo
COPY . /repo
