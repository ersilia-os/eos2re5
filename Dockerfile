FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN sudo apt update
RUN sudo apt install python2.7
RUN sudo apt install gfortran-7
RUN conda create -n eos2re5-py27 python=2.7 -y
RUN conda install -n eos2re5-py27 -c rdkit rdkit=2018.09.3 -y
RUN conda install -n eos2re5-py27 -c conda-forge openbabel=3.0.0 -y 
RUN conda install -n eos2re5-py27 -c rmg mopac=2017 -y
RUN conda install -n eos2re5-py27 pip -y
RUN conda activate eos2re5-py27 && pip install scikit-learn==0.17.1 && pip install scipy==1.2.3 && conda deactivate


WORKDIR /repo
COPY . /repo
