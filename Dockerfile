FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN sudo apt update
RUN sudo apt install python2.7

RUN python2.7 get-pip.py
RUN python2.7 -m pip install scikit-learn==0.17.1

RUN conda create -n py27 python=2.7
RUN conda install -n py27 -c rdkit rdkit=2018.09.3
RUN conda install -n py27 -c anaconda scipy

WORKDIR /repo
COPY . /repo
