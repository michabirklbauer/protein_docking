# Dockerfile for running Streamlit App
# author: Micha Birklbauer

FROM ubuntu:20.04

LABEL maintainer="micha.birklbauer@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    libopenbabel-dev \
    libopenbabel6 \
    openbabel \
    python3-openbabel \
    pymol \
    python3-pymol \
    python3-distutils \
    python3-lxml \
    python3-rdkit \
    python3-pip

RUN pip3 install numpy
RUN pip3 install plip --no-deps
RUN pip3 install biopandas
RUN pip3 install matplotlib
RUN pip3 install streamlit

COPY streamlit_app.py .

CMD  ["streamlit", "run", "streamlit_app.py"]
