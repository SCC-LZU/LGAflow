FROM ubuntu:latest

WORKDIR /root

SHELL ["bash", "-c"]
RUN apt update && apt install -y wget python3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
chmod +x Miniconda3-latest-Linux-x86_64.sh && \
./Miniconda3-latest-Linux-x86_64.sh -b -f
ENV PATH=/root/miniconda3/bin:$PATH
RUN conda install -c conda-forge -n base mamba
ADD ./busco.yaml busco.yaml
RUN mamba env create -f busco.yaml
ENV PATH /root/miniconda3/envs/busco/bin:$PATH
