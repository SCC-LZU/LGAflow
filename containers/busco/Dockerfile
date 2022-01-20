FROM ubuntu:latest

WORKDIR /root

SHELL ["bash", "-c"]
RUN apt update && apt install -y wget python3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
chmod +x Miniconda3-latest-Linux-x86_64.sh && \
./Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local/miniconda3
ENV PATH=/usr/local/miniconda3/bin:$PATH
RUN conda install -c conda-forge -n base mamba
RUN mamba install -y -c conda-forge -c bioconda -n base busco=3.0.2
RUN rm -f Miniconda3-latest-Linux-x86_64.sh
ENV RSTUDIO_WHICH_R=/usr/local/miniconda3/bin/R \
    AUGUSTUS_SCRIPTS_PATH=/usr/local/miniconda3/bin/ \
    AUGUSTUS_CONFIG_PATH=/usr/local/miniconda3/config/ \
    AUGUSTUS_BIN_PATH=/usr/local/miniconda3/bin/



