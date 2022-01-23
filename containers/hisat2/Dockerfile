FROM ubuntu:latest

RUN apt update && apt install -y wget unzip python 
RUN cd /usr/local && \
    wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1.zip && \
    unzip hisat2-2.2.1.zip
ENV PATH /usr/local/hisat2-2.2.1:$PATH

#build samtools
RUN apt update && apt install -y build-essential wget python3
RUN apt install -y libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev
RUN wget https://download.fastgit.org/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && \
    tar -xjf samtools-1.14.tar.bz2 && \
    cd samtools-1.14 && \
    ./configure --prefix /usr/local && \
    make -j$(nproc) && \
    make install
