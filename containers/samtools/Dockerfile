FROM ubuntu:20.04


#compile from source
RUN apt update && apt install -y build-essential wget
RUN apt install -y libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev
RUN wget https://download.fastgit.org/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && \
    tar -xjf samtools-1.14.tar.bz2 && \
    cd samtools-1.14 && \
    ./configure && \
    make && \
    make install