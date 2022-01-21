FROM ubuntu:latest

RUN apt update && apt install -y wget unzip python
RUN cd /usr/local && \
    wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1.zip && \
    unzip hisat2-2.2.1.zip
ENV PATH /usr/local/hisat2-2.2.1:$PATH
