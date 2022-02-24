FROM ubuntu:latest

RUN apt update && apt install -y perl wget build-essential gawk
RUN wget https://download.fastgit.org/EVidenceModeler/EVidenceModeler/archive/refs/tags/v1.1.1.tar.gz -O evidencemodeler-v1.1.1.tar.gz &&\
    tar -xvf evidencemodeler-v1.1.1.tar.gz && \
    mv EVidenceModeler-1.1.1 /usr/local/bin/
RUN export PERL_MM_USE_DEFAULT=1 && \
cpan URI::Escape
RUN cpan DBI
ENV EVM_HOME /usr/local/bin/EVidenceModeler-1.1.1
