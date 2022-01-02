FROM alpine

RUN apk add --update perl
RUN wget https://download.fastgit.org/EVidenceModeler/EVidenceModeler/archive/refs/tags/v1.1.1.tar.gz -O evidencemodeler-v1.1.1.tar.gz &&\
    tar -xvf evidencemodeler-v1.1.1.tar.gz && \
    mv EVidenceModeler-1.1.1 /usr/local/bin/
ENV EVM_HOME /usr/local/bin/EVidenceModeler-1.1.1
