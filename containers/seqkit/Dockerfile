FROM alpine:latest

RUN apk add --upgrade bash
RUN wget http://app.shenwei.me/data/seqkit/seqkit_linux_amd64.tar.gz && \
    tar -xvf seqkit_linux_amd64.tar.gz && \
    mv seqkit /usr/local/bin/