FROM alpine:latest
RUN apk add --upgrade bash
RUN wget https://download.fastgit.org/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
    mv bedtools.static.binary /usr/local/bin/bedtools && \
    chmod a+x /usr/local/bedtools