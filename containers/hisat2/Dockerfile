FROM alpine:latest

RUN apk add --upgrade perl python2
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1.zip && \
    unzip hisat2-2.2.1.zip && \
    rm -f hisat2-2.2.1.zip
ENV PATH $PWD/hisat2-2.2.1/:$PATH
