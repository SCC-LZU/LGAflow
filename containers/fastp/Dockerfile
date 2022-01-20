FROM alpine:latest

# pre-built binary
RUN apk add --upgrade bash
RUN wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp && \
    cp ./fastp /usr/local/bin/fastp
    
# https://github.com/OpenGene/fastp#get-fastp