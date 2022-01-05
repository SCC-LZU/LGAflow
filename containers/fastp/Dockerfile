FROM alpine:latest

# pre-built binary
RUN wget http://opengene.org/fastp/fastp && \
    cp ./fastp /usr/local/bin/fastp && \
    chmod a+x ./fastp
# https://github.com/OpenGene/fastp#get-fastp