FROM chrishah/busco-docker:latest
RUN apt-get remove -y ncbi-blast+
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz && \
tar -xzvf ncbi-blast-2.2.30+-x64-linux.tar.gz && \
cp -f ncbi-blast-2.2.30+/bin/* /usr/bin/




