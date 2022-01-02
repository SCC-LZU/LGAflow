FROM ubuntu:20.04

RUN apt update && apt install -y build-essential wget
RUN apt install -y python3 python3-pip perl && \
    pip install h5py && \
    export PERL_MM_USE_DEFAULT=1
# intall rmblast
RUN wget http://www.repeatmasker.org/rmblast-2.11.0+-x64-linux.tar.gz && \
    tar -xvf rmblast-2.11.0+-x64-linux.tar.gz && \
    mv rmblast-2.11.0 /usr/local/ && \
    mv /usr/local/rmblast-2.11.0 /usr/local/rmblast

# install TRF
COPY ./trf409.linux64 /root

# install repeatmasker
RUN wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz && \
    mv RepeatMasker-4.1.2-p1.tar.gz /usr/local && \
    cd /usr/local && \
    tar -xvf RepeatMasker-4.1.2-p1.tar.gz && \
    rm RepeatMasker-4.1.2-p1.tar.gz
# using latest Dfam
RUN cd /usr/local/RepeatMasker/Libraries && \
    wget https://www.dfam.org/releases/current/families/Dfam.h5.gz && \
    gunzip Dfam.h5.gz && \
    rm Dfam.h5.gz
# using Repbase
# COPY ./RepBaseRepeatMaskerEdition-20181026.tar.gz /usr/local/RepeatMasker/
# RUN cd /usr/local/RepeatMasker && \
#     gunzip RepBaseRepeatMaskerEdition-20181026.tar.gz && \
#     tar xvf RepBaseRepeatMaskerEdition-20181026.tar && \
#     rm RepBaseRepeatMaskerEdition-20181026.tar
RUN cd /usr/local/RepeatMasker && \
    perl ./configure -libdir /usr/local/RepeatMasker/Libraries/ -rmblast_dir /usr/local/rmblast -trf_prgm /root/trf409.linux64 -default_search_engine rmblast
# install RepeatModeler
RUN cd /root && \
    wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz && \
    tar -xvf RECON-1.08.tar.gz && \
    cd /root/RECON-1.08/src && \
    make && \
    make install && \
    sed -i '3c\ $path = /root/RECON-1.08/bin' ../script/recon.pl && \
    mv /root/RECON-1.08 /usr/local/RECON && \
    cd /root && \
    wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz && \
    tar -xvf RepeatScout-1.0.6.tar.gz && \
    cd /root/RepeatScout-1.0.6 && \
    make && \
    make install && \
    cd /root && \
    wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.2a.tar.gz && \
    tar -xvf RepeatModeler-2.0.2a.tar.gz && \
    mv RepeatModeler-2.0.2a /usr/local/ && \
    cd /usr/local/RepeatModeler-2.0.2a && \
    cpan JSON File::Which URI LWP && \
ENV PATH /usr/local/RepeatModeler-2.0.2a:/usr/local/RepeatMasker:/usr/local/rmblast/bin:$PATH

