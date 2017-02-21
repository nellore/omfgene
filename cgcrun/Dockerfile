FROM biocontainers/biocontainers:latest
RUN conda install samtools=1.3.1
RUN curl -O http://invisible-island.net/datafiles/release/mawk.tar.gz; tar xvzf mawk.tar.gz; cd mawk*; ./configure; make; PATH=/data/mawk-1.3.4-20161120:$PATH
