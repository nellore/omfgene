# Specifying the ubuntu base image
FROM ubuntu
    
# Name and email of the person who maintains the file
MAINTAINER Rosalind Franklin rosalind.franklin@sbgenomics.com
  
# Set working directory as "/"
WORKDIR /
  
# Updating ubuntu and installing other necessary software
RUN apt-get update --yes \
&& apt-get install wget build-essential zlib1g-dev libncurses5-dev vim --yes
   
# Pulling SAMTools from its repository, unpacking the archive and installing
RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 \
&& tar jxf samtools-1.2.tar.bz2 \
&& cd samtools-1.2 \
&& make \
&& make install
  

RUN wget http://invisible-island.net/datafiles/release/mawk.tar.gz \
&& tar xvzf mawk.tar.gz \
&& cd mawk* \
&& ./configure \
&& make \
&& PATH=/data/mawk-1.3.4-20161120:$PATH

RUN bash -c 'mkdir -pv /data/cgc_outputs' && chmod 777 /data/cgc_outputs/ && cd /data/cgc_outputs && wget "https://raw.githubusercontent.com/nellore/omfgene/master/omfgene.awk" 
VOLUME ["/data/cgc_outputs"]
