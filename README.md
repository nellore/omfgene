# omfgene

`omfgene.awk` is an awk script that filters a paired-end RNA-seq SAM/BAM for mates that either a) align to different chromosomes or b) align to the same chromosome > 500k bases apart. Such alignments may be evidence of gene fusions. The output of `omfgene.awk` is in the format
```
mate 1 chrom <TAB> mate 1 alignment start coordinate rounded to nearest 100 \
    <TAB> mate 2 chrom <TAB> mate 2 alignment start coordinate rounded to nearest 100
```
. We ran 
```
samtools view <a TCGA RNA-seq BAM> | mawk -f omfgene.awk | sort | uniq -c \
    | gzip ><a TCGA RNA-seq BAM>.discord.tsv.gz
```
across TCGA RNA-seq BAMs that were previously aligned to hg19 using [Seven Bridges' Cancer Genomics Cloud](https://cgc.sbgenomics.com/), where `mawk` is [a fast implementation](http://invisible-island.net/mawk/) of awk. This was done by

1. creating the Docker image omfgene using `Dockerfile` in `cgcrun/`. We ran the following sequence of commands
```
cd /path/to/omfgene/cgcrun
docker build -t biomawk .
docker run biomawk
docker login cgc-images.sbgenomics.com
docker commit $(docker ps -a | head -n2 | tail -n1 | cut -d' ' -f1) cgc-images.sbgenomics.com/anellor1/omfgene:latest
docker push cgc-images.sbgenomics.com/anellor1/omfgene:latest
```
Note the `docker login` command above required entering a username (ours was `anellor1`) and password, which was our CGC auth token.
2. using the CGC tool editor to create the tool `omfgene`, which was set up to execute the command beginning with `samtools view` above on an `m1.small` Amazon EC2 instance.
3. 
