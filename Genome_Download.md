# Downloading and manipulating a genome
Let's download a high quality reference genome for B. dorsalis, and then we'll do some basic editing of it. We'll download [this one](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023373825.1/).
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/373/825/GCF_023373825.1_ASM2337382v1/GCF_023373825.1_ASM2337382v1_genomic.fna.gz
```
One of the tools we can use for manipulating this genome file will require an uncompressed version, so let's do that too.
```
gunzip -c GCF_023373825.1_ASM2337382v1_genomic.fna.gz > GCF_023373825.1_ASM2337382v1_genomic.fna
```
Now, let's download setk
