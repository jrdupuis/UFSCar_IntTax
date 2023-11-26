# Read mapping and SNP calling using GATK
## Getting some data
For this we're going to need a reference genome as well as some shotgun data to map to the reference and call SNPs from. We'll use the modified B. dorsalis PacBio reference genome for the former.

For the shotgun reads, we'll need to download these from multiple individuals, which adds repetitiveness to this task. To download from NCBI, we're going to need to use SRAtoolkit which is available [here](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)

Here's a list of SRA accessions that we can deal with 
```
head SRR_list
```

Downloading these SRAs is a two step process. First we need to do a `prefetch` based on the SRR numbers, and then we actually get the fastq file using `fastq-dump`

```
for f in `cat SRR_list`; do ./sratoolkit.3.0.7-centos_linux64/bin/prefetch $f.sra; done
for f in `cat SRR_list`; do ./sratoolkit.3.0.7-centos_linux64/bin/fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip $f/$f.sra; done
```

To make things fast to run through this process, let's deal with just four individuals, and subsample them down to just 1M reads

```
for f in `cat list_4`; do zcat "$f"_pass_1.fastq.gz | head -n 4000000 > $f.1Mreads.R1.fastq; done
for f in `cat list_4`; do zcat "$f"_pass_2.fastq.gz | head -n 4000000 > $f.1Mreads.R2.fastq; done
```
Can you follow what's going on in that line? Why are we `head`ing that number of lines?
Now you should have 4 sets of reads that look like this:
```
SRR22045704.1Mreads.R1.fastq
SRR22045704.1Mreads.R2.fastq
```

The rest of the steps for calling SNPs with GATK require a few other pieces of software, which will let us try several different ways for installing software on a cluster.

First up is fastp, which can do read trimming among other things. Fastp comes as a pre-built binary, so no installation is required. We just need to download the files to the cluster and make sure they are accessible with `chmod`.
```
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
```

BWA2 is a common read mapping software. It also comes as a prebuilt binary, but this time the download is a tar.bz2 file, so we need to uncompress it.
```
https://github.com/bwa-mem2/bwa-mem2/releases
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xvf bwa-mem2-2.2.1_x64-linux.tar.bz2 
```

BamUtil is useful for modifying outputs of read mapping (bam/sam files), and this time we're going to create a conda environment to be able to run bamUtil. A conda environments is a useful way to set up one or multiple softwares that are installed with python. Each time you want to use the tools in an environment, you need to activate it similarly to how you load a module. These steps can be done in a few different ways, but here we'll create an empty environment in a specific location, then activate it, and finally install bamUtils into that environment.
```
module load ccs/conda/python
conda create --prefix /scratch/jdu282/bamUtil_env python=3.9
conda activate /scratch/jdu282/bamUtil_env
conda install -c bioconda bamutil
```
When you have a conda environment activated, the shell prompt will change to show that. It's often useful to check what all is installed in a specific environment and we can do that with the following command:
```
conda list
```

Finally, we need GATK itself, which as of version 4 contains the tool picard as well (which used to be installed separately).
```
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip
unzip  gatk-4.4.0.0.zip
```
GATK uses a java wrapper to call commands, so instead of calling a .jar file, as you often do with java tools (e.g. `java --jar [file.jar]`), we can just call the wrapper as long as java is functioning (thus, we need to load the java module on MCC). Some of these steps also use samtools, so let's load that module too.
```
module load ccs/java/jdk-17.0.2
module load samtools-1.12-gcc-9.3.0-zo3utt7
```

module load ccs/conda/python

	## step 6 takes ~30 min with 1M reads


## download vcftools
VCFtools is a useful software for filtering and manipulating vcf files. It doesn't come as a pre-built binary, so we're going to have to install it like in the old days using configure, make, and make install. First let's download the software using git clone.
```
git clone https://github.com/vcftools/vcftools.git
```
Now follow the instructions in the vcftools documentation.
```
./autogen.sh 
./configure --prefix=/scratch/jdu282/vcftools
make
make install
```
Did it work? Sometimes it can be hard to figure out what error messages mean during installs. But in the error message above you should probably see a couple things mentioned: pkgconf and zlib. These are various compilers and other nitty gritty computer science dependencies of vcftools. Once you get used to seeing some of these names, you can often figure out if a cluster has other versions of the dependencies that might work better with vcftools. Let's search around and see if anything works.
```
module spider pkgconf
module spider zlib
module load pkgconf-1.7.4-intel-2021.3.0-eangrvu
module load zlib-1.2.11-gcc-8.4.1-b4szoqs
```

# For dealing with 4 individual, 1M reads datasets
/scratch/jdu282/vcftools/bin/vcftools --vcf 9_merged.vcf --remove-indels --max-missing 1.0 --minDP 2 --recode --out 9_merged_0miss_minDP2.vcf
