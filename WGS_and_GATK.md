# Read mapping and SNP calling using GATK
### Getting some data
For this we're going to need a reference genome as well as some shotgun data to map to the reference and call SNPs from. We'll use the modified *B. dorsalis* PacBio reference genome for the former.

For the shotgun reads, we'll need to download these from multiple individuals, which adds repetitiveness to this task. To download from NCBI, we're going to need to use SRAtoolkit which is available [here](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)

We're going to use the data from this [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA893460/), which has 237 individuals of shotgun data for *B. dorsalis*. I've made a full list of the SRR accessions here: `/pscratch/jdu282_brazil_bootcamp2023/data/Bdor_pop_WGS_SRA_list`. Here's [the paper](https://onlinelibrary.wiley.com/doi/full/10.1111/eva.13507) that documents this Bioproject. 

Downloading these SRAs is a two step process. First we need to do a `prefetch` based on the SRR numbers, and then we actually get the fastq file using `fastq-dump`. So set up a job submission file and add these two commands:
```
for f in `cat [list of SRR accessions]`; do ./sratoolkit.3.0.7-centos_linux64/bin/prefetch $f.sra; done
for f in `cat [list of SRR accessions]`; do ./sratoolkit.3.0.7-centos_linux64/bin/fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip $f/$f.sra; done
```
Don't submit that with the full list as above. It's a pretty slow process (to do that full list, I would give it like 48 hours to be safe), and we don't need everyone downloading from NCBI at the same time. **So, as a class, let's split up the list and have everyone submit a unique set of SRA downloads!**

### Let's move on with some pre-downloaded data
While those jobs are running, let's deal with a few individuals that I've already downloaded, which are available here:`/pscratch/jdu282_brazil_bootcamp2023/data/Bdor_pop_WGS_SRA_data` You can see we have paired-end data for four individuals here, with the following SRA numbers and metadata:
```
SRR22045704	Mayotte	Dembéni	IO	Bactrocera dorsalis	-12.844351	45.166149	May_Dem_6	B104	Illumina PE150	11/28/19
SRR22045705	Grande Comore	Mandzissani	IO	Bactrocera dorsalis	-11.887344	43.406058	Com_Gra_10	B25	Illumina PE150	9/20/20
SRR22045731	Réunion	BasTerHaut	IO	Bactrocera dorsalis	-21.322553	55.485116	Reu_Bas_12	B130	Illumina PE150	2/8/19
SRR22045735	Malawi	Zomba	SSA	Bactrocera dorsalis	-15.383333	35.333333	Malaw_Zom1-10	MWZB1-MWZB10	Illumina PE150	2011
```
I pulled the metadata information from the supplementary file associated with the paper. These authors didn't do a great job of documenting specimen data with NCBI accession information, so it makes it much more difficult to match up what SRA data came from what specimen/locality/etc. **As a class, we can chat about proper protocols for data management of this type of data**

Since the SRA numbers aren't very informative, let's rename them quick before we start manipulating them at all. This is another data management technique that I use in all my projects: if specimens have species data or geographic origin data, get each specimen's sequence files renamed as early as you can in the process. This avoids the situation, for example, of using SRA accession numbers through all your data filtering steps, and so if you want to relate species/geography data to those intermediate steps, you have to cross reference some spreadsheet instead of just being able to look at the sample name/ID. This is done within reason, i.e., I don't need all that metadata information for each specimen (I can't locate a place by its lat/long coordinates off the top of my head, and some columns are the same for all specimens (e.g. Illumina PE150)). So I'd rename these specimens in this manner: `SRR22045704` becomes `Mayotte_Dembeni_B104` where B104 is the unique specimen identifier referenced in the paper. I would normally include the species name in this code (with format `[species]_[country]_[locality]_[specimenID]`), but since these are all *B. dorsalis*, it's not needed here. 

**_task_** Come up with sensible specimen names to use, soft link these data files to your own directory and rename them in the process, and finally, create a list of the four specimens using your newly created specimen IDs. NOTE, when you rename them, it is often good practice to retain the `R1` and `R2` parts of the file names, as these are commonly used to refer to read 1 and read 2 of paired-end data. 

### Subsetting these reads, to speed up the process
These are raw reads straight off the sequencer. **_task_** How many reads were sequenced per individual?

To make things fast to run through this process, let's subsample them down to just 1M reads. This will make this whole process doable in one session, but will obviously affect the final SNP dataset that we obtain at the end. Ok, on to subsetting:
```
for f in `cat [whatever you called your list of 4 individual names]`; do zcat "$f"_pass_1.fastq.gz | head -n 4000000 > $f.1Mreads.R1.fastq; done
for f in `cat [whatever you called your list of 4 individual names]`; do zcat "$f"_pass_2.fastq.gz | head -n 4000000 > $f.1Mreads.R2.fastq; done
```
Can you follow what's going on in that line? What is `zcat` doing? Why are we `head`ing that number of lines? And what does the `>` do?

Now you should have 4 sets of reads that look like this:
```
Mayotte_Dembeni_B104.1Mreads.R1.fastq
Mayotte_Dembeni_B104.1Mreads.R2.fastq
```

### Let's do some installs!
The rest of the steps for calling SNPs with GATK require a few other pieces of software, which will let us try several different ways for installing software on a cluster. We're going to go through this together, but you can always find info about how to install these on their websites, often easiest if you just google "[software] installation".

1. First up is **fastp**, which can do read trimming among other things. Fastp comes as a pre-built binary, so no installation is required. We just need to download the files to the cluster and make sure they are accessible with `chmod`.
```
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
```

2. **BWA2** is a common read mapping software. It also comes as a prebuilt binary, but this time the download is a tar.bz2 file, so we need to uncompress it.
```
https://github.com/bwa-mem2/bwa-mem2/releases
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xvf bwa-mem2-2.2.1_x64-linux.tar.bz2 
```

3. **BamUtil** is useful for modifying outputs of read mapping (bam/sam files), and this time we're going to create a conda environment to be able to run bamUtil. A conda environments is a useful way to set up one or multiple softwares that are installed with python. Each time you want to use the tools in an environment, you need to activate it similarly to how you load a module. These steps can be done in a few different ways, but here we'll create an empty environment in a specific location, then activate it, and finally install bamUtils into that environment.
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
Note, that conda is python-based tool. So you need to have python in your PATH whenever you want to use conda. Thus, in a job file, you'd need to include the `module load ccs/conda/python` part as well as (and before) your `conda activate` part. Otherwise, you'll get an error of `-bash: conda: command not found`.

4. Finally, we need **GATK** itself, which as of version 4 contains the tool picard as well (which used to be installed separately).
```
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip
unzip  gatk-4.4.0.0.zip
```
GATK uses a java wrapper to call commands, so instead of calling a .jar file, as you often do with java tools (e.g. `java --jar [file.jar]`), we can just call the wrapper as long as java is functioning (thus, we need to load the java module on MCC). Some of these steps also use samtools, so let's load that module too.
```
module load ccs/java/jdk-17.0.2
module load samtools-1.12-gcc-9.3.0-zo3utt7
```

## Let's do it!
_**data management note**_ This is a bit of a long process, so when I run it, I label my steps and all my files, so they are in-order and make sense on the cluster. I also make a directory called `RawData` where I store the raw data files or soft links to them. As you'll see, each step's output gets put in its own numbered directory. We will talk about each of these steps as we go along, but this approach leads to a directory structure like this in the end, where files (job submission file, output, and stdout/stderr) all group together by step:
```
[jdu282@mcc-login001 dorsalis_wgs]$ ls -l
total 2855570
-rw-r--r-- 1 jdu282 users       6095 Nov 22 14:19 1_fastp.17377445.err
-rw-r--r-- 1 jdu282 users          0 Nov 22 14:14 1_fastp.17377445.out
drwxr-xr-x 2 jdu282 users       4096 Nov 22 14:19 1_fastp_out
-rw-r----- 1 jdu282 users       1104 Nov 22 14:13 1_fastp.sh
-rw-r--r-- 1 jdu282 users      18644 Nov 22 14:29 2_bwa2.17377446.err
-rw-r--r-- 1 jdu282 users        429 Nov 22 14:29 2_bwa2.17377446.out
drwxr-xr-x 2 jdu282 users       4096 Nov 22 14:29 2_bwa2_out
-rw-r----- 1 jdu282 users       1458 Nov 22 14:13 2_bwa2.sh
-rw-r--r-- 1 jdu282 users      21172 Nov 22 14:31 3_picard.17377447.err
-rw-r--r-- 1 jdu282 users         68 Nov 22 14:31 3_picard.17377447.out
drwxr-xr-x 2 jdu282 users       4096 Nov 22 14:31 3_picard_out
-rw-r----- 1 jdu282 users        991 Nov 22 14:13 3_picard.sh
-rw-r--r-- 1 jdu282 users       1822 Nov 22 14:33 4_bamUtil_clip.17377448.err
-rw-r--r-- 1 jdu282 users          0 Nov 22 14:32 4_bamUtil_clip.17377448.out
drwxr-xr-x 2 jdu282 users       4096 Nov 22 14:33 4_bamUtil_clip_out
-rw-r--r-- 1 jdu282 users        980 Nov 22 14:13 4_bamUtil_clip.sh
-rw-r--r-- 1 jdu282 users          0 Nov 22 14:33 5_add_RG_bam.17377449.err
-rw-r--r-- 1 jdu282 users          0 Nov 22 14:33 5_add_RG_bam.17377449.out
drwxr-xr-x 2 jdu282 users       4096 Nov 22 14:35 5_add_RG_bam_out
-rw-r--r-- 1 jdu282 users       1060 Nov 22 14:13 5_add_RG_bam_out.sh
-rw-r--r-- 1 jdu282 users        862 Nov 22 14:13 6_gatk_hapcaller_header
drwxr-xr-x 2 jdu282 users       4096 Nov 22 15:16 6_gatk_hapcaller_out
-rw-r--r-- 1 jdu282 users        368 Nov 22 12:42 6_gatk_hapcaller_submit.sh
-rw-r--r-- 1 jdu282 users      35842 Nov 22 15:30 7_GenomicsDBImport.17377478.err
-rw-r--r-- 1 jdu282 users        738 Nov 22 15:30 7_GenomicsDBImport.17377478.out
drwxr-xr-x 8 jdu282 users       4096 Nov 22 15:30 7_GenomicsDBImport_out
-rw-r--r-- 1 jdu282 users       1215 Nov 22 15:19 7_GenomicsDBImport.sh
drwxr-xr-x 2 jdu282 users       4096 Nov 22 15:41 8_GenotypeGVCFs
-rw-r--r-- 1 jdu282 users      37105 Nov 22 15:41 8_genotypeGVCFs.17377480.err
-rw-r--r-- 1 jdu282 users        275 Nov 22 15:41 8_genotypeGVCFs.17377480.out
-rw-r--r-- 1 jdu282 users       1098 Nov 22 15:36 8_genotypeGVCFs.sh
-rw-r--r-- 1 jdu282 users       2880 Nov 22 16:04 9_merged_0miss_minDP2.vcf.log
-rw-r--r-- 1 jdu282 users     224849 Nov 22 16:04 9_merged_0miss_minDP2.vcf.recode.vcf
-rw-r--r-- 1 jdu282 users  120307671 Nov 22 15:43 9_merged.vcf
-rw-r--r-- 1 jdu282 users   14758973 Nov 22 15:41 9_merged.vcf.gz
-rw-r--r-- 1 jdu282 users     327547 Nov 22 15:41 9_merged.vcf.gz.tbi
-rw-r--r-- 1 jdu282 users     133788 Nov 22 14:19 fastp.json
-rw-r--r-- 1 jdu282 users        721 Nov 22 14:26 genome.dict
lrwxrwxrwx 1 jdu282 users         71 Nov 22 14:23 genome.fasta -> ../dorsalis_genomes/GCF_023373825.1_ASM2337382v1_chromosomes_renamed.fa
-rw-r--r-- 1 jdu282 users 1013707204 Nov 22 14:24 genome.fasta.0123
-rw-r--r-- 1 jdu282 users      17005 Nov 22 14:24 genome.fasta.amb
-rw-r--r-- 1 jdu282 users        272 Nov 22 14:24 genome.fasta.ann
-rw-r--r-- 1 jdu282 users 1647274325 Nov 22 14:26 genome.fasta.bwt.2bit.64
-rw-r--r-- 1 jdu282 users        293 Nov 22 14:26 genome.fasta.fai
-rw-r--r-- 1 jdu282 users  126713402 Nov 22 14:24 genome.fasta.pac
lrwxrwxrwx 1 jdu282 users          6 Nov 22 14:13 list -> list_4
-rw-r--r-- 1 jdu282 users         49 Nov 22 14:13 list_4
-rw-r--r-- 1 jdu282 users       2844 Nov 22 14:13 list_full
-rw-r--r-- 1 jdu282 users        324 Nov 22 15:26 list.sample_map
-rw-r--r-- 1 jdu282 users       2805 Nov 22 16:00 out.log
drwxr-xr-x 2 jdu282 users       4096 Nov 22 14:13 RawData
-rw-r--r-- 1 jdu282 users         72 Nov 22 14:32 seq.list
-rw-r--r-- 1 jdu282 users        390 Nov 22 15:34 seq.path.list
```
This approach requires a few special items in the job submission file. My job submission header for most of these files looks like this, with notes for a few of them:
```
#!/bin/bash
#SBATCH --time 2:00:00     
#SBATCH --job-name=1_fastp    # Job name, this forms the basis of what I call the job submission file, and all subsequent files
#SBATCH --nodes=1        
#SBATCH --ntasks=20       
#SBATCH --account=coa_jdu282_brazil_bootcamp2023  
#SBATCH --partition=normal
#SBATCH --mail-type ALL    
#SBATCH --mail-user julian.dupuis@uky.edu   
#SBATCH --mem=8g
#SBATCH -e %x.%j.err     # using %x in the -e and -o options applies the job name (from above) to these output files. %j does the same for job number
#SBATCH -o %x.%j.out
```
Below, I'm pasting the actual body of the job submission script for each step. In general, the options we're using here follow the "gatk best practices". See info [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices) and [here](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).

#### Step 1 fastp
Fastp does general read trimming (based on quality score, kmer analysis, etc.), and generates trimmed reads and a html output. 
```
mkdir 1_fastp_out
for f in `cat list`; do fastp --in1 ./RawData/$f.1Mreads.R1.fastq --in2 ./RawData/$f.1Mreads.R1.fastq --out1 ./1_fastp_out/$f.R1.fastpOut.fq.gz --out2 ./1_fastp_out/$f.R2.fastpOut.fq.gz  --detect_adapter_for_pe --trim_poly_g --thread 20 --overrepresentation_analysis --length_required 35 --compression 9 -h ./1_fastp_out/$f.html; done
```
Check out the html output to see what fastp can figure out for you.

#### Step 2 bwa2
BWA2, as we used before, will do the actual mapping of reads, and do it very fast for this application. We will also create a few index files for the reference genome here (hint, in my version of this, I created a soft link for the reference genome), and in the end, use samtools to order (via the sort function) the resulting sam/bam files. 
```
module load ccs/java/jdk-17.0.2
module load samtools-1.12-gcc-9.3.0-zo3utt7

bwa-mem2 index ./genome.fasta
samtools faidx ./genome.fasta
gatk CreateSequenceDictionary R=genome.fasta O=genome.dict  #really a picard tool here

mkdir 2_bwa2_out
for f in `cat list`; do echo "starting bwa on $f" ; bwa-mem2 mem -t 20 -M ./genome.fasta ./1_fastp_out/$f.R1.fastpOut.fq.gz ./1_fastp_out/$f.R2.fastpOut.fq.gz > ./2_bwa2_out/$f.paired.sam; done

module load samtools-1.12-gcc-9.3.0-zo3utt7
for f in `cat list`; do echo "starting samtools on $f"; samtools view -b ./2_bwa2_out/$f.paired.sam | samtools sort --threads 20 -o ./2_bwa2_out/$f.bam;  done
```
You can see some of my for loops have a bit more progress checking built into them (`do echo "starting samtools on $f";`). Not a necessity, but is sometimes nice to be able to track what specimens might be failing in a big list of data. 

#### step 3 Markdups
Here we're using a tool from picard (as mentioned above, now built into the gatk installation, as of version 4) to mark, and in our case, remove duplicate reads (those originating from the same template molecule. See info about these [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-). 
```
mkdir 3_picard_out
for f in `cat list`; do gatk MarkDuplicates I=./2_bwa2_out/$f.bam O=./3_picard_out/$f.dedup.bam M=./3_picard_out/$f.dupstats.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.1; done
```

#### step 4 

## Let's look at the output
### download vcftools
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
