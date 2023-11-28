## Analysing a GBS dataset with Stacks

Let's take a look at this [bioproject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA987164) from [this paper](https://www.sciencedirect.com/science/article/abs/pii/S1055790323001926).<br>
What about the metadata for this dataset? If you click on the SRA link (link to 798 SRA files), and then click on "send to" > "file" > and "run info", that opens up the information associated with each specimen through their BioSample and SRA submission. But not all pertinent info is there (country of origin, collection locality, etc.). 
But if you go to the publication and look at Table S1, the authors included SRR numbers and pertinent locality information, so we're in luck.

This data is downloaded (just as we did before with sratoolkit) and is here `/pscratch/jdu282_brazil_bootcamp2023/data/Bdor_GBS`. So create a directory for this analysis in your `scratch` or in your named dir in `pscratch`, and let's soft link all those files into your directory of choice, so we don't have to bog down NCBI with duplicate downloads. Here's a quick way to do so:
```
for f in `ls /pscratch/jdu282_brazil_bootcamp2023/data/Bdor_GBS`; do ln -s /pscratch/jdu282_brazil_bootcamp2023/data/Bdor_GBS/$f; done
```
Can you follow what's going on there? If not, break it down to what is happening in the first part (`for f in blah`) and see if you can deduce it.

There are also 2 list files provided `Bdor_GBS_list` for the full list of individuals and `Bdor_GBS_list_small` with just a handful of specimens for fast run times.

## Installing Stacks
Stacks requires a configure/make/make install, as stated in the [manual](https://catchenlab.life.illinois.edu/stacks/manual/#install). But thankfully it's a painless install with default compilers on MCC.
```
wget https://catchenlab.life.illinois.edu/stacks/source/stacks-2.65.tar.gz
./configure --prefix=/scratch/jdu282/stacks-2.65
make
make install
```

## Running Stacks
Running stacks is a relatively simple process, as they provide wrappers that will do either de novo assembly or reference-based in more-or-less 1 line of code. We've got a reference genome here, so reference-based it is. In general, there are three steps in the reference-based process:
1. demutiplexing
2. mapping
3. calling SNPs
However, for this dataset, what was available on NCBI wasn't raw sequence data, but demultiplexed, filtered data per individual. So we don't need to do step #1 here (although I've included an example of that step in the script preceeded with ##!##). There are several separate scripts that make up a full Stacks run (`ustacks`, `cstacks`, `sstacks`, `gstacks`, etc.), and these wrappers just automate the process of calling each individual step. You could run them all individually, and in some cases (complicated datasets, or if you want to call SNPs in one dataset based on the catalog generated from another dataset) you need to understand what each step is doing. But for today, we'll use the wrapper to make this process quick and simple (relatively speaking...).

Here's what my job script looks like:

```
#!/bin/bash
#SBATCH --time 24:00:00     # Time requested to run the job (format hours:minutes:seconds or days-hours:minutes)
#SBATCH --job-name=stacks    # Job name
#SBATCH --nodes=1        # Number of nodes to allocate. Same as SBATCH -N
#SBATCH --ntasks=24       # Number of cores to allocate. Same as SBATCH -n
#SBATCH --account=coa_jdu282_brazil_bootcamp2023  # Project allocation account name; use coa_jdu282_brazil_bootcamp2023
#SBATCH --partition=normal
#SBATCH --mail-type ALL    # Send email when job starts/ends/fails; other value option: ALL, NONE, BEGIN, END, FAIL, REQUEUE
#SBATCH --mail-user julian.dupuis@uky.edu   # your email to receive notifications about submitted job 
#SBATCH --mem=32g
#SBATCH -e %x.%j.err
#SBATCH -o %x.%j.out

##!## 
##!## step 1. demultiplexing (this is example code from a different dataset, but is data from a novaseq lane, which the sequencer was able to demultiplex into 3 subpools)
##!## mkdir 1_demultiplexed_noRadcheck
##!## process_radtags -1 ./0_RawData/Hemi_subpool1_R1.fastq.gz -2 ./0_RawData/Hemi_subpool1_R2.fastq.gz -b ./0_RawData/Hemi_subpool1_barcodes --renz_1 mluCI --renz_2 nlaIII --inline_null -w 0.15 -s 10 -c -q -r --filter_illumina  -E phred33 -D -o ./1_demultiplexed_noRadcheck --disable-rad-check
##!## process_radtags -1 ./0_RawData/Hemi_subpool2_R1.fastq.gz -2 ./0_RawData/Hemi_subpool2_R2.fastq.gz -b ./0_RawData/Hemi_subpool2_barcodes --renz_1 mluCI --renz_2 nlaIII --inline_null -w 0.15 -s 10 -c -q -r --filter_illumina  -E phred33 -D -o ./1_demultiplexed_noRadcheck --disable-rad-check
##!## process_radtags -1 ./0_RawData/Hemi_subpool3_R1.fastq.gz -2 ./0_RawData/Hemi_subpool3_R2.fastq.gz -b ./0_RawData/Hemi_subpool3_barcodes --renz_1 mluCI --renz_2 nlaIII --inline_null -w 0.15 -s 10 -c -q -r --filter_illumina  -E phred33 -D -o ./1_demultiplexed_noRadcheck --disable-rad-check
##!## 

# step 2. mapping
bwa-mem2 index ./genome.fasta

mkdir 2_bwa2_mem
for f in `cat list`; do bwa-mem2 mem -t 24 genome.fasta RawData/fastq/$f.fastq.gz > 2_bwa2_mem/$f.sam; done

module load samtools-1.12-gcc-9.3.0-zo3utt7

for f in `cat list`; do samtools sort -O BAM -@ 24 2_bwa2_mem/$f.sam > 2_bwa2_mem/$f.bam; done

# step 3. ref_map wrapper
mkdir 3_stacks_out
ref_map.pl --samples 2_bwa2_mem -o ./3_stacks_out -T 24 --popmap ./list_popmap -X "populations: -p 1 -r 0.01 --write-random-snp --ordered-export --vcf"

#populations -P ./3_Stacks_out -t 24 -M ./list_popmap  -p 1 -r 0.01 --write-random-snp --ordered-export --vcf
```

As I like to do with GATK steps, you can see I number the output directories here. Mapping (step 2) is just as we've done it before with a for loop, so only requires your list of specimen IDs (matching the names of the read files) and the reads themselves. Stacks itself (step 3) makes use of an additional file, which is called a population map (or `popmap`). This is basically a fancier version of our list of specimens that also contains a population identifier for each individual. So it's a 2-column file like this:
```
carambolae_FrenchGuiana_ms04723	1
carambolae_FrenchGuiana_ms04724	1
dorsalis_Hawaii_ms05004	1
dorsalis_Hawaii_ms05015	1
dorsalis_Vietnam_ms05029	1
dorsalis_Vietnam_ms05030	1
```
where "1" indicates that all 6 of these individuals belong to the same population. This can be used to apply population-specific filters (see [here](https://catchenlab.life.illinois.edu/stacks/comp/populations.php)), and we'll talk about that as a group. In this case, I'm simplifying our lives and just treating the whole dataset like one big population. To easily create the `popmap` file, we can use our list of individuals and do this:
```
sed "s/[a-zA-Z]*$/\t1/" list > list_popmap
```

You'll see that at the end of that job script, there is a commented-out command to call `populations`. `populations` is the final stage of both the de-novo and reference-based Stacks pipelines. It basically takes the catalog that Stacks generates and generates final SNP files of various formats and with various filters. We will talk as a group about my philosophy for filtering SNPs in stacks versus after the fact, but I often end up running populations a few times, either to generate additional data formats for the SNP dataset, or to apply different filters to see what their effect is on the final dataset. So I included a line to run `populations` by itself here (this is one of the individual scripts, similar to `cstacks` or `gstacks` that I mentioned earlier). 

