## Structure (Pritchard et al. 2000)
STRUCTURE is a very commonly used program in population genetics to assess population clustering/structure. See more info [here](https://web.stanford.edu/group/pritchardlab/structure.html). 

Although premade binaries are available for structure, they are compiled for a 32 bit system, and MCC is 64 bit, so we need to make from scratch. I've installed structure here `/pscratch/jdu282_brazil_bootcamp2023/programs/structure_kernel_src/` and the executable is called `structure`, but you can also download it from [here](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html) (wget the "c code" option under "Download source code", and `make` within that directory). Structure needs 3 input files, a `.str` file containing genotype data, and two parameter files: `mainparams` and `extraparams`. The software download contains examples of these files. I've included example files for the swallowtail butterfly microsatellite [dataset](https://onlinelibrary.wiley.com/doi/full/10.1111/jeb.12931) in this directory. We will walk through the contents of these files together. 

## Running structure
As discussed, structure runs a Bayesian MCMC process where it evaluates different population clustering situations with the individuals in the dataset. It replicates this hundreds of thousands to millions of times within the MCMC process, replicates those MCMC processes usually tens of times, and does all of this independently for each value of K (the number of unique genetic clusters). So we actually need to run tens to potentially hundreds of processes, depending on the values of K evaluated and the number of replicates of the MCMC process. The way I do this is to submit an individual job for each replicate for a given value of K, so that I'm not waiting for a for loop to go through each replicate individually before starting on the next. 

`structure_qsub_header` contains a pretty standard job submission header (note where I'm writing stdout and stderr, and that I'm not including an email address in this job header):
```
#!/bin/bash
#SBATCH --time 24:00:00     # Time requested to run the job (format hours:minutes:seconds or days-hours:minutes)
#SBATCH --job-name=structure  # Job name
#SBATCH --nodes=1        # Number of nodes to allocate. Same as SBATCH -N
#SBATCH --ntasks=1       # Number of cores to allocate. Same as SBATCH -n
#SBATCH --account=coa_jdu282_brazil_bootcamp2023  # Project allocation account name; use coa_jdu282_brazil_bootcamp2023
#SBATCH --partition=normal
#SBATCH --mem=1g
#SBATCH -e ./err/%x.%j.err  # Error file for this job.
#SBATCH -o ./out/%x.%j.out  # Output file for this job.

export PATH=$PATH:/pscratch/jdu282_brazil_bootcamp2023/programs/structure_kernel_src/
```

`structure.sh` contains the following:
```
qsub_folder=qsubs
output_folder=StructureResults

mkdir $qsub_folder
mkdir $output_folder
mkdir err
mkdir out

for k in {1..4}; 
do 
for r in {1..5}; 
do
echo k_$k.rep_$r 
cat structure_qsub_header > $qsub_folder/k_$k.rep_$r.qsub
echo "structure -D $RANDOM -K $k -o $output_folder/$k.$r.output" >> $qsub_folder/k_$k.rep_$r.qsub
sleep 1
sbatch $qsub_folder/k_$k.rep_$r.qsub &
done
sleep 1
done
```
You can see that it sets up a directory structure for the job submission files (`$qsub_folder`), the output and error dirs (`out` and `err`), and then is a double for loop to generate jobs for a set of K values (`for k in `), and a set of replicates per K (`for r in `). We'll talk about the use of the `-D` option and the `sleep` commands, and why those are important.  

This command (and these associated file) can be submitted with `bash structure.sh`. Make sure all required files (`mainparams`, `extraparams`) are in the same directory as `structure.sh` and `structure_qsub_header`. This set of jobs should take <10 minutes to run, and we'll walk through the outputs. We'll make use of [StructureSelector](https://lmme.ac.cn/StructureSelector/index.html) for processing/summarizing the run and to get access to the output file `ClumppIndFile.output` for the best value of K for plotting in excel.
