# Downloading and manipulating a genome, and also exploring PATHs and bash_profile
Let's download a high quality reference genome for B. dorsalis, and then we'll do some basic editing of it. We'll download [this one](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023373825.1/).
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/373/825/GCF_023373825.1_ASM2337382v1/GCF_023373825.1_ASM2337382v1_genomic.fna.gz
```
One of the tools we can use for manipulating this genome file will require an uncompressed version, so let's do that too.
```
gunzip -c GCF_023373825.1_ASM2337382v1_genomic.fna.gz > GCF_023373825.1_ASM2337382v1_genomic.fna
```
Let's take a look at this genome. Can you predict what each of these commands will do?
```
ls -lh GCF_023373825.1_ASM2337382v1_genomic.fna
wc -l GCF_023373825.1_ASM2337382v1_genomic.fna
grep ">" GCF_023373825.1_ASM2337382v1_genomic.fna
grep ">" GCF_023373825.1_ASM2337382v1_genomic.fna | wc -l
zcat GCF_023373825.1_ASM2337382v1_genomic.fna.gz | grep ">" | wc -l
```
(the last one is a bit repetitive, but using `zcat` and piping is a nice way to access information that's in a compressed file)

### Downloading seqtk
Now, let's download [seqtk](https://github.com/lh3/seqtk), and read how to install it. Do this in your personal folder in pscratch on MCC. 
```
git clone https://github.com/lh3/seqtk.git
cd seqtk
ls -l
make
```
If you `ls -l` in this directory, you should see something like this:
```
[jdu282@mcc-login001 seqtk]$ ls -l
total 258
-rw-r--r-- 1 jdu282 users  18784 Nov 26 14:04 khash.h
-rw-r--r-- 1 jdu282 users   8634 Nov 26 14:04 kseq.h
-rw-r--r-- 1 jdu282 users   1086 Nov 26 14:04 LICENSE
-rw-r--r-- 1 jdu282 users    284 Nov 26 14:04 Makefile
-rw-r--r-- 1 jdu282 users    809 Nov 26 14:04 NEWS.md
-rw-r--r-- 1 jdu282 users   1780 Nov 26 14:04 README.md
-rwxr-xr-x 1 jdu282 users 309312 Nov 26 14:04 seqtk
-rw-r--r-- 1 jdu282 users  65772 Nov 26 14:04 seqtk.c
```
You can see there is one additional file in there that wasn't in there before the `make` command, and it's the name of the software and has high accessibility in the first column (`-rwxr-xr-x`). The next step is usually to test it out and see if it works (i.e. if it's installed correctly). To make sure any software is installed correctly, my first course of action is to try to execute the software and see whether it produces something sensible or some random error. My go-to commands are the following derivations of these, here using seqtk as the example command:
```
./seqtk        # just trying it out without any options
./seqtk --help        # usually help will bring up all the options if something is installed correctly
./seqtk -h      # some softwares allow both abbreviated and long versions of options, some only one or the other
```
Which one of these will work will depend on the software in question; as you can see, seqtk doesn't like either help command, but does like being called without options. Even if you've never used the software you're installing you should be able to tell whether the output of those commands is something that makes sense (i.e. we know this software deals with sequence/genome manipulation and you see options specific to that kind of thing, as compared to some error message that has nothing to do with sequencing). 

What if you try the above command without `./`? You'll probably see something like this: 
```
[jdu282@mcc-login001 seqtk]$ seqtk
-bash: seqtk: command not found
```
Even though you're in the directory where you just installed seqtk, your cluster account does not automatically know to look in there for the executable file, which is why we needed the `./` above. Likewise, you'll get the same error if you do this
```
cd ../
./seqtk
```
Since we moved up one directory, seqtk is no longer in the current directory as we are specifying, so the cluster returns an error. 

There's a few options to solve this issue, and increase the flexibility as to where you can install software on the cluster. First, you could use the full path everytime you call any software. So we just need to get the `pwd` of the location we downloaded seqtk (Here I downloaded a backup version for everyone to use):
```
[jdu282@mcc-login001 seqtk]$ pwd
/pscratch/jdu282_brazil_bootcamp2023/programs/seqtk
```
And we can now cd to anywhere, and be able to call seqtk using `/pscratch/jdu282_brazil_bootcamp2023/programs/seqtk/seqtk`. Give it a try by cd'ing somewhere else and calling that command. 

### Knowing your PATH
The other option is to add seqtk to our PATH. $PATH is an environment that the cluster system (or your own personal computer) looks in to find executables and software. We can see what's in our PATH with `echo $PATH`. The output to that might look complicated, but it's just a long list of directory locations, separated by `:` showing all the locations your system will look for a command every time you type that command. We can add other locations to PATH with the following:
```
export PATH=$PATH:[path to where you just downloaded seqtk]/
```
This adds that location to your PATH for this session only, and so now you can cd anywhere and access seqtk by just typing `seqtk`. Give it a try. A couple things to NOTE:<br>
1. This export only works for your current session. So if you log off on the cluster and log back in, your PATH will no longer contain this location (see below, bash_profile for permanent solution)
2. You need to export the directory location of where the executable command is, not the actual executable command. i.e. for my backup install above, this export command would work `export PATH=$PATH:/pscratch/jdu282_brazil_bootcamp2023/programs/seqtk` while this one would not work `export PATH=$PATH:/pscratch/jdu282_brazil_bootcamp2023/programs/seqtk/seqtk`. This is a common problem when installing and trying out new software. You need to add the location of the executable to your PATH, not the executable itself. 

### Bash_profile
The permanent solution to NOTE 1 above is to add that export path command to a file in your home directory. Problem is, that file is normally invisible to keep you from accidentally deleting it, and in the case of MCC, you have to generate it manually. <br>
1. Go to your home directory `cd ~` (using tilda `~` is a nice shortcut for your home directory, so is `$HOME`; both are shorter than doing `cd /home/[username]`.
2. Create your bash_profile file: `nano bash_profile`
3. Paste in the following:
```
# .bash_profile

# User specific environment and startup programs
	
# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

PATH=$PATH:$HOME/.local/bin:$HOME/bin

export PATH
```
Now below that, you can paste in the `export PATH` command that you used before. So it should look something like this:
```
# .bash_profile

# User specific environment and startup programs
	
# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

PATH=$PATH:$HOME/.local/bin:$HOME/bin

export PATH

export PATH=$PATH:/pscratch/jdu282_brazil_bootcamp2023/programs/seqtk
```
Now everytime you log into this cluster, the system looks in this file and adds all those locations with `export PATH` lines to your PATH. This essentially creates a permanent solution for adding locations to PATH. NOTE, when you edit this file, the change doesn't go into effect immediately. As I said, the system looks here everytime you login. So, once you edit this file, you need to log out and log back in to make the change effective. 

### Using seqtk
Ok, now that you've got it installed and in your PATH, let's use this tool. You should have seen from above (hint, `grep ">" GCF_023373825.1_ASM2337382v1_genomic.fna | wc -l`) that this genome file has 587 scaffolds, but if you omit the `wc -l` from that command, you can see there are 6 named chromosomes, the mitochondrial genome, and 580 random scaffolds. Knowing that this species has 6 autosomes gives us a pretty good hint that these 580 random scaffolds are probably lower quality or poorly mapping/assembling regions. Let's test that out by calculating the size of each of these entries in this fasta file. We can do that with samtools, and MCC has a module for samtools. 
```
module load samtools-1.12-gcc-9.3.0-zo3utt7
samtools faidx GCF_023373825.1_ASM2337382v1_genomic.fna
ls -l       # what's new in here?
nano GCF_023373825.1_ASM2337382v1_genomic.fna.fai
```
This is called an index file, and is often a requirement for mapping reads or manipulating genomes and associated data (we'll use these frequently for this class). You can see what the column headers mean [here](http://www.htslib.org/doc/faidx.html), but the main thing for right now is that column 1 is the scaffold name and column 2 is how long that scaffold is. So you can see those first 6 scaffolds (the named chromosomes) are all much longer than the rest of them, and the mitogenome is ~16k bp, which is standard for most insects. 

This is a pretty good genome, and for that reason, we might want to just limit ourselves to the named chromosomes and ignore the rest of the assembly since it's not going to give us any location information (gene proximity). Seqtk can help us do that. Unfortunately, the documentation for seqtk is not great, but what we can use is the `subseq` command, which will extract sequences from a multi-fasta that are found in a file we provide. So from the output above (listing the chromosomes in the assembly file), if you create a file named "Bdor_chromosomes" with the following contents:
```
NC_064303.1 Bactrocera dorsalis isolate Fly_Bdor chromosome 1, ASM2337382v1, whole genome shotgun sequence
NC_064304.1 Bactrocera dorsalis isolate Fly_Bdor chromosome 2, ASM2337382v1, whole genome shotgun sequence
NC_064305.1 Bactrocera dorsalis isolate Fly_Bdor chromosome 3, ASM2337382v1, whole genome shotgun sequence
NC_064306.1 Bactrocera dorsalis isolate Fly_Bdor chromosome 4, ASM2337382v1, whole genome shotgun sequence
NC_064307.1 Bactrocera dorsalis isolate Fly_Bdor chromosome 5, ASM2337382v1, whole genome shotgun sequence
NC_064308.1 Bactrocera dorsalis isolate Fly_Bdor chromosome 6, ASM2337382v1, whole genome shotgun sequence
```
We can then extract just those sequences and write them to a new file with:
```
seqtk subseq GCF_023373825.1_ASM2337382v1_genomic.fna Bdor_chromosomes > GCF_023373825.1_ASM2337382v1_chromosomes.fa
```
Actually, most multi-fasta file interpreters will only look at the first part of each scaffold name, so we could even have just this in our file, and it would work the same:
```
NC_064303.1
NC_064304.1
NC_064305.1
NC_064306.1
NC_064307.1
NC_064308.1
```
Actually, we don't even need seqtk for this task (although it can do other processes which make it useful). We could have used samtools as well: `samtools faidx GCF_023373825.1_ASM2337382v1_genomic.fna NC_064303.1 NC_064304.1 NC_064305.1 NC_064306.1 NC_064307.1 NC_064308.1` Give these a try to make sure I'm not lying to you!

Finally, how about we simplify those chromosome names, so they are sensible ("Chromosome1") instead of the NCBI reference number ("NC_064303.1"). As a final test, can you look at the poor help documentation of seqtk and figure out what this command is doing step by step?
```
seqtk rename GCF_023373825.1_ASM2337382v1_chromosomes.fa Chromosome > temp && awk '{print $1}' temp > GCF_023373825.1_ASM2337382v1_chromosomes_renamed.fa && rm temp
```
