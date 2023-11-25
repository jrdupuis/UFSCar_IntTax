## Analysing a GBS dataset with Stacks

Let's take a look at this [bioproject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA987164) from [this paper](https://www.sciencedirect.com/science/article/abs/pii/S1055790323001926)
What about the metadata for this dataset? If you click on the SRA link (link to 798 SRA files), and then click on "send to" > "file" > and "run info", that opens up the information associated with each specimen through their BioSample and SRA submission. But not all pertinent info is there (country of origin, collection locality, etc.). 
But if you go to the publication and look at Table S1, the authors included SRR numbers and pertinent locality information, so we're in luck.

## Installing Stacks
Stacks requires a configure/make/make install, as stated in the [manual](https://catchenlab.life.illinois.edu/stacks/manual/#install).
```
wget https://catchenlab.life.illinois.edu/stacks/source/stacks-2.65.tar.gz
./configure --prefix=/scratch/jdu282/stacks-2.65
make
make install
```

# make pop map from list
sed "s/[a-zA-Z]*$/\t1/" list > list_popmap
