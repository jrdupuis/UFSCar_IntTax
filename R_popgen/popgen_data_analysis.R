########
########
######## Population genetic data analysis in R
######## Session 1: basics of R and adegenet
########
########


# In general, functions in R follow this format:
function()
# where "function" is the name of whatever function/operation you're doing and within the parentheses you spell out all of the parameters, inputs/outputs, and options for the given function.


# if you're not sure what options a function has, you can do the following:
??function
# This will open up a web browser and get you to an info page about functions (you will have to click on the individual function in question in the browser to get to the info page; if there are multiple functions with the same name (i.e. from different libraries), you'll have to pick the right one from the list that pops up in your browser).


# Install a few packages that will be useful for population genetic analyses. Packages (or sometimes called libraries) contain specific functions that are useful for specific tasks. These ones that we're loading are all made to do various population genetic analyses/calculations/statistics. You can often find out more about packages and their functions by looking at their CRAN or github websites. For example, for adegenet, see https://cran.r-project.org/web/packages/adegenet/index.html and https://github.com/thibautjombart/adegenet (you can find manuals, tutorials, etc. at those websites).
install.packages("adegenet")
install.packages("poppr")
install.packages("pegas")
install.packages("hierfstat")

# Once you have a package installed on your computer (the "install.packages()" function), then you will just need to load the library with library() every time you start up R. i.e. once you've ran the above scripts once on your computer, you'll never have to do those again (unless you need to update the package). 
library(adegenet)
library(poppr)
library(pegas)
library(hierfstat)

# adegenet contains some nice toy datasets that we can use as examples.
# load one of adegenet's toy datasets called nancycats, which is a bunch of microsatellite data for cat populations from Nancy, France (hence the fun name)
data("nancycats")

# once you have run that data() function, the nancycats dataset is loaded into memory, and you can access it by simply entering
nancycats
# There you should see the following output:
					# /// GENIND OBJECT /////////
					
					 # // 237 individuals; 9 loci; 108 alleles; size: 150.5 Kb
					
					 # // Basic content
					   # @tab:  237 x 108 matrix of allele counts
					   # @loc.n.all: number of alleles per locus (range: 8-18)
					   # @loc.fac: locus factor for the 108 columns of @tab
					   # @all.names: list of allele names for each locus
					   # @ploidy: ploidy of each individual  (range: 2-2)
					   # @type:  codom
					   # @call: genind(tab = truenames(nancycats)$tab, pop = truenames(nancycats)$pop)
					
					 # // Optional content
					   # @pop: population of each individual (group size range: 9-23)
					   # @other: a list containing: xy 
# This tells you that nancycats is what's called a GENIND object in R (i.e. population GENetic data for INDividuals), and it gives you some basic info (number of individuals, loci, alleles, etc.). Within this object you have both allele/genotype data and some other info. You can access different sets of info stored within this object by referring to nancycats with the "@" symbol and whatever other info you want to see. For example, if you type out:
nancycats@ploidy
# that will return a list of the ploidy of all individuals in the dataset (there are 237 individuals, and they're all diploid, so relatively boring).

# The @pop subobject contains a list of the populations that each individual belongs to. So the following will give you a list of 237 population assignments with rather arbitrary population names ("P01", "P02", etc.)
nancycats@pop

# The @tab subobject contains counts of all alleles in all individuals, so if you entered "nancycats@tab" that would return a 237 row x 108 column block of text. We can access just the first 5 lines of that by also using the head() function:
head(nancycats@tab)
# From that output, you can see that that the rows are individuals (with individual identifier numbers that start with "N") and columns are locus.allele combinations. So "fca8" is a locus name, and it is a dinucleotide repeat microsatellite, with alleles of size 117 bp, 119 bp, 121 bp, etc. all the way to 149 bp. After locus fca8 is locus fca23, and on and on. 

# Another way to subset that big block of data from @tab is to use square brackets to pull out specific rows and/or columns by their number (e.g. row 1, row 2, column 1). You can also combine that square brackete subsetting with calling individual loci from the genind object. 
nancycats[10:18, loc="fca8"]@tab  
# 10 refers to number of rows being displayed, 18 to the number of columns displayed, and then the only locus being called is fca8 from the @tab subobject (note, the row and column headers are included in those 10 and 18 numbers).
# These are allele counts, so if you look at the row for individual N224, you can see it has 0s for all entries except a 2 in column fca8.135, so this individual has two copies of the 135 bp allele for this locus (it's a homozygote). Individual N7 is a heterozygote for alleles 137 and 141. 


## The native R function summary() is a useful function for lots of data handling/analysis in R. If we use it on nancycats, we can see it automatically knows that it's dealing with population genetic data (because it's a genind object) and it'll actually calculate a few basic statistics for us:
summary(nancycats)
# lots of that info should be familiar to you at this point.

# Another useful thing is to take the output of an R function and assign it to a new object name. This will allow us to identify individual parts of, for instance, that summary() function. So the following is taking the same summary() function, but assigning it to a new object called "summ"
summ <- summary(nancycats)
# As you can see, that last command doesn't return any output from R in the R console, but if we call that new name, we'll get the same output as before:
summ
# From this point, you can see that summary contains several different statistics, such as number of individuals, group (population) sizes, heterozygosities, etc. Often when R outputs this genearlized stats, it also assigns shorter names to those sets of numbers so we can access them easily. To tell what the shorter names are we can run:
names(summ)
# And you can use those names to call individual sets of numbers, for instance the observed heterozygosity:
summ$Hobs
# Those will show more decimal points than before (the summary() function simplifies the numbers for us), but it's the same. So what this allows us to do is access each of those sets of numbers individually. What if we wanted to compare the observed and expected heterozygosities? Well we could do something like this:
plot(summ$Hobs,summ$Hexp, xlab="Hobs", ylab="Hexp", main="Expected heterozygosity as a function of observed heterozygosity per locus")
	# NOTE: plotting and graphing stuff in R is an addicting art of getting the right code to do what you want. We aren't going to focus too much on the all the intricacies of plotting here, but google is your best friend in making really sharp looking graphs in R.


# If you dig around in the manual/tutorials for adegenet, you'll see that in addition to the "@" subobjects within genind objects, there are additionall accessors within that object to access specific info about the object's data:
				## Additional accessors:
				# nInd: returns the number of individuals in the object; only for genind.
				# nLoc: returns the number of loci.
				# nAll: returns the number of alleles for each locus.
				# nPop: returns the number of populations.
				# tab: returns a table of allele numbers, or frequencies (if requested), with optional replacement of missing values; replaces the former accessor ’truenames’.
				# indNames: returns/sets labels for individuals; only for genind.
				# locNames: returns/sets labels for loci.
				# alleles: returns/sets alleles.
				# ploidy: returns/sets ploidy of the individuals; when setting values, a single value can be provided, in which case constant ploidy is assumed.
				# pop: returns/sets a factor grouping individuals; only for genind.
				# strata: returns/sets data defining strata of individuals; only for genind.
				# hier: returns/sets hierarchical groups of individuals; only for genind.
				# other: returns/sets misc information stored as a list.
# Some of these are useful for exploring new datasets. i.e. if I hadn't told you that individual identifiers started with "N" and loci with "fca", differentiating the loci from the individuals can be difficult. You can check with:				
indNames(nancycats)
head(indNames(nancycats), 10)  # note, we can also add an additional option to the head command to show a specific number of entries instead of 5.
locNames(nancycats)

# Another useful way to look at individual's alleles/genotypes is to take the data from the genind object and turn it into a simple data frame ("df" in R speak) with the following function:
??genind2df
# So we can take nancycats, turn it into a df and separate the alleles for each individual with a "/". This will be more similar to how we'll look at genotype data later on.
tail(genind2df(nancycats, sep = "/"))


# Let's look at a few functions within the library poppr that are useful for some additional summary statistics. By using "??" you can see these have more options than some of the previous functions (especially info_table() ). 
??locus_table
??info_table
# Let's run these two on nancycats. The first is useful for characteristics of the loci.
locus_table(nancycats)
# And the latter can look at primarily missing data or at ploidy level (which again is boring in this dataset). So let's look at missingness by setting the "type" option for info_table()
info_table(nancycats, type = "missing")
# There you can see the missingness per population and locus. Some of those combinations have pretty high missingness, which can be problematic. You'll have noticed that many of the options for this function are for plotting, so let's re-do that command and turn the plot funtion on
info_table(nancycats, type = "missing", plot=T)
# That's nice! colors!


# Imagine trying to calculate FST or something for locus fca45. Population P17 has 100% missing data, so it's going to be pretty tough (i.e. impossible) to get its observed/expected heterozygosities. With this type of dataset, we'd often remove loci with missing data higher than a certain threshold. So let's do that with something like the following
nancycats %>% missingno("loci", cutoff = 0.05) %>% info_table(plot = TRUE, scale = FALSE)
# This has a few things going on. First, the "%>%" is what's known as a pipe in coding speak, and it takes the output of one function and passes it to another function. So in this case we're taking nancycats and instead of referring to it as an option within the parentheses of funtion(), we're piping it to the function missingno(). That function is removing any loci with average missingness higher than a 5% (0.05) cutoff, and then it creates a new genind object with those loci removed. Finally, we're then piping that output (the new genind object) to info_table() again to create a new missingness plot. We could also just create the new genind object and assign it to a new object name, and double check that the new object makes sense
nancycats_noMiss <- nancycats %>% missingno("loci", cutoff = 0.05)
nancycats_noMiss

# We can also use missingno() in the same fashoin to remove genotypes with high missing data, say anything higher than 25%:
nancycats %>%  missingno("geno", cutoff = 0.25) %>%  info_table(plot = TRUE)


# The library poppr also has a function called poppr(), which is pretty useful for general stats.
poppr(nancycats)
# what's all that mean? Let's check the help page and select the function we're using (hint it's the poppr::poppr link)
??poppr
# We don't need to get into the nitty gritty details of all this, but you can see a description of the various columns under "Value" on the help page. Some of this is similar to the previous summary functions we used, but some is new. For example, the Shannon-Weiner Diversity Index is a useful general measure of "genetic diversity", which is a bit more sophisticated than plain heterozygosity. It might be interesting to compare those two sets of values (i.e. are heterozygosity and the Shannon-Weiner index correlated?). We can figure out the shorthand names for these datasets similarly to above using names(), and do this comparison like so:
names(poppr(nancycats)) 
# more conveniently than the last time we did this, the column names are already the shorthand codes for calling individual column's data 
plot(poppr(nancycats)$H, poppr(nancycats)$Hexp)
# If you look back at the original output of poppr(), you'll see that the Shannon-Weiner value that's at 5.5 is actually from a "total" calculation in the last row, so that's not really what we are looking for, so let's take those two sets of numbers, and remove the last entry with the head() function. The "-1" in this function is basically telling head to start from the end of the list of numbers and show all but the last one. To make sure it works, let's write both the original list and the edited list to objects, then we can return both of those to compare their contents:
original_list <- poppr(nancycats)$H
no_last_entry <- head(poppr(nancycats)$H, -1)
original_list
no_last_entry
# you could also use the length() function to count the number of entries in these objects
length(original_list)
length(no_last_entry)
# Now we can re-do the plot with those head() functions added in to get rid of those total calculations
plot(head(poppr(nancycats)$H,-1) , head(poppr(nancycats)$Hexp, -1))

# what about population size and genetic diversity? This time, let's write the pertinent data into new objects and call those directly with plot to make the plot command a bit easier to read
Hdiv <- head(poppr(nancycats)$H, -1)
Hdiv
length(Hdiv)
N <- head(poppr(nancycats)$N,-1)
length(N)
plot(N,Hdiv)
# There's a nice correlation! The more individuals you sample, the more alleles you're likely to sample, so you find more genetic diversity.


#####
# To practice what you've just learned, let's work through same kind of commands with a new dataset. Create new commands below for this dataset, and try to summarize it in the following general ways: What is this dataset and how big is it (# individuals, loci, etc.), what is the genetic diversity of different populations/loci, and does the data have a lot of missingness? The dataset to look at is another example dataset in adegenet, this time with cows! 
data(microbov)
microbov





########
########
######## Population genetic data analysis in R
######## Session 2: more serious popgen
########
########

# Alright, we're going to be using many of the same packages as last week, so let's load them up:
library(adegenet)
library(poppr)
library(pegas)
library(hierfstat)

# HWE! Let's think about one of the big assumptions in population genetics, which is Hardy-Weinberg Equilibrium. The library pegas gives us a nice easy way to calculate HWE with hw.test(). Let's go back to our feline friend nancycats and see what we can do with that function
data(nancycats)
??hw.test
# Not very many options, but you'll see the main one is B which is a Monte Carlo procedure. The math behind that procedure is a bit more than we need to deal with here, but basically what that does is give you another way to assess the statistical significance of the HWE Chi-Squared GOF test, like we went over in lecture, but in this case it basically permutes the dataset to create an expected distribution given the alleles in said data. Let's ignore it for the meantime and set B=0
hw.test(nancycats, B = 0)
# In that output, you should see a few familiar things. We can see it's testing HWE for each locus (makes sense, as we covered that in lecture), then it gives us a chi^2 value (this is the calculated value), degrees of freedom (df), and a Pr(chi^2) which is basically the P-value. Compared to the examples we did in lecture, both the chi-squared and df values should look really big. And if you look up a chi-squared critical table online, you'll notice that even when the degrees of freedom are really high (like above 100), the critical value is way lower than the calculated values we're seeing here. E.g., for an alpha of 0.05 and df of 100, the critical value is only 77 and almost all of our calculated values are way above that (for a more moderate example, when df=55, the critical value is 34). 
# But if we think about what we just did, maybe we can figure out why those numbers are so ridiculously high. What exactly is in nancycats, as far as the number of loci and number of populations? Can you remember how to do that from last week?

# So if you remember something like nancycats@pop, you'll remember that there are 17 populations in this dataset. And we just calculated HWE for the entire dataset, not just a single population. So that explains why both the df and calculated chi-squared critical values are so high. We already saw that the genetic diversity of the populations in nancycats differs, so we can assume there are different allele frequencies in those 17 populations. So what we need to do is to separate the populations of nancycats and test for HWE for each of them. Let's start off with just testing for HWE in one population using the seppop() function from adegenet. 
??seppop
hw.test(seppop(nancycats)$P01, B=0)
# So here the df and chi-squared calculated values are a bit closer in magnitude to what we saw in lecture. And if you check a chi-squared critical table, you'll find critical values more in-line with these calculated values. E.g. for alpha=0.05, df=6, then the critical value=1.635 and for df=10, crit=3.9. But let's go back to our lectures on HWE and think about these values. Following the method we used in lecture, what would be your statistical conclusion for for the first couple loci (fca8 and fca23) given those chi^2 and df values? Do those conclusions match up with the P-values calculated by hw.test()? 
# So what's going on? Let's chat as a group when you get to this point.

				######## JRD CHEAT SHEET
				# For fca8, our lecture approach would say calc > crit (5.13 > 1.635), so we'd reject the null of HWE (fca8 is not in HWE), but it's P-value indicates that the test statistic is not significant (therefore, the locus is in HWE)
				# For fca23, our lecture approach would also say calc > crit (26.04 > 3.9), so we'd reject the null of HWE (fca23 is not in HWE). Here, the P-value comes to the same conclusion
				# STATS stuff: exact test vs. approximate test (different ways to calculate the expected distribution with which we are comparing our test stat)
				# change B to 1000 and chat about that


# Let's ignore the math that goes into that, and just trust the P-values that hw.test() provides us. 
# If we wanted to go ahead and test for HWE for all loci, we could use a pipe command like we used last week. i.e. we could use seppop() without specifying a population name and then pipe the output to hw.test():
seppop(nancycats) %>% lapply(hw.test, B = 0)
# This is not super digestible, so we can use some other R functions to make a nice plot to show us what locus/population combinations have significant P-values. We don't need to go into too much detail about these plotting steps, but see if you can understand what they are doing
(nanhwe.pop <- seppop(nancycats) %>% lapply(hw.test, B = 0))			 # re-do the hw.test() and write it to a new object
(nanhwe.mat <- sapply(nanhwe.pop, "[", i = TRUE, j = 3))			 # Take the third column with all rows
alpha  <- 0.05			 # sets a variable equal to our alpha
newmat <- nanhwe.mat			 # renaming the matrix we created so we can replace anything with a P-value>0.05 below
newmat[newmat > alpha] <- 1			 # replacing the P-values of anything more than our alpha (which is 0.05)
library("lattice")			 # required for some plotting, which was actually installed as a dependency with poppr
levelplot(t(newmat), ylab="locus", xlab="population")				 # and make a plot
# In this plot, we've re-written all non-significant P-values as a value of 1, and are only retaining the original P-values for those <0.05. 

# Let's calculate FST. There are a lot of packages to do this, so let's explore a few. First, we can use a simple function wc() in the library hierfstat.
??wc
wc(nancycats)
# From the ??wc manual page, you can see under Value that the output of this function also creates some additional data that isn't shown by the default output. 
wc(nancycats)$per.loc

# The library pegas can also calculate FST (as well as the other F-statistics) with the Fst() function.
??Fst
Fst(nancycats)
# Well that didn't work, and the error message isn't terribly informative. If you look at the ??Fst manual page, you might notice that it requires an object of class "loci", and we're trying to give it an object of class "genind", so that's probably our issue. But pegas gives us a way around that with the following command:
??as.loci
Fst(as.loci(nancycats))

# From both of those manual pages, you can see both hierfstat and pegas are calculating Weir and Cockerham's FST, so theoretically those outputs should be the same. Let's see if we can isolate the FST values from the rest of the output from those two functions and compare them in a scatterplot like we did with heterozygosity last week. Try this out on your own and when everyone's to this point we can work through it together on the screen.
							###### JRD CHEAT SHEET for isolating Fst values on their own and plotting
							Fst1 <- wc(nancycats)
							Fst1$per.loc
							names(Fst1$per.loc)
							Fst1$per.loc$FST
							Fst2 <- Fst(as.loci(nancycats))
							Fst2
							names(Fst2)
							Fst2[,2]
							plot(Fst1$per.loc$FST,Fst2[,2])
					

# Let's test out another library, called mmod and its function Fst_Hedrick()
# install.packages("mmod")
library(mmod)
Gst_Hedrick(nancycats)
# And let's compare that to one of the sets of per-locus FSTs from the last step
Gst_Hedrick(nancycats)$per.locus
Fst1 <- wc(nancycats)
plot(Fst1$per.loc$FST,Gst_Hedrick(nancycats)$per.locus)
# Well, that's not quite as nice as before. So what's going on? Let's chat as a group about this. 
				#### JRD CHEATSHEET
				nc.diff_stats <- diff_stats(nancycats, phi_st=T)
				diff_stats(nancycats)$per.locus[,2]
				plot(Fst1$per.loc$FST, diff_stats(nancycats)$per.locus[,3])
				plot(diff_stats(nancycats)$per.locus[,4], diff_stats(nancycats)$per.locus[,3])
				with(nc.diff_stats, pairs(per.locus[,3:6], upper.panel=panel.smooth)) # some estimates correlate better than others


# What about a pairwise calculation of FST? For this, let's look at the genet.dist() function within hierfstat.
??genet.dist
# Pick a calculation and specify the method below:
genet.dist(nancycats, method = )
# Again, lots of libraries will calculate pairwise genetic distance/differentiation, so let's compare this with another calculation in adegenet. To do so, we'll have to tweak nancycats a little bit. Right now it's a genind object, so it storing the popgen information for individuals. But we want to calculate pairwise distance/FST between populations, so we need to turn the genind object into a genpop object, which stores info for populations rather than individuals. Let's check out genind2genpop()
??genind2genpop
# You can see the only really important option here is "pop", which tells the function what individuals belong to which populations. Thankfully, because we're dealing with a well-put-together genind object, we already have that in the @pop accession. So we can check that out, but in reality it's not really necessary as by defaul the genind2genpop() function will look in @pop for the population delimiters. But we can test that out by trying the funtion with and without specifying the pop option.
nancycats@pop
genind2genpop(nancycats, pop=nancycats@pop)
genind2genpop(nancycats) 

# We'll use the dist.genpop() function in adegenet for this calculation.
??dist.genpop
# And let's write that genpop object to a new name so we can refer to it in the dist.genpop() function.
nancyPop <- genind2genpop(nancycats, pop=nancycats@pop)
dist.genpop(nancyPop,method=1)

# Let's compare those two matrices of genetic distances, by first writing them to new objects and then plotting them.
genetDistAdegenet <- dist.genpop(nancyPop,method=1)
genetDistHierfstat <- genet.dist(nancycats, method = "Fst")
plot(genetDistAdegenet, genetDistHierfstat, xlab="adegenet", ylab="hierfstat")
# Pretty nice relationship despite the variation. Wonder if it's a significant correlation?
model <- lm(genetDistHierfstat ~ genetDistAdegenet) 	# Note, the linear model function (??lm and look for stats::lm) expects the dependent variable first (i.e. y ~ x instead of x ~ y like when we plotted the data), so pay attention to the order within those functions 
summary(model)
# Look for the p-value in there to see it's a significant relationship, and the R-squared value shows that the correlation is relatively strong. We can also add that linear model to the plot and the R squared value
abline(model, col="red", lty=2)  
legend("topleft",legend=paste("R-squared is", format(summary(model)$r.squared,digits=3)))


# One nice thing about nancycats, is that we have some spatial data to go along with the genetic data. See the @other accession. 
nancycats@other
# Not quite latitude and longitude, but if you looked into the metadata for nancycats, you'd find that these values are geographic locations of each population. So let's see if there's a relationship between geographic distance and genetic distance, i.e. isolation by distance. Remember that in IBD plots, we're looking at pairwise comparisons of the geographic distance vs. genetic distance for each possible pair of populations, so what we need is basically two matrices, one for genetic distance and one for geographic distance. We just made a couple matrixces for genetic distances, so now for geographic distance, we could calculate this using several packages. Since this isn't actual lat/long data, let's just use the dist() function in the stats package, which calculates a raw distance value from x, y coordinates
??dist
geogDist <- dist(nancycats$other$xy)
geogDist
# The important thing for this kind of calculation is making sure the matrices are formatted the same, i.e. same column/row headers, whether or not the diagonal is displayed, etc. Looks good in this case, so now let's look at the mantel.randtest() function
??mantel.randtest
ibd <- mantel.randtest(genetDistAdegenet,geogDist)
ibd
plot(geogDist, genetDistAdegenet)
dist_lm <- lm(as.vector(genetDistAdegenet) ~ as.vector(geogDist))
abline(dist_lm, col="red", lty=2)
# Let's chat about this as a group.


#### Let's change directions and go back to the cow dataset you played around with last week. One nice thing about this dataset is we have different ways to break up the individuals into different "populations". If we check out the @other tab we can see what we have to play with.
data(microbov)
microbov
microbov@other
# coun, breed, and spe refer to country of origin, breed (obviously), and species, respectively. So for all the calculations we've done which are per-population, we could treat any of these divisions as the "population". Adegenet provides an easy way to deal with these different "population levels", which it refers to as strata. The strata() function has some useful info about strata
??strata
# To define strata and switch strata, there are a few steps required, but once you get the hang of them they are pretty intuitive. Before we get into that, let's check what the @pop accession is for microbov by default. And also see if anything is stored as strata within the microbov object
microbov@pop
popNames(microbov) # two different ways to check the @pop accession
strata(microbov)
# So the population right now is set as the breed, and nothing is defined as strata for the object. To switch that to say the the species, we first need to take the @other info from microbov and make it into a data frame object, and then we can assign that data frame to the strata accession within the genind object
strata_df <- data.frame(other(microbov))
head(strata_df)
strata(microbov) <- strata_df
head(strata(microbov))
head(microbov@strata) # two different ways to access the @strata accessor
# Once we have the data frame stored in @strata, we can switch between strata by setting the @pop to different vectors within @strata
setPop(microbov) <- ~spe
microbov@pop
# And can always switch the @pop back to what it was before
setPop(microbov) <- ~breed
microbov@pop


### Optional depending on how you feel about doing these steps: For next week, use these strata to explore population-specific statistics for each of the three strata. 



####### Ok, now onto another big task in doing this kind of analysis. Dealing with data that isn't already loaded into R! Let's first walk through what typical microsatellite data would look like as a group. We'll focus on .str format, but realize there are a bunch of other formats for microsatellite genotypes (genepop, FSTAT, etc.). I'll provide you a few .str files to save somewhere on your computer. Then we'll read them into R and explore. These data are microsatellite genotypes from swallowtail butterflies in Alberta (AB) and British Columbia (BC) Canada, which I generated during my PhD!

# First, we'll need to show R where the files are that we want to read, and we do so by setting the "working directory" or "wd" with setwd(). For my computer that looks like:
setwd("/Users/julian/Library/CloudStorage/OneDrive-UniversityofKentucky/MBP_2023/UK/teaching/UFSCar_2023/R_popgen/")
# Next we'll use adegenet's read.structure() function to read these files into R.
??read.structure
# Lots of options to specify, but this function will actually walk you through some of these options if you don't supply any options. See if you can answer the questions that read.structure is asking you.
ABBC <- read.structure("ABBC_8Jan2015.str")
ABBC
# Alternatively, it's a little quicker to just specify that info in the function options.
ABBC2 <- read.structure("ABBC_8Jan2015.str", n.ind=781, n.loc=10, onerowperind=T,  row.marknames=NULL, col.lab=1, col.pop=2, NA.char="-9", ask=F)

# Let's explore this dataset a bit:
ABBC2
ABBC2@pop
indNames(ABBC2)

# The other file I gave you contains some additional information about these individuals and populations. Let's read it into R using read.table() and take a look.
ABBC_other <- read.table("ABBC_8Jan2015_other.txt", header=T)
head(ABBC_other)
# You can see it has the same individual names, their population codes, a region code, and lat and long. It'd be great to be able to use some of this info like we did for nancycats and microbov, so let's take this data and add it to the @other accessor for the ABBC2 object. First we can double check the @other tab is empty, then we can write this data to @other and double check that it worked. 
ABBC2
ABBC2@other
other(ABBC2) <- ABBC_other    # And then repeat the 2 previous lines

# We can set the strata by referring to the 2 columns in @other with pop and region. Remember the square brackets allow us to refer to rows and columns (in that order) within a dataframe, so pop and region are columns 2 and 3, and we can apply both by using a colon (i.e. take columns 2 to 3 from ABBC_other and set them as strata).
strata(ABBC2) <- ABBC_other[,2:3]
# And double check that it worked!
ABBC2
strata(ABBC2) 
# So now we can change the strata like we did for microbov.
setPop(ABBC2) <- ~region
ABBC2@pop


#### For next week, explore this swallowtail butterfly dataset with the stats that we have learned so far. How does genetic diversity differ between populations and regions? Are the populations differeniated? etc. etc. If the R bug has bitten you, you could even try to make an IBD plot for this dataset using the lat long information! Save any plots you make and we can go over them next week. 


########
########
######## Population genetic data analysis in R
######## Session 3: popgen on microsatellite data
########
########

# You know the drill
library(adegenet)
library(poppr)
library(pegas)
library(hierfstat)

# First, we'll be assessing some population structure with the swallowtail dataset, so let's load that up and get @other set up with strata. And let's use the updated version of the other file
setwd("/Users/julian/Library/CloudStorage/OneDrive-UniversityofKentucky/MBP_2023/UK/teaching/UFSCar_2023/R_popgen/")
ABBC2 <- read.structure("ABBC_8Jan2015.str", n.ind=781, n.loc=10, onerowperind=T,  row.marknames=NULL, col.lab=1, col.pop=2, NA.char="-9", ask=F)
ABBC_other <- read.table("ABBC_8Jan2015_other2.txt", header=T)
other(ABBC2) <- ABBC_other
strata(ABBC2) <- ABBC_other[,2:5]
	# useful just to check what's in @pop
	ABBC2@pop
	# and remember we can change that with
	setPop(ABBC2) <- ~pop_name

# Let's start with a Principal Components Analysis (PCA). Remember this is an unsupervised method for assessing population structure. First, we have to deal with the missing data in this dataset (missing alleles are noted as "-9"), since PCAs will not work with missing data. A standard way to deal with missing data in ordination methods is to replace any missing data values with the mean allele frequency for that locus. First we can count how many missing values are in this dataset
sum(is.na(ABBC2$tab))
# replace the missing values with mean frequencies
X <- tab(ABBC2, freq = TRUE, NA.method = "mean")
# So what is X?
class(X)
dim(X)
sum(is.na(X))

# Next, let's do the actual PCA on this object we just created. Some of the options in this statistical analyses are above our paygrade in this class, so I won't be going through all the options for each of these.
??dudi.pca
pca1 <- dudi.pca(X, scale = F, scannf = F, nf = 3)
pca1
# Lots of info in there. Of particular interest are:
	# $eig		which are the eigenvalues of the analysis (the amount of variance explained by each PC)
	# $li		the PCs themselves, which are synthetic variables summarizing the genetic diversity
	# $c1		Loadings of individual alleles, representing contributions to the PCs

# Let's think about these PCs a bit. So remember that the PCs are basically imaginary axes of variation in multivariate data. So these axes represent imaginary ranges of variation when you consider all the loci at once. PC1 explains the most variation, PC2 the second most, etc. etc.
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

# next up, let's look at the PCs, both in raw form, and then plotted
head(pca1$li)
s.label(pca1$li)
# and let's add the eigenvalue info to the plot
add.scatter.eig(pca1$eig[1:20], 3,1,2)
# So the PC values themselves are basically just where individuals are falling along the axes of variation for each PC. So you can basically treat pairs of those PC values as X,Y coordinates (that's what plotted in the plot). You can even look at different pairs than just 1,2
s.label(pca1$li,1,3)
# By default the first plot we made plotted PC1 vs. PC2 (pca1$li,1,2), and now we're plotting PC1 vs. PC3, so the placement of individuals on the Y axis has changed, but their placement on the X axis hasn't (although that's hard to tell with those big name boxes). Let's make that plot a bit more informative by putting population names on there with ellipses and lines to connect each population's individuals
s.class(pca1$li, pop(ABBC2))
# Now it's a bit easier to see the differences between PCs
s.class(pca1$li, pop(ABBC2),xax=1,yax=3)

# Here's where you can lose your mind in the plotting options in R. A simple example:
col <- funky(49)
s.class(pca1$li, pop(ABBC2),xax=1,yax=2, col=transp(col,.6), axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
title("PCA of ABBC swallowtails\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

# You can see that many of these populations are grouping together, so let's explore that regional information in the other info file. Let's set the strata to region_name and replot the PCA with the updated population strata
setPop(ABBC2) <- ~region_name
pop(ABBC2) # just to double check
s.class(pca1$li, pop(ABBC2))
# This dataset actually contains 2 species (zel and mach) from multiple regions, and in some areas these 2 species hybridize (particularly in parts of the Foothills region)


		# What about IBD in this dataset? The use of actual lat long values in R is a bit complicated, so this is just for funsies.
		library(hierfstat)
		setPop(ABBC2) <- ~pop
		genD <- genet.dist(ABBC2, method = "Fst") # calculate genetic distance
		latlong <- ABBC_other[c('lat','long')] # access lat long coordinates from the other file
		latlong_unique <- unique(latlong) # remove duplicates (since we're calculating distances for pairs of populations, not individuals)
		library(geosphere)
		library(sna)
		geoD_temp1 <- distm(latlong_unique[c('lat','long')],fun=distVincentyEllipsoid) # calculate a geographic distance from those lat longs
		geoD_temp2 <- geoD_temp1/1000 # turn it from meters to kilometers
		geoD <- as.dist(upper.tri.remove(geoD_temp2, remove.val=NA)) # remove upper triangle of matrix to get it into .dist class
		ibd <- mantel.randtest(genD,geoD, nrepet=10000) # calculate the ibd
		ibd
		plot(geoD,genD)
		dist_lm <- lm(as.vector(genD) ~ as.vector(geoD))
		abline(dist_lm, col="red", lty=2)

#####
# Ok, next onto a different type of data. Let's look at a vcf file containing SNPs! This is a ddRADseq dataset from a tropical tephritid fruit fly called the melon fly, Zeugodacus cucurbitae (called "Bactrocera cucurbitae" when this dataset was generated, hence the "bcur" file name).
library(vcfR)
vcf <- read.vcfR("bcur_world_only_0.2miss_0.01maf.recode.vcf")
snps <- vcfR2genind(vcf)
snps

# Again, we've got additional data for this species
snps_other <- read.table("bcur_other.txt", header=T)
other(snps) <- snps_other
strata(snps) <- snps_other
setPop(snps) <- ~pop    # vs. ~global
popNames(snps)

## Your mission, should you choose to accept it, is to explore this SNP dataset as we did for the nancycats, microbov, and ABBC swallowtail butterfly microsatellite datasets. ALTERNATIVELY, you could load in another SNP dataset (like one you've generated from WGS or RAD/GBS data) and analyze it yourself!
