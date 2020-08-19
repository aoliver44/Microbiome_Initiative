## Paired-End DADA2 Amplicon Analysis Pipeline
This pipeline is just one of a billion that is pretty good at analyzing 16SrRNA amplicons. I think its strength is that it is pretty clear what is going on under the hood (except demultiplexing...which i still use Qiime2 for). **It will run on the HPC2 cluster, and with some modifications should run on HPC3 just fine**

### Dependencies
* anaconda/2.7-4.3.1
* qiime2-2018.4
* R/3.5.1
    * library(dada2) (i run dada2_1.11.1)
        * You may have to install this, just like in RStudio
    * library(utils)
* pigz/2.3.1


### Step 1: Import and demultiplex
In your working directory, put a file called metadata.tsv, which should be formated exactly how Qiime2 requires a metadata file. You can find more information on their website. Next, make a sub-directory called raw_data. Place **ONLY** your sequence data in this folder called raw_data. The files should be named:
```
forward.fastq.gz
reverse.fastq.gz
barcodes.fastq.gz
```
Run the script titled Amplicon_1. This script will demultiplex your files and assess quality of the reads. It will generate a pdf of the first 16 samples quality graphs.

### Step 2: Everything else
After looking at the quality graphs, determine the quality values you wish to filter and trim your sequences. Then modify this line in Amplicon_2.sh accordingly:
```
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,260), trimLeft=5, trimRight=5, rm.phix=TRUE, compress=TRUE, multithread=8);
```
*more information can be found on the dada2 website*

Now run Amplicon_2.sh. This should take some time, depending on the amount of reads. The outputs will be several ".rds" files, which are great "save points". If something goes wrong in the pipeline, you dont lose your progress, you can just start at the last successful ".rds" file. Also output is a taxa file and an OTU table. 

### Step 3: Rarefy and make THE OTU table (in RSTUDIO probably)
R Dependencies (other versions may work):
* R version 3.4.4
* zoo v1.8-5
* tidyverse_1.3.0
* EcolUtils v0.1

Finally its time to make the final version of the OTU table. In this step, you will determine what cutoff should be made (# of sequences to rarefy to), and also we will attempt to make downstream analyses easier with a taxonomic "cheat".

``` 
##########################
### GENERATE BASIC ENV ###
##########################

library(tidyverse)
library(EcolUtils)
library(zoo)
setwd("~/path/to/taxa/and/otu")

# Import OTU file and rarefy
qiime_raw <- read.csv("OTU_table.csv", row.names = 1, check.names = FALSE, sep = ",")

rownames(qiime_raw) <- sapply(strsplit(basename(rownames(qiime_raw)), '_'), getElement, 1);

qiime_clean <- qiime_raw
# If your OTU file is HUGE, saving it and loading it as an
# rds file can save A LOT of time.
saveRDS(object = qiime_clean, file = "qiime_clean.rds")
# qiime_clean <- readRDS(file = "qiime_clean.rds")

# Import and clean up the taxonomy names
taxa <- read.csv("Species_taxa.csv", check.names = FALSE, sep = ",")
taxa <- data.frame(lapply(taxa, function(x) {gsub("\\[.*?\\]", "", x)}))
```

So that above is importing your OTU and Taxa table into R. Here is also that "cheat". What this block of code attempts to do is find taxonomic assignments that have nothing in them (i.e. f_Enterobacteriacea;g_;s_...notice how genus and species are blank). It will then look back to the last "known" taxonomic assignment (in this example, Enterobacteriaceae), and assign it through to species. It will also make sure each OTU is unique. So again in this example, the final phylogeny will be f_Enterobacteriacea;g_Enterobacteriacea;s_Enterobacteriacea_2 (the two being some number to make it unique). 

I think this is a GOLDEN way to do things, because when you go to do downstream analysis at higher taxonomic levels, often time your answer ends up being g_;s_. Cool cool cool. Not helpful. But with my trick, you at least have an indication of the highest taxonomy that was assigned. It also helps flesh out taxabarplots more. Many times Clostridiales dont get assigned pass order, and would then get lumped into the "other" catagory in a Genus-level barplot. Here they dont! You just have to modify the figure legend to make sure it is obvious that there is an Order in the genus-level barplot (i.e. Clostridiales(o))
 Viola:

```
# what this next part does is attempt to bring the last known phylgenetic group to unknown sections
#taxa <- taxa %>% separate(., col = Taxon, into = c("L1","L2","L3","L4","L5","L6","L7"), sep = "; ", remove = T, extra = "drop")
taxa <- data.frame(lapply(taxa, function(x) {gsub("[a-z]__", "", x)}))
taxa <- taxa %>% mutate_all(na_if,"")
taxa <- data.frame(t(apply(taxa,1,function(x) na.locf(x))))

# next up, you will replace the sequence string in OTU file
# with the taxa string in the taxa file
taxa$taxonomy <- paste0("k__", taxa$Kingdom,";","p__", taxa$Phylum,";","c__", taxa$Class,";","o__", taxa$Order,";","f__", taxa$Family,";","g__", taxa$Genus,";","s__", taxa$Species)
colnames(taxa)[1] <- "Sequence"
names(qiime_clean) <- taxa$taxonomy[match(names(qiime_clean),taxa$Sequence)]
names(qiime_clean) <- make.unique(colnames(qiime_clean))

metadata <- read.csv("metadata.tsv", check.names = FALSE, sep = "\t")
colnames(metadata)[1] <- "sample_id"
metadata <- metadata %>% filter(., sample_id !="#q2:types")
```

Finally it is time to rarefy:
```
# Before you rarefy, filter out chloroplasts and mitochondira
qiime_clean <- qiime_clean %>% select(-contains(c("Chloroplast", "Mitochondria"))

# look to see the level you want to rarefy at
barplot(sort(rowSums(qiime_clean)), ylim = c(0, max(rowSums(qiime_clean))), 
        xlim = c(0, NROW(qiime_clean)), col = "Blue", main = "Reads per Sample") 
sort(rowSums(qiime_clean))

# rarefy to 2000 sequences
set.seed(seed = 999)
rare_perm_otu <- rrarefy.perm(qiime_clean, sample = 2000, n = 15, round.out = T)
alpha_rare <- rare_perm_otu[rowSums(rare_perm_otu) >= 2000-(2000*.1), colSums(rare_perm_otu) >= 1]
barplot(sort(rowSums(alpha_rare)), ylim = c(0, max(rowSums(rare_perm_otu))), 
        xlim = c(0,NROW(rare_perm_otu)), col = "Blue", main = "Reads after rarefaction")
sort(rowSums(alpha_rare))


# merge OTU with taxa and metadata (the file i use most)
otu_merged <- merge(metadata, as.data.frame(alpha_rare), by.x = "sample_id", by.y = "row.names", all.x = F) 
```