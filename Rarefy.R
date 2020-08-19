##########################
### GENERATE BASIC ENV ###
##########################

library(tidyverse)
library(EcolUtils)
library(zoo)
setwd("~/path/to/files/")

# Import OTU file and rarefy...loading RDS object bc this part take like 10 min
# and im impatient

qiime_raw <- read.csv("OTU_table.csv", row.names = 1, check.names = FALSE, sep = ",")
rownames(qiime_raw) <- sapply(strsplit(basename(rownames(qiime_raw)), '_'), getElement, 1);
qiime_clean <- qiime_raw
saveRDS(object = qiime_clean, file = "qiime_clean.rds")

#qiime_clean <- readRDS(file = "qiime_clean.rds")

# Import and clean up the taxonomy
taxa <- read.csv("Species_taxa.csv", check.names = FALSE, sep = ",")
taxa <- data.frame(lapply(taxa, function(x) {gsub("\\[.*?\\]", "", x)}))

# what this next part does is attempt to bring the last known phylgenetic group to unknown sections
#taxa <- taxa %>% separate(., col = Taxon, into = c("L1","L2","L3","L4","L5","L6","L7"), sep = "; ", remove = T, extra = "drop")
taxa <- data.frame(lapply(taxa, function(x) {gsub("[a-z]__", "", x)}))
taxa <- taxa %>% mutate_all(na_if,"")
taxa <- data.frame(t(apply(taxa,1,function(x) na.locf(x))))

taxa$taxonomy <- paste0("k__", taxa$Kingdom,";","p__", taxa$Phylum,";","c__", taxa$Class,";","o__", taxa$Order,";","f__", taxa$Family,";","g__", taxa$Genus,";","s__", taxa$Species)
colnames(taxa)[1] <- "Sequence"
names(qiime_clean) <- taxa$taxonomy[match(names(qiime_clean),taxa$Sequence)]
names(qiime_clean) <- make.unique(colnames(qiime_clean))

metadata <- read.csv("metadata.tsv", check.names = FALSE, sep = "\t")
colnames(metadata)[1] <- "sample_id"
metadata <- metadata %>% filter(., sample_id != "#q2:types")

# filter out, pre rarefied table, of chloroplasts and mitochondira
qiime_clean <- qiime_clean %>% select(-contains(c("Chloroplast", "Mitochondria")))

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


# merge OTU with taxa and metadata
otu_merged <- merge(metadata, as.data.frame(alpha_rare), by.x = "sample_id", by.y = "row.names", all.x = F) 