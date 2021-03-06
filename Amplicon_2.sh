#!/bin/bash

#--------------------------SBATCH settings------

#SBATCH --job-name=Amp_2      ## job name
##SBATCH -A YOUR_LAB_ACCOUNT     ## account to charge
#SBATCH -p free          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=16    ## number of cores the job needs
##SBATCH --mail-user=MYEMAIL@uci.edu ## your email address
##SBATCH --mail-type=begin,end,fail ##type of emails to receive
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=slurm-%J.out ##output info file

#--- If there is a doube hash (##) before SBATCH then it is deactivated---

##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

#========Begin commands for job======

cd demultiplexed_seqs/*/data
module load R/3.6.2

echo "
library(dada2);
library(ggplot2);
library(utils);
path <- getwd();

# showing the program where to look for the files
fnFs <- sort(list.files(path, pattern='L001_R1_001.fastq', full.names = TRUE));
fnRs <- sort(list.files(path, pattern='L001_R2_001.fastq', full.names = TRUE));
sample.names <- sapply(strsplit(basename(fnFs), '_'), getElement, 1);

# making the filtered files to write to
filtFs <- file.path(path, 'filtered', paste0(sample.names, '_F_filt.fastq.gz'));
filtRs <- file.path(path, 'filtered', paste0(sample.names, '_R_filt.fastq.gz'));

# filtering, good idea to choose your values from Figaro!
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(296,167), maxEE=c(2,5), trimLeft=5, trimRight=5, rm.phix=TRUE, compress=TRUE, multithread=16);
saveRDS(out, file = 'out.rds')
write.csv(out, file = 'Filter-Trim.stats.csv')

# predicting error rate for the reads
errF <- dada2::learnErrors(filtFs, multithread=16, randomize=TRUE);
saveRDS(errF, file = 'errF.rds')
errR <- dada2::learnErrors(filtRs, multithread=16, randomize=TRUE);
saveRDS(errR, file = 'errR.rds')
Error_plot_F <- dada2::plotErrors(errF, nominalQ=TRUE);
Error_plot_R <- dada2::plotErrors(errR, nominalQ=TRUE);
ggsave('Error_plot_F.pdf', plot = Error_plot_F);
ggsave('Error_plot_R.pdf', plot = Error_plot_R);

# dereplicating reads to reduce complexity
derepFs <- dada2::derepFastq(filtFs, verbose=TRUE);
saveRDS(derepFs, file = 'derepFs.rds')
derepRs <- dada2::derepFastq(filtRs, verbose=TRUE);
saveRDS(derepRs, file = 'derepRs.rds')

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada2::dada(derepFs, err=errF, multithread=16, pool='pseudo');
saveRDS(dadaFs, file = 'dadaFs.rds')
dadaRs <- dada2::dada(derepRs, err=errR, multithread=16, pool='pseudo');
saveRDS(dadaRs, file = 'dadaRs.rds')

mergers <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE);
saveRDS(mergers, file = 'mergers.rds')
seqtab <- dada2::makeSequenceTable(mergers);

seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method='consensus', multithread=16, verbose=TRUE);
saveRDS(seqtab.nochim, file = 'seqtab_nochim.rds')

getN <- function(x) sum(getUniques(x));
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim));
colnames(track) <- c('input', 'filtered', 'denoisedF', 'denoisedR', 'merged', 'nonchim');
rownames(track) <- sample.names
write.csv(track, file = 'Dada2_stats_full.csv');

# Assign taxonomy
taxa <- dada2::assignTaxonomy(seqtab.nochim, '/dfs3b/whitesonlab/rdp_database/rdp_train_set_16.fa.gz', multithread=TRUE, minBoot = 70);
taxa <- addSpecies(taxa, '/dfs3b/whitesonlab/rdp_database/rdp_species_assignment_16.fa.gz')
saveRDS(taxa, file = 'taxa.rds')
write.csv(seqtab.nochim, 'OTU_table.csv');
write.csv(taxa, 'Species_taxa.csv');

dev.off();
 " | R --vanilla --no-save

gzip *.fastq

cp OTU_table.csv ../../../OTU_table.csv
cp Species_taxa.csv ../../../Species_taxa.csv
cp Dada2_stats_full.csv ../../../Dada2_stats_full.csv

# If you want to use the silva classifier:
# taxa <- dada2::assignTaxonomy(seqtab.nochim, '~/tax/silva_nr_v128_train_set.fa.gz', multithread=TRUE);
# taxa <- dada2::addSpecies(taxa, '~/tax/silva_species_assignment_v128.fa.gz');