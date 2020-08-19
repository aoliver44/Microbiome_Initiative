#!/bin/bash
#$ -N Denoise
#$ -q bio,abio*,pub*,free*
#$ -pe openmp 8
#$ -R y
#$ -m beas
#$ -cwd
#$ -ckpt restart


module load R/3.5.1

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

# filtering
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,260), trimLeft=5, trimRight=5, rm.phix=TRUE, compress=TRUE, multithread=8);
saveRDS(out, file = 'out.rds')
write.csv(out, file = 'Filter-Trim.stats.csv')

# predicting error rate for the reads
errF <- dada2::learnErrors(filtFs, multithread=8, randomize=TRUE);
saveRDS(errF, file = 'errF.rds')
errR <- dada2::learnErrors(filtRs, multithread=8, randomize=TRUE);
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

dadaFs <- dada2::dada(derepFs, err=errF, multithread=8, pool='pseudo');
saveRDS(dadaFs, file = 'dadaFs.rds')
dadaRs <- dada2::dada(derepRs, err=errR, multithread=8, pool='pseudo');
saveRDS(dadaRs, file = 'dadaRs.rds')

mergers <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE);
saveRDS(mergers, file = 'mergers.rds')
seqtab <- dada2::makeSequenceTable(mergers);

seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method='consensus', multithread=8, verbose=TRUE);
saveRDS(seqtab.nochim, file = 'seqtab_nochim.rds')

getN <- function(x) sum(getUniques(x));
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim));
colnames(track) <- c('input', 'filtered', 'denoisedF', 'denoisedR', 'merged', 'nonchim');
rownames(track) <- sample.names
write.csv(track, file = 'Dada2_stats_full.csv');

# Assign taxonomy
taxa <- dada2::assignTaxonomy(seqtab.nochim, '/dfs5/bio/whitesonlab/rdp_database/rdp_train_set_16.fa.gz', multithread=TRUE, minBoot = 60);
taxa <- addSpecies(taxa, '/dfs5/bio/whitesonlab/rdp_database/rdp_species_assignment_16.fa.gz')
saveRDS(taxa, file = '../../../taxa.rds')
write.csv(seqtab.nochim, '../../../OTU_table.csv');
write.csv(taxa, 'Species_taxa.csv');

dev.off();
 " | R --vanilla --no-save
 
# If you want to use the silva classifier:
# taxa <- dada2::assignTaxonomy(seqtab.nochim, '~/tax/silva_nr_v128_train_set.fa.gz', multithread=TRUE);
# taxa <- dada2::addSpecies(taxa, '~/tax/silva_species_assignment_v128.fa.gz');