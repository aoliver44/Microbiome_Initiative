#!/bin/bash

#--------------------------SBATCH settings------

#SBATCH --job-name=Amp_1      ## job name
##SBATCH -A YOUR_LAB_ACCOUNT     ## account to charge
#SBATCH -p free          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=6    ## number of cores the job needs
##SBATCH --mail-user=MYEMAIL@uci.edu ## your email address
##SBATCH --mail-type=begin,end,fail ##type of emails to receive
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=slurm-%J.out ##output info file

#--- If there is a doube hash (##) before SBATCH then it is deactivated---

##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

#========Begin commands for job======

module load qiime2/2020.8
source activate qiime2-2020.8


qiime tools import \
   --type EMPPairedEndSequences \
   --input-path raw_data \
   --output-path imported_data.qza

qiime demux emp-paired \
  --m-barcodes-file metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --i-seqs imported_data.qza \
  --o-per-sample-sequences demux.qza \
  --o-error-correction-details correction_details.qza \
  --p-rev-comp-barcodes \
  --p-rev-comp-mapping-barcodes

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# unzip the demux.qza and cd into it and then the data folder.
# unzip all the fastq.gz by gunzip *.fastq.gz
# place this script into the folder with all the unziped fastqs.

unzip demux.qza -d demultiplexed_seqs
cd demultiplexed_seqs/*/data
gunzip *.fastq.gz

#Clean as you go:
#rm ../../imported_data.qza
#rm ../../demux.qza
#rm ../../demux.qzv

conda deactivate
module purge
module load R/3.6.2

echo "
library(dada2);
library(ggplot2);
library(utils);
path <- getwd();
list.files(path);

fnFs <- sort(list.files(path, pattern='L001_R1_001.fastq', full.names = TRUE));
fnRs <- sort(list.files(path, pattern='L001_R2_001.fastq', full.names = TRUE));

sample.names <- sapply(strsplit(basename(fnFs), '_'), getElement, 1);

QC_forward <- dada2::plotQualityProfile(fnFs[1:16]);
QC_reverse <- dada2::plotQualityProfile(fnRs[1:16]);

ggsave('QC-Forward.pdf', plot = QC_forward);
ggsave('QC-Reverse.pdf', plot = QC_reverse);
dev.off();
" | R --no-save --vanilla

cp QC-Forward.pdf ../../../
cp QC-Reverse.pdf ../../../

# Setting up the HPC to run Dada2 through R is a little tricky. If i never used R 3.6.2 before,
# i had to fake install something like ggplot. i.e. install.packages('ggplot2'). This
# then asked for a CRAN mirror (i just picked the US (138 i think?)) and then it asked
# me if i wanted to install a personal 3.6.2 library. THIS IS VITAL. You hit yes. And i think
# yes again. Then you can install Dada2. Devtools is already installed.

# Install Dada2

#library(devtools)
#install.packages(ggplot2)