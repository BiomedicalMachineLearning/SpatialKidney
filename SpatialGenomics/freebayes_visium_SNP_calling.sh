#PBS -N 'frBayes_Visium_kidney5'
#PBS -S /bin/bash
#PBS -r n
#PBS -j oe
#PBS -l walltime=72:00:00
#PBS -A UQ-IMB
#PBS -l select=1:ncpus=8:mem=120GB


cd /90days/uqhnguy8/Spatial/QuanLabData/Visium/Visium5/Visium5_Kidney

source /sw/RCC/Anaconda/Anaconda3/5.2.0/etc/profile.d/conda.sh /etc/profile.d/conda.sh
conda activate /90days/uqhnguy8/Spatial/stLearn
freebayes -f genome.fa possorted_genome_bam.bam >FreeBayes_Visium5.vcf
