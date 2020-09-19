#################
#This script has the workflow for formating and processing SNP from Visium data and those from GWAS Catalog database found associated with Kidney diseases 
#################

#download GWAS SNP catalog to get Kidney_disease_ralated_SNPs.txt, then process to get the bed format 
library(data.table)
dat <-fread("Kidney_disease_ralated_SNPs.txt")
colnames(dat)
dat_short <-dat[,c(12:13,22,14)]
head(dat_short)
write.table(dat_short, "Kidney_disease_ralated_SNPs_bed4.txt", sep="\t", col.names=F, row.names=F, quote=F)

#format kidney disease files to a bed file, add "chr"  

awk '{FS="\t"}{OFS="\t"}{print "chr"$1, $2, $2+1, $3, $4"."}' Kidney_disease_ralated_SNPs_bed4.bed >Kidney_disease_ralated_SNPs_bed4_V2.bed

grep -v chrNA Kidney_disease_ralated_SNPs_bed4_V2.bed >Kidney_disease_ralated_SNPs_bed4_V2_noNA.bed #remove chrNA 

sed 's/,\ /\_/g' Kidney_disease_ralated_SNPs_bed4_V2_noNA.bed >Kidney_disease_ralated_SNPs_bed4_V2_noNA_v2.bed #replace ","

#sort by chr  
bedtools sort -i Kidney_disease_ralated_SNPs_bed4_V2_noNA_v2.bed >Kidney_disease_ralated_SNPs_bed4_V2_noNA_v2_sorted.bed 

#format vcf file 
grep -v "#" FreeBayes_Vis5_A.vcf | awk '{OFS="\t"}{print $1, $2, $2+1, $4, $5, $6}' >FreeBayes_Vis5_A_bed5.vcf

grep -v "#" FreeBayes_Vis5_B.vcf | awk '{OFS="\t"}{print $1, $2, $2+1, $4, $5, $6}' >FreeBayes_Vis5_B_bed5.vcf

grep -v "#" FreeBayes_Vis5_C.vcf | awk '{OFS="\t"}{print $1, $2, $2+1, $4, $5, $6}' >FreeBayes_Vis5_C_bed5.vcf

grep -v "#" FreeBayes_Vis5_D.vcf | awk '{OFS="\t"}{print $1, $2, $2+1, $4, $5, $6}' >FreeBayes_Vis5_D_bed5.vcf

#find the overlapping 

module load  bedtools

bedtools window -a Kidney_disease_ralated_SNPs_bed4_V2_noNA_v2_sorted.bed -b FreeBayes_Vis5_A_bed5.vcf -w 10 >Vis5_A_SNP_overlaping_kidneyDiseaseSNPs.bed   

bedtools window -a Kidney_disease_ralated_SNPs_bed4_V2_noNA_v2_sorted.bed -b FreeBayes_Vis5_B_bed5.vcf -w 10 >Vis5_B_SNP_overlaping_kidneyDiseaseSNPs.bed   

bedtools window -a Kidney_disease_ralated_SNPs_bed4_V2_noNA_v2_sorted.bed -b FreeBayes_Vis5_C_bed5.vcf -w 10 >Vis5_C_SNP_overlaping_kidneyDiseaseSNPs.bed   

bedtools window -a Kidney_disease_ralated_SNPs_bed4_V2_noNA_v2_sorted.bed -b FreeBayes_Vis5_D_bed5_chr.vcf -w 10 >Vis5_D_SNP_overlaping_kidneyDiseaseSNPs.bed   
