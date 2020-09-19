cd /30days/uqhnguy8/Visium/Visium5

#Format kidney disease files 

awk '{FS="\t"}{OFS="\t"}{print "chr"$1, $2, $2+1, $3, $4"."}' Kidney_disease_ralated_SNPs_bed4.bed >Kidney_disease_ralated_SNPs_bed4_V2.bed

grep -v chrNA Kidney_disease_ralated_SNPs_bed4_V2.bed >Kidney_disease_ralated_SNPs_bed4_V2_noNA.bed

awk '{print $4}' Kidney_disease_ralated_SNPs_chr_withNA.bed| sed 's/:/ /' | awk '{OFS="\t"}{print $1, $2, $2+1, ".", "."}' >Kidney_disease_ralated_SNPs_chr_withNA_v2.bed 


cat Kidney_disease_ralated_SNPs_chr_withNA_v2.bed  >>Kidney_disease_ralated_SNPs_bed4_V2_noNA.bed


sed 's/,\ /\_/g' Kidney_disease_ralated_SNPs_bed4_V2_noNA.bed >Kidney_disease_ralated_SNPs_bed4_V2_noNA_v2.bed


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

awk '{OFS="\t"} {print "chr"$1, $2, $3, $4, $5, $6}' FreeBayes_Vis5_D_bed5.vcf >FreeBayes_Vis5_D_bed5_chr.vcf

bedtools window -a Kidney_disease_ralated_SNPs_bed4_V2_noNA_v2_sorted.bed -b FreeBayes_Vis5_D_bed5_chr.vcf -w 10 >Vis5_D_SNP_overlaping_kidneyDiseaseSNPs.bed   
