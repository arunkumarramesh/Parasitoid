## full DGN analysis code, not meant to be run as a script

## first obtain the DGN data from SRA
prefetch --option-file SRR_Acc_List_Bergman.txt -O 2bergman/&
prefetch --option-file SRR_Acc_List_Clark.txt -O 2clark/&
prefetch --option-file SRR_Acc_List_Ages.txt -O 2ages/&
prefetch --option-file SRR_Acc_List_dpgp2.txt -O 2dpgp2&
prefetch --option-file SRR_Acc_List_dpgp3.txt -O 2dpgp3&
prefetch --option-file SRR_Acc_List_dgrp_sra2.txt -O 2dgrp&
prefetch --option-file SRR_Acc_List_nuzhdin1.txt -O 2nuzhdin&
prefetch --option-file SRR_Acc_List_pool.txt -O 2pool&

## Unpack sra files
for file in *.sra; do fastq-dump --gzip --split-files $file; done&
for file in *.fastq.gz; do mv $file ${file/.fastq.gz/.fq.gz}; done

## trim reads, pe and se separately
for file in *_1.fq.gz; do java -jar trimmomatic-0.36.jar PE -phred33 -threads 10 $file ${file/_1.fq.gz/_2.fq.gz} ${file/_1.fq.gz/_1.paired.fq.gz} ${file/_1.fq.gz/_1.unpaired.fq.gz} ${file/_1.fq.gz/_2.paired.fq.gz} ${file/_1.fq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:20:10:1:true TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 2>Trim3err; done&

for file in *_1.fq.gz; do java -jar trimmomatic-0.36.jar SE -phred33 -threads 10 $file  ${file/_1.fq.gz/_1.paired.fq.gz} ILLUMINACLIP:TruSeq3-SE.fa:2:20:10:1:true TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 2>Trim3err; done&

## map reads using bwa mem, pe and se separately
for file in *_1.paired.fq.gz; do bwa mem -t 10 dmel-all-chromosome-r5.13.fasta $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} 2>paired_map_err | samtools view -Sb - > ${file/_1.paired.fq.gz/_paired_wasp_mapped.bam} 2>paired_map_err; done&

for file in *_1.paired.fq.gz; do bwa mem -t 5 dmel-all-chromosome-r5.13.fasta $file 2>paired_map_err | samtools view -Sb - > ${file/_1.paired.fq.gz/_paired_wasp_mapped.bam} 2>paired_map_err; done&

## reorder reads
for file in *paired_wasp_mapped.bam; do picard-tools ReorderSam R=dmel-all-chromosome-r5.13.fasta I=$file O=${file/.bam/_reorder.bam} 2>reorder_err_1; done&

## add read groups
for file in *_reorder.bam ; do picard-tools AddOrReplaceReadGroups I=$file O=${file/.bam/_readgroup.bam} LB=species PL=illumina PU=1 SM=$file 2>readgroup_err; done&

## read filtering, only keep those with quality > 20 and for pe reads only keep properly paired reads
for file in *_readgroup.bam; do samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b $file >${file/.bam/_mq20.bam} 2> mq_john_err; done&
for file in *_readgroup.bam; do samtools view -q 20 -b $file >${file/.bam/_mq20.bam} 2> mq_john_err; done&

## sort reads
for file in *_mq20.bam; do samtools sort -@ 5 $file >${file/.bam/_sort.bam} 2>sort_err; done&

## mark duplicate reads
for file in *_sort.bam ; do picard-tools MarkDuplicates I=$file O=${file/.bam/_marked.bam} M=${file/.bam/_metrics.txt} 2>mark_err; done&

## build bam index
for file in *_marked.bam ; do picard-tools BuildBamIndex I=$file 2>index_err; done&

## realign reads near indels, some DGN require quality scores to be fixed
for file in *marked.bam; do java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R dmel-all-chromosome-r5.13.fasta --fix_misencoded_quality_scores -o ${file/.bam/.intervals} -I $file >realign_wasp_target_out_gen0 2> realign_wasp_target_err_gen0 ; done&

for file in *marked.bam; do java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R dmel-all-chromosome-r5.13.fasta -o ${file/.bam/.intervals} -I $file >realign_wasp_target_out_gen0 2> realign_wasp_target_err_gen0 ; done&

for file in *_marked.bam; do java -jar GenomeAnalysisTK.jar -T IndelRealigner --fix_misencoded_quality_scores --targetIntervals ${file/.bam/.intervals} -R dmel-all-chromosome-r5.13.fasta -I $file -o ${file/.bam/_realign.bam} >realign_wasp_target_out_gen0 2> realign_wasp_err; done&

for file in *_marked.bam; do java -jar GenomeAnalysisTK.jar -T IndelRealigner --targetIntervals ${file/.bam/.intervals} -R dmel-all-chromosome-r5.13.fasta -I $file -o ${file/.bam/_realign.bam} >realign_wasp_target_out_gen0 2> realign_wasp_err; done&

## build bam index
for file in *_sort.bam ; do picard-tools BuildBamIndex I=$file 2>index_err; done&

## FIRST DO LECTIN REGION ANALYSES
## extract lectin region files
for file in *_sort.bam ; do samtools view -b $file "2L:3716399-3719774" >${file/.bam/_lectin.bam}; done&

## call variants in lectin-24A using CRISP
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 2L:3716799-3718774 --VCF CRISP_dpgp_lectin24a.vcf

## annotate variants using snpEff and run dgn_lectin.R using output, conducts all the population genetics analyses
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_lectin24a.vcf > CRISP_dpgp_lectin24a_annot.vcf

## from dgn_lectin.R we get a lectin genotype and meta data files to create a custom made vcf file
paste lectin_gt lectin_gt | awk '{ n=NF/2; pad=""; for(i=1; i<=n; i++) { printf "%s%s/%s", pad, $i, $(i+n); pad=" "; } printf "\n"; }' > lectin_gt2
paste -d ' ' lectin_meta lectin_gt2 >lectin_meta_gt2
sed -i 's/ /\t/g' lectin_meta_gt2
head -n 1 lectin.vcf >lectin_header
## add in vcf header information
nano lectin_header
cat lectin_header lectin_meta_gt2 >lectin_masked.vcf
sed -i 's/NA\/NA/.\/./g'  lectin_masked.vcf
bgzip -c lectin_masked.vcf > lectin_masked.vcf.gz
tabix -p vcf lectin_masked.vcf.gz

## get per site fst for lectin region
vcftools --vcf lectin_masked.vcf --weir-fst-pop B_POP  --weir-fst-pop CO_POP --weir-fst-pop EA_POP --weir-fst-pop EB_POP --weir-fst-pop ED_POP --weir-fst-pop EF_POP --weir-fst-pop EG_POP --weir-fst-pop ER_POP --weir-fst-pop FR_POP --weir-fst-pop GA_POP --weir-fst-pop GU_POP --weir-fst-pop I_POP  --weir-fst-pop N_POP  --weir-fst-pop NG_POP --weir-fst-pop RAL_POP    --weir-fst-pop RG_POP --weir-fst-pop SB_POP --weir-fst-pop SD_POP --weir-fst-pop SF_POP --weir-fst-pop SP_POP --weir-fst-pop T_POP  --weir-fst-pop UG_POP --weir-fst-pop W_POP  --weir-fst-pop ZI_POP --weir-fst-pop ZS_POP --weir-fst-pop ZW_POP --out lectin_dpgp_stop_gain_masked_vcftools_POP

## remove indels from vcf file
sed '/indel/d' lectin_masked.vcf >lectin_masked_noindel.vcf
bgzip -c lectin_masked_noindel.vcf > lectin_masked_noindel.vcf.gz
tabix -p vcf lectin_masked_noindel.vcf.gz

# create a multi-sample fasta file for all DGN samples, Popover3 file containing sample analysed for lectin from dgn_lectin.R 
cat Popover3 | while read line; do samtools faidx dmel-all-chromosome-r5.13.fasta 2L:3716880-3717728 | bcftools/bcftools consensus -M N -s $line -p ${line/$/_} lectin_masked_noindel.vcf.gz >>lectin_masked.fa; done
sed -i 's/2L.*//' lectin_masked.fa

##reverse complement seq with https://www.bioinformatics.org/sms2/rev_comp.html, and get tree using MEGA

## pairwise lectin-24a coding region fst estimate using output from dgn_lectin.R 
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop asia.txt --weir-fst-pop ocenia.txt  --out res 2> asia_ocenia_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop asia.txt --weir-fst-pop america.txt  --out res 2> asia_america_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop asia.txt --weir-fst-pop europe_north_africa.txt  --out res 2> asia_europe_north_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop asia.txt --weir-fst-pop southern_africa.txt  --out res 2> asia_southern_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop asia.txt --weir-fst-pop central_africa.txt  --out res 2> asia_central_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop asia.txt --weir-fst-pop west_africa.txt  --out res 2> asia_west_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop asia.txt --weir-fst-pop east_africa.txt  --out res 2> asia_east_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop ocenia.txt --weir-fst-pop america.txt  --out res 2> ocenia_america_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop ocenia.txt --weir-fst-pop europe_north_africa.txt  --out res 2> ocenia_europe_north_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop ocenia.txt --weir-fst-pop southern_africa.txt  --out res 2> ocenia_southern_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop ocenia.txt --weir-fst-pop central_africa.txt  --out res 2> ocenia_central_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop ocenia.txt --weir-fst-pop west_africa.txt  --out res 2> ocenia_west_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop ocenia.txt --weir-fst-pop east_africa.txt  --out res 2> ocenia_east_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop america.txt --weir-fst-pop europe_north_africa.txt  --out res 2> america_europe_north_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop america.txt --weir-fst-pop southern_africa.txt  --out res 2> america_southern_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop america.txt --weir-fst-pop central_africa.txt  --out res 2> america_central_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop america.txt --weir-fst-pop west_africa.txt  --out res 2> america_west_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop america.txt --weir-fst-pop east_africa.txt  --out res 2> america_east_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop southern_africa.txt  --out res 2> europe_north_africa_southern_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop central_africa.txt  --out res 2> europe_north_africa_central_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop west_africa.txt  --out res 2> europe_north_africa_west_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop east_africa.txt  --out res 2> europe_north_africa_east_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop southern_africa.txt --weir-fst-pop central_africa.txt  --out res 2> southern_africa_central_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop southern_africa.txt --weir-fst-pop west_africa.txt  --out res 2> southern_africa_west_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop southern_africa.txt --weir-fst-pop east_africa.txt  --out res 2> southern_africa_east_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop central_africa.txt --weir-fst-pop west_africa.txt  --out res 2> central_africa_west_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop central_africa.txt --weir-fst-pop east_africa.txt  --out res 2> central_africa_east_africa_gfst
vcftools --chr 2L --from-bp 3716880 --to-bp 3717728 --vcf lectin_masked.vcf --weir-fst-pop west_africa.txt --weir-fst-pop east_africa.txt  --out res 2> west_africa_east_africa_gfst


## combine fst files into one
for file in *gfst; do sed -i -n '/^Weir/p' $file; done
for file in *gfst; do sed -n '/mean/p' $file >${file/_gfst/_gfst_mean} ; done
for file in *gfst; do sed -n '/weighted/p' $file >${file/_gfst/_gfst_weighted} ; done
for file in *_gfst_mean; do sed -i 's/Weir and Cockerham mean Fst estimate: //' $file ; done
for file in *_gfst_weighted; do sed -i 's/Weir and Cockerham weighted Fst estimate: //' $file ; done

## input for mapping in dgn_lectin.R 
grep "" *_gfst_mean >pairwise_mean_fst
grep "" *_gfst_weighted >pairwise_weighted_fst
sed -i -e 's/:/\t/' -e 's/_gfst_mean//' -e 's/_africa_/_africa\t/' -e 's/a_/a\t/' pairwise_mean_fst
sed -i -e 's/:/\t/' -e 's/_gfst_weighted//' -e 's/_africa_/_africa\t/' -e 's/a_/a\t/' pairwise_weighted_fst


## CONTINUE ANALYSES FOR FOR REMAINING GENOMIC READS (NOT LECTIN SPECIFIC BAMS)
## run CRISP to call variants in autosomes

CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 2L:1-6000000 --VCF CRISP_dpgp_2L_1-6000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 2L:6000001-12000000 --VCF CRISP_dpgp_2L_6000001-12000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 2L:12000001-18000000 --VCF CRISP_dpgp_2L_12000001-18000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 2L:18000001-23011544 --VCF CRISP_dpgp_2L_18000001-23011544.vcf&

CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3L:1-6000000 --VCF CRISP_dpgp_3L_1-6000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3L:6000001-12000000 --VCF CRISP_dpgp_3L_6000001-12000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3L:12000001-18000000 --VCF CRISP_dpgp_3L_12000001-18000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3L:18000001-24543557 --VCF CRISP_dpgp_3L_18000001-24543557.vcf&

CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 2R:1-5500000 --VCF CRISP_dpgp_2R_1-5500000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 2R:5500001-11000000 --VCF CRISP_dpgp_2R_5500001-11000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 2R:11000001-16000000 --VCF CRISP_dpgp_2R_11000001-16000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 2R:16000001-21146708 --VCF CRISP_dpgp_2R_16000001-21146708.vcf&

CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3R:1-5000000 --VCF CRISP_dpgp_3R_1-5000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3R:5000001-10000000 --VCF CRISP_dpgp_3R_5000001-10000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3R:10000001-15000000 --VCF CRISP_dpgp_3R_10000001-15000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3R:15000001-20000000 --VCF CRISP_dpgp_3R_15000001-20000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3R:20000001-24000000 --VCF CRISP_dpgp_3R_20000001-24000000.vcf&
CRISP --bams file.list --ref dmel-all-chromosome-r5.13.fasta -p 1 --regions 3R:24000001-27905053 --VCF CRISP_dpgp_3R_24000001-27905053.vcf&


## use snpEff to annotated variants
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_2L_12000001-18000000.vcf > CRISP_dpgp_2L_12000001-18000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_2L_1-6000000.vcf   > CRISP_dpgp_2L_1-6000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_2L_18000001-23011544.vcf > CRISP_dpgp_2L_18000001-23011544_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_2L_6000001-12000000.vcf  > CRISP_dpgp_2L_6000001-12000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_2R_11000001-16000000.vcf > CRISP_dpgp_2R_11000001-16000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_2R_1-5500000.vcf   > CRISP_dpgp_2R_1-5500000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_2R_16000001-21146708.vcf > CRISP_dpgp_2R_16000001-21146708_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_2R_5500001-11000000.vcf  > CRISP_dpgp_2R_5500001-11000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3L_12000001-18000000.vcf > CRISP_dpgp_3L_12000001-18000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3L_1-6000000.vcf   > CRISP_dpgp_3L_1-6000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3L_18000001-24543557.vcf > CRISP_dpgp_3L_18000001-24543557_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3L_6000001-12000000.vcf  > CRISP_dpgp_3L_6000001-12000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3R_10000001-15000000.vcf > CRISP_dpgp_3R_10000001-15000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3R_1-5000000.vcf   > CRISP_dpgp_3R_1-5000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3R_15000001-20000000.vcf > CRISP_dpgp_3R_15000001-20000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3R_20000001-24000000.vcf > CRISP_dpgp_3R_20000001-24000000_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3R_24000001-27905053.vcf > CRISP_dpgp_3R_24000001-27905053_ann.vcf&
java -Xmx4G -jar snpEff.jar dmel-all-chromosome-r5.13 CRISP_dpgp_3R_5000001-10000000.vcf  > CRISP_dpgp_3R_5000001-10000000_ann.vcf&

## create individual folders for further parallel analyses

mkdir CRISP_dpgp_2L_12000001-18000000_folder/
mkdir CRISP_dpgp_2L_1-6000000_folder/
mkdir CRISP_dpgp_2L_18000001-23011544_folder/
mkdir CRISP_dpgp_2L_6000001-12000000_folder/
mkdir CRISP_dpgp_2R_11000001-16000000_folder/
mkdir CRISP_dpgp_2R_1-5500000_folder/
mkdir CRISP_dpgp_2R_16000001-21146708_folder/
mkdir CRISP_dpgp_2R_5500001-11000000_folder/
mkdir CRISP_dpgp_3L_12000001-18000000_folder/
mkdir CRISP_dpgp_3L_1-6000000_folder/
mkdir CRISP_dpgp_3L_18000001-24543557_folder/
mkdir CRISP_dpgp_3L_6000001-12000000_folder/
mkdir CRISP_dpgp_3R_10000001-15000000_folder/
mkdir CRISP_dpgp_3R_1-5000000_folder/
mkdir CRISP_dpgp_3R_15000001-20000000_folder/
mkdir CRISP_dpgp_3R_20000001-24000000_folder/
mkdir CRISP_dpgp_3R_24000001-27905053_folder/
mkdir CRISP_dpgp_3R_5000001-10000000_folder/
mkdir CRISP_dpgp_X_11000001-16500000_folder/
mkdir CRISP_dpgp_X_1-5500000_folder/
mkdir CRISP_dpgp_X_16500001-22422827_folder/
mkdir CRISP_dpgp_X_5500001-11000000_folder/

## move vcf files and scripts and other metadata to process vcfs into subfolders

cp CRISP_dpgp_2L_12000001-18000000_ann.vcf CRISP_dpgp_2L_12000001-18000000_folder/ &
cp CRISP_dpgp_2L_1-6000000_ann.vcf CRISP_dpgp_2L_1-6000000_folder/ &
cp CRISP_dpgp_2L_18000001-23011544_ann.vcf CRISP_dpgp_2L_18000001-23011544_folder/ &
cp CRISP_dpgp_2L_6000001-12000000_ann.vcf CRISP_dpgp_2L_6000001-12000000_folder/&
cp CRISP_dpgp_2R_11000001-16000000_ann.vcf CRISP_dpgp_2R_11000001-16000000_folder/ &
cp CRISP_dpgp_2R_16000001-21146708_ann.vcf CRISP_dpgp_2R_16000001-21146708_folder/ &
cp CRISP_dpgp_3L_12000001-18000000_ann.vcf CRISP_dpgp_3L_12000001-18000000_folder/ &
cp CRISP_dpgp_3L_18000001-24543557_ann.vcf CRISP_dpgp_3L_18000001-24543557_folder/ &
cp CRISP_dpgp_3R_10000001-15000000_ann.vcf CRISP_dpgp_3R_10000001-15000000_folder/ &
cp CRISP_dpgp_3R_15000001-20000000_ann.vcf CRISP_dpgp_3R_15000001-20000000_folder/ &
cp CRISP_dpgp_3L_1-6000000_ann.vcf CRISP_dpgp_3L_1-6000000_folder/ &
cp CRISP_dpgp_3R_24000001-27905053_ann.vcf CRISP_dpgp_3R_24000001-27905053_folder/&
cp CRISP_dpgp_2R_5500001-11000000_ann.vcf CRISP_dpgp_2R_5500001-11000000_folder/&
cp CRISP_dpgp_3R_1-5000000_ann.vcf CRISP_dpgp_3R_1-5000000_folder/ &
cp CRISP_dpgp_3R_20000001-24000000_ann.vcf CRISP_dpgp_3R_20000001-24000000_folder/ &
cp CRISP_dpgp_3R_5000001-10000000_ann.vcf CRISP_dpgp_3R_5000001-10000000_folder/ &
cp CRISP_dpgp_3R_10000001-15000000_ann.vcf CRISP_dpgp_3R_10000001-15000000_folder/&
cp CRISP_dpgp_3L_6000001-12000000_ann.vcf CRISP_dpgp_3L_6000001-12000000_folder/ &

cp CRISP_dpgp_parse.R CRISP_dpgp_2L_12000001-18000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_2L_1-6000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_2L_18000001-23011544_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_2L_6000001-12000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_2R_11000001-16000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_2R_1-5500000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_2R_16000001-21146708_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_2R_5500001-11000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3L_12000001-18000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3L_1-6000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3L_18000001-24543557_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3L_6000001-12000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3R_10000001-15000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3R_1-5000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3R_15000001-20000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3R_20000001-24000000_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3R_24000001-27905053_folder/
cp CRISP_dpgp_parse.R CRISP_dpgp_3R_5000001-10000000_folder/

cp *_tracts.txt  CRISP_dpgp_2L_12000001-18000000_folder/
cp *_tracts.txt  CRISP_dpgp_2L_1-6000000_folder/
cp *_tracts.txt  CRISP_dpgp_2L_18000001-23011544_folder/
cp *_tracts.txt  CRISP_dpgp_2L_6000001-12000000_folder/
cp *_tracts.txt  CRISP_dpgp_2R_11000001-16000000_folder/
cp *_tracts.txt  CRISP_dpgp_2R_1-5500000_folder/
cp *_tracts.txt  CRISP_dpgp_2R_16000001-21146708_folder/
cp *_tracts.txt  CRISP_dpgp_2R_5500001-11000000_folder/
cp *_tracts.txt  CRISP_dpgp_3L_12000001-18000000_folder/
cp *_tracts.txt  CRISP_dpgp_3L_1-6000000_folder/
cp *_tracts.txt  CRISP_dpgp_3L_18000001-24543557_folder/
cp *_tracts.txt  CRISP_dpgp_3L_6000001-12000000_folder/
cp *_tracts.txt  CRISP_dpgp_3R_10000001-15000000_folder/
cp *_tracts.txt  CRISP_dpgp_3R_1-5000000_folder/
cp *_tracts.txt  CRISP_dpgp_3R_15000001-20000000_folder/
cp *_tracts.txt  CRISP_dpgp_3R_20000001-24000000_folder/
cp *_tracts.txt  CRISP_dpgp_3R_24000001-27905053_folder/
cp *_tracts.txt  CRISP_dpgp_3R_5000001-10000000_folder/
cp *_tracts.txt  CRISP_dpgp_X_11000001-16500000_folder/
cp *_tracts.txt  CRISP_dpgp_X_1-5500000_folder/
cp *_tracts.txt  CRISP_dpgp_X_16500001-22422827_folder/
cp *_tracts.txt  CRISP_dpgp_X_5500001-11000000_folder/

cp SraRunTable*  CRISP_dpgp_2L_12000001-18000000_folder/
cp SraRunTable*  CRISP_dpgp_2L_1-6000000_folder/
cp SraRunTable*  CRISP_dpgp_2L_18000001-23011544_folder/
cp SraRunTable*  CRISP_dpgp_2L_6000001-12000000_folder/
cp SraRunTable*  CRISP_dpgp_2R_11000001-16000000_folder/
cp SraRunTable*  CRISP_dpgp_2R_1-5500000_folder/
cp SraRunTable*  CRISP_dpgp_2R_16000001-21146708_folder/
cp SraRunTable*  CRISP_dpgp_2R_5500001-11000000_folder/
cp SraRunTable*  CRISP_dpgp_3L_12000001-18000000_folder/
cp SraRunTable*  CRISP_dpgp_3L_1-6000000_folder/
cp SraRunTable*  CRISP_dpgp_3L_18000001-24543557_folder/
cp SraRunTable*  CRISP_dpgp_3L_6000001-12000000_folder/
cp SraRunTable*  CRISP_dpgp_3R_10000001-15000000_folder/
cp SraRunTable*  CRISP_dpgp_3R_1-5000000_folder/
cp SraRunTable*  CRISP_dpgp_3R_15000001-20000000_folder/
cp SraRunTable*  CRISP_dpgp_3R_20000001-24000000_folder/
cp SraRunTable*  CRISP_dpgp_3R_24000001-27905053_folder/
cp SraRunTable*  CRISP_dpgp_3R_5000001-10000000_folder/
cp SraRunTable*  CRISP_dpgp_X_11000001-16500000_folder/
cp SraRunTable*  CRISP_dpgp_X_1-5500000_folder/
cp SraRunTable*  CRISP_dpgp_X_16500001-22422827_folder/
cp SraRunTable*  CRISP_dpgp_X_5500001-11000000_folder/

cp TableS*  CRISP_dpgp_2L_12000001-18000000_folder/
cp TableS*  CRISP_dpgp_2L_1-6000000_folder/
cp TableS*  CRISP_dpgp_2L_18000001-23011544_folder/
cp TableS*  CRISP_dpgp_2L_6000001-12000000_folder/
cp TableS*  CRISP_dpgp_2R_11000001-16000000_folder/
cp TableS*  CRISP_dpgp_2R_1-5500000_folder/
cp TableS*  CRISP_dpgp_2R_16000001-21146708_folder/
cp TableS*  CRISP_dpgp_2R_5500001-11000000_folder/
cp TableS*  CRISP_dpgp_3L_12000001-18000000_folder/
cp TableS*  CRISP_dpgp_3L_1-6000000_folder/
cp TableS*  CRISP_dpgp_3L_18000001-24543557_folder/
cp TableS*  CRISP_dpgp_3L_6000001-12000000_folder/
cp TableS*  CRISP_dpgp_3R_10000001-15000000_folder/
cp TableS*  CRISP_dpgp_3R_1-5000000_folder/
cp TableS*  CRISP_dpgp_3R_15000001-20000000_folder/
cp TableS*  CRISP_dpgp_3R_20000001-24000000_folder/
cp TableS*  CRISP_dpgp_3R_24000001-27905053_folder/
cp TableS*  CRISP_dpgp_3R_5000001-10000000_folder/
cp TableS*  CRISP_dpgp_X_11000001-16500000_folder/
cp TableS*  CRISP_dpgp_X_1-5500000_folder/
cp TableS*  CRISP_dpgp_X_16500001-22422827_folder/
cp TableS*  CRISP_dpgp_X_5500001-11000000_folder/

## modify scripts to accept correct input files
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_2L_18000001-23011544/' CRISP_dpgp_2L_18000001-23011544_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_2R_11000001-16000000/' CRISP_dpgp_2R_11000001-16000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_2R_1-5500000/' CRISP_dpgp_2R_1-5500000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_2R_16000001-21146708/' CRISP_dpgp_2R_16000001-21146708_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_2R_5500001-11000000/' CRISP_dpgp_2R_5500001-11000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3L_12000001-18000000/' CRISP_dpgp_3L_12000001-18000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3L_1-6000000/' CRISP_dpgp_3L_1-6000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_2L_12000001-18000000/' CRISP_dpgp_2L_12000001-18000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_2L_1-6000000/' CRISP_dpgp_2L_1-6000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_2L_6000001-12000000/' CRISP_dpgp_2L_6000001-12000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3L_18000001-24543557/' CRISP_dpgp_3L_18000001-24543557_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3L_6000001-12000000/' CRISP_dpgp_3L_6000001-12000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3R_10000001-15000000/' CRISP_dpgp_3R_10000001-15000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3R_1-5000000/' CRISP_dpgp_3R_1-5000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3R_15000001-20000000/' CRISP_dpgp_3R_15000001-20000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3R_20000001-24000000/' CRISP_dpgp_3R_20000001-24000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3R_5000001-10000000/' CRISP_dpgp_3R_5000001-10000000_folder/CRISP_dpgp_parse.R
sed -i 's/CRISP_dpgp_3R_24000001-27905053/CRISP_dpgp_3R_24000001-27905053/' CRISP_dpgp_3R_24000001-27905053_folder/CRISP_dpgp_parse.R

## run the vcf process script in parallel across multiple folders, an example of this script is in Github. Purpose of the CRISP_dpgp_parse.R is to filter snps (e.g. missing data, allele frequency, no indels, biallelic)

cd CRISP_dpgp_2L_18000001-23011544_folder/
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_2R_11000001-16000000_folder/
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_2R_1-5500000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_2R_16000001-21146708_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_2R_5500001-11000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3L_12000001-18000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3L_1-6000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_2L_12000001-18000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_2L_1-6000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_2L_6000001-12000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3L_18000001-24543557_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3L_6000001-12000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3R_10000001-15000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3R_1-5000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3R_15000001-20000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3R_20000001-24000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3R_5000001-10000000_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&
cd ../CRISP_dpgp_3R_24000001-27905053_folder
rm CRISP_dpgp_*_parse_ann_*
rm CRISP_dpgp_*_parse_geno_*
Rscript CRISP_dpgp_parse.R&

## process the output in each folder separately, combine files that were output in 1000 snp chunks from CRISP_dpgp_parse.R

cd ../CRISP_dpgp_2L_1-6000000_folder/
ls CRISP_dpgp_2L_1-6000000_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2L_1-6000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2L_1-6000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_2L_1-6000000_parse_ann_all

cd ../CRISP_dpgp_2L_6000001-12000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2L_6000001-12000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2L_6000001-12000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_2L_12000001-18000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2L_12000001-18000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2L_12000001-18000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_2L_18000001-23011544_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2L_18000001-23011544_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2L_18000001-23011544_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_2R_1-5500000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2R_1-5500000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2R_1-5500000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_2R_5500001-11000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2R_5500001-11000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2R_5500001-11000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_2R_11000001-16000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2R_11000001-16000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2R_11000001-16000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_2R_16000001-21146708_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2R_16000001-21146708_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_2R_16000001-21146708_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3L_1-6000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3L_1-6000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3L_1-6000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3L_12000001-18000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3L_6000001-12000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3L_6000001-12000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3L_12000001-18000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3L_12000001-18000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3L_12000001-18000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3L_18000001-24543557_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3L_18000001-24543557_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3L_18000001-24543557_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3R_1-5000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_1-5000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_1-5000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3R_5000001-10000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_5000001-10000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_5000001-10000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3R_510000001-15000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_510000001-15000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_10000001-15000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3R_15000001-20000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_15000001-20000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_15000001-20000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3R_20000001-24000000_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_20000001-24000000_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_20000001-24000000_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all

cd ../CRISP_dpgp_3R_24000001-27905053_folder/
ls CRISP_dpgp_*_parse_geno_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_24000001-27905053_parse_geno_all' | tr '\n' ' ' >cat_geno.sh
chmod 777 cat_geno.sh
./cat_geno.sh
ls CRISP_dpgp_*_parse_ann_* |sort -V  | sed -e '1i cat' -e '$a>CRISP_dpgp_3R_24000001-27905053_parse_ann_all' | tr '\n' ' ' >cat_ann.sh
chmod 777 cat_ann.sh
./cat_ann.sh
sed -i '/V1/d' CRISP_dpgp_*_parse_ann_all


## now in parent folder, obtain all genotype calls from subfolders and combine them
cd ../
find . -name "*_parse_geno_all" -type f -print0 | xargs -0 cp -t .
head -n 1 header_CRISP_dpgp_X_5500001-11000000_ >header
ls  CRISP_dpgp_*_parse_geno_all | sort -V | sed -e '1i cat header' -e '$a>CRISP_dpgp_Genotypes' | tr '\n' ' ' >cat_Geno.sh
chmod 777 cat_Geno.sh
./cat_Geno.sh
ls CRISP_dpgp_*_parse_geno_all | sort -V | sed -e '1i cat' -e '$a>tester' | tr '\n' ' ' >cat_Geno2.sh
chmod 777 cat_Geno2.sh
./cat_Geno2.sh

##  in parent folder, obtain all annotations calls from subfolders and combine them, need to add chromosome names for annotation files, missed this in script
find . -name "*_parse_ann_all" -type f -print0 | xargs -0 cp -t .
sed -i 's/^/2L /' 	CRISP_dpgp_2L_12000001-18000000_parse_ann_all
sed -i 's/^/2L /' 	CRISP_dpgp_2L_1-6000000_parse_ann_all
sed -i 's/^/2L /' 	CRISP_dpgp_2L_18000001-23011544_parse_ann_all
sed -i 's/^/2L /' 	CRISP_dpgp_2L_6000001-12000000_parse_ann_all
sed -i 's/^/2R /' 	CRISP_dpgp_2R_11000001-16000000_parse_ann_all
sed -i 's/^/2R /' 	CRISP_dpgp_2R_1-5500000_parse_ann_all
sed -i 's/^/2R /' 	CRISP_dpgp_2R_16000001-21146708_parse_ann_all
sed -i 's/^/2R /' 	CRISP_dpgp_2R_5500001-11000000_parse_ann_all
sed -i 's/^/3L /' 	CRISP_dpgp_3L_12000001-18000000_parse_ann_all
sed -i 's/^/3L /' 	CRISP_dpgp_3L_1-6000000_parse_ann_all
sed -i 's/^/3L /' 	CRISP_dpgp_3L_18000001-24543557_parse_ann_all
sed -i 's/^/3L /' 	CRISP_dpgp_3L_6000001-12000000_parse_ann_all
sed -i 's/^/3R /' 	CRISP_dpgp_3R_10000001-15000000_parse_ann_all
sed -i 's/^/3R /' 	CRISP_dpgp_3R_1-5000000_parse_ann_all
sed -i 's/^/3R /' 	CRISP_dpgp_3R_15000001-20000000_parse_ann_all
sed -i 's/^/3R /' 	CRISP_dpgp_3R_20000001-24000000_parse_ann_all
sed -i 's/^/3R /' 	CRISP_dpgp_3R_24000001-27905053_parse_ann_all
sed -i 's/^/3R /' 	CRISP_dpgp_3R_5000001-10000000_parse_ann_all
ls  CRISP_dpgp_*_parse_ann_all | sort -V | sed -e '1i cat' -e '$a>CRISP_dpgp_Annotations' | tr '\n' ' ' >cat_Anno.sh
chmod 777 cat_Anno.sh
./cat_Anno.sh
sed -i '/V1/d' CRISP_dpgp_Annotations

## process the collected genotype files into a vcf compatible coding
cut -d ' ' -f 2- tester >tester_dat
sed -i 's/NA/./g' tester_dat
paste tester_dat tester_dat | awk '{ n=NF/2; pad=""; for(i=1; i<=n; i++) { printf "%s%s/%s", pad, $i, $(i+n); pad=" "; } printf "\n"; }' > tester_dat2
sed -i 's/^/. T . . . GT /' tester_dat2
cut -d' ' -f 1 tester | paste -d' ' - tester_dat2 >tester_dat3
sed -i 's/_/ /' tester_dat3
sed -i 's/_/ /' tester_dat3
sed -i 's/ /\t/g'  tester_dat3
# paste in a vcf header into header_tester
nano header_tester
# sample names taken from a subfolder output
cp CRISP_dpgp_2L_1-6000000_folder/header_CRISP_dpgp_2L_1-6000000_ .
head -n 1 header_CRISP_dpgp_2L_1-6000000_ >header_tester
nano header_tester
sed -i 's/ /\t/g' header_tester
cat header_tester tester_dat3 >tester_wh.vcf

## from original (CRISP called) vcfs, get the snp annotations
for file in *_ann.vcf; do tail -n +28 $file | cut -f 1-6 >>vcf_meta_data; done&
sed -i 's/\t/_/' vcf_meta_data
sed -i 's/\t/_/' vcf_meta_data

## get a list of all snps that remain after filtering
cut -d '_' -f 1-2 CRISP_dpgp_Genotypes | tail -n +2 >snp_ids
sed -i 's/$/_/' snp_ids

## subset snp annotations for filtered snps only
grep -F -f snp_ids vcf_meta_data >vcf_meta_data_snps
sed -i 's/_/\t/' vcf_meta_data_snps
sed -i 's/_/\t/' vcf_meta_data_snps

## get meta data for filtered snps, w/o effect annotations
cut -f 1-9 tester_wh.vcf >tester_wh_meta

## now this script combine effect annotations with filtered custom made vcf meta data, annotations are in third column, output file in dpgp_masked
Rscript meta_add_refalt.R

## extracts just gt calls
tail -n +9 tester_wh.vcf | cut -f 10- >tester_wh_gt
## extracts header column
head -n 8 tester_wh.vcf >tester_wh_header

## combine dpgp_masked containing effect annotations and other meta data from  meta_add_refalt.R and all gt calls
paste dpgp_masked tester_wh_gt >dpgp_masked_gt
## identify and extarct unique calls
sort -V dpgp_masked_gt | uniq  >dpgp_masked_gt_uniq
cat tester_wh_header dpgp_masked_gt_uniq >dpgp_masked.vcf

## to read sample names into R
head -n 8 dpgp_masked.vcf | tail -n +8 | sed 's/#CHROM/CHROM/g' >dpgp_masked2_header
head -n 7 tester_wh_header >tester_wh_header2


## split filtered vcf file into parts for faster processing by vcftools
split -n l/15 dpgp_masked.vcf
cp xaa xaa_h

## add vcf header to part files
cat tester_wh_header xab >xab_h
cat tester_wh_header xac >xac_h
cat tester_wh_header xad >xad_h
cat tester_wh_header xae >xae_h
cat tester_wh_header xaf >xaf_h
cat tester_wh_header xag >xag_h
cat tester_wh_header xah >xah_h
cat tester_wh_header xai >xai_h
cat tester_wh_header xaj >xaj_h
cat tester_wh_header xak >xak_h
cat tester_wh_header xal >xal_h
cat tester_wh_header xam >xam_h
cat tester_wh_header xan >xan_h
cat tester_wh_header xao >xao_h

## dpgp_fst_vcftools_part1 is the script to calculate fst, modify file to accept part input files
sed 's/xaa/xab/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part2.R
sed 's/xaa/xac/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part3.R
sed 's/xaa/xad/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part4.R
sed 's/xaa/xae/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part5.R
sed 's/xaa/xaf/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part6.R
sed 's/xaa/xag/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part7.R
sed 's/xaa/xah/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part8.R
sed 's/xaa/xai/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part9.R
sed 's/xaa/xaj/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part10.R
sed 's/xaa/xak/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part11.R
sed 's/xaa/xal/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part12.R
sed 's/xaa/xam/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part13.R
sed 's/xaa/xan/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part14.R
sed 's/xaa/xao/g' dpgp_fst_vcftools_part1.R >dpgp_fst_vcftools_part15.R

# run fst script
Rscript dpgp_fst_vcftools_part1.R >>part1_out 2>part1_err&
Rscript dpgp_fst_vcftools_part2.R >>part2_out 2>part2_err&
Rscript dpgp_fst_vcftools_part3.R >>part3_out 2>part3_err&
Rscript dpgp_fst_vcftools_part4.R >>part4_out 2>part4_err&
Rscript dpgp_fst_vcftools_part5.R >>part5_out 2>part5_err&
Rscript dpgp_fst_vcftools_part6.R >>part6_out 2>part6_err&
Rscript dpgp_fst_vcftools_part7.R >>part7_out 2>part7_err&
Rscript dpgp_fst_vcftools_part8.R >>part8_out 2>part8_err&
Rscript dpgp_fst_vcftools_part9.R >>part9_out 2>part9_err&
Rscript dpgp_fst_vcftools_part10.R >>part10_out 2>part10_err&
Rscript dpgp_fst_vcftools_part11.R >>part11_out 2>part11_err&
Rscript dpgp_fst_vcftools_part12.R >>part12_out 2>part12_err&
Rscript dpgp_fst_vcftools_part13.R >>part13_out 2>part13_err&
Rscript dpgp_fst_vcftools_part14.R >>part14_out 2>part14_err&
Rscript dpgp_fst_vcftools_part15.R >>part15_out 2>part15_err&

## combine fst output
cat fst_new_header dpgp_vcftools_POP.weir.fst_xa* >dpgp_vcftools_POP.weir.fst

## get lectin region specific fst
awk -F "\t" '{ if ($1 == "2L" &&  $2 < 3718050 &&  $2 > 3716812) { print } }' dpgp_vcftools_POP.weir.fst >dpgp_vcftools_POP_lectinregion.weir.fst

## make a haploid vcf file for snps where fst was calculated
sed -e 's/\/.//g'  -e 's/\./NA/g' dpgp_masked.vcf > dpgp_masked_hap.vcf 
sed  -e '/nan/d' dpgp_vcftools_POP.weir.fst |  cut -f 1-2 | sed -e 's/\t/_/g' -e 's/$/_/g' -e '/POS/d' >fst_snps
sed -i -e 's/\t/_/' -e 's/\t/_/' dpgp_masked_hap.vcf
grep -F -f fst_snps dpgp_masked_hap.vcf >dpgp_masked_hap_fst.vcf
sed -i -e 's/_/\t/' -e 's/_/\t/' dpgp_masked_hap_fst.vcf
cat tester_wh_header dpgp_masked_hap_fst.vcf >dpgp_masked_hap_fst_header.vcf

## get stop codons for filtered vcf file
grep '\sstop_gained\s\|\sstop_lost\s' dpgp_masked.vcf >dpgp_stops_3
cat tester_wh_header dpgp_stops_3 >dpgp_stops3.vcf

## process vcf containing stop codons for R compatibility
sed -e 's/\/.//g'  -e 's/\./NA/g' dpgp_stops3.vcf >dpgp_stops2.vcf
sed -i 's/4NA1/4.1/' dpgp_stops2.vcf
sed -i 's/#CHROM/CHROM/' dpgp_stops2.vcf

## calculate pairwise fst for stop codons for samples from three geographic regions
vcftools --vcf dpgp_stops3.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop southern_africa.txt  --out ENA_SA_stop 
vcftools --vcf dpgp_stops3.vcf --weir-fst-pop america.txt  --weir-fst-pop southern_africa.txt  --out NA_SA_stop 
vcftools --vcf dpgp_stops3.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop america.txt   --out ENA_NA_stop 

## also repeat for the coding deletion vcf that was created separately
vcftools --vcf coding_del.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop southern_africa.txt  --out ENA_SA_codingdel
vcftools --vcf coding_del.vcf --weir-fst-pop america.txt  --weir-fst-pop southern_africa.txt  --out NA_SA_codingdel 
vcftools --vcf coding_del.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop america.txt   --out ENA_NA_codingdel

## calculate pairwise fst across all stop codons for samples from three geographic regions
vcftools --vcf dpgp_stops3.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop southern_africa.txt  --out ENA_SA_stop_filtered  --positions stop_codon_list
vcftools --vcf dpgp_stops3.vcf --weir-fst-pop america.txt  --weir-fst-pop southern_africa.txt  --out NA_SA_stop_filtered  --positions stop_codon_list
vcftools --vcf dpgp_stops3.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop america.txt   --out ENA_NA_stop_filtered  --positions stop_codon_list
grep 'Weir and Cockerham mean Fst estimate' *_stop_filtered.log | sed 's/_stop_filtered.log:Weir and Cockerham mean Fst estimate://' >all_stop_fst.txt

## get stop codon positions for recombination analyses
awk '{print $1 ":" ($2) ".." ($2)}' dpgp_stops_pos >dpgp_stops_r

### get recombination rate estimates for regions containing stop codons
perl RRC-open-v2.3.pl -M dpgp_stops_r
sed 's/:/\t/' dpgp_stops_r.rrc >dpgp_stops_r2.rrc


## get stop codon position and positions of 100bp upstream and downstream
sed -e '/#/d' -e '/CHROM/d' dpgp_stops3.vcf | cut -f 1-2 >dpgp_stops_pos
awk '{print $1 ":" ($2 - 100) "-" ($2 + 100)}' dpgp_stops_pos >dpgp_stops_regions

## split location files into 2 since too much for a single BLAST search
head -n 3000 dpgp_stops_regions >dpgp_stops_regions_part1
tail -n +3001 dpgp_stops_regions >dpgp_stops_regions_part2

## extract sequence from reference genome
samtools faidx dmel-all-chromosome-r5.13.fasta -r dpgp_stops_regions_part1 >dpgp_stops_dmel_part1.fa
samtools faidx dmel-all-chromosome-r5.13.fasta -r dpgp_stops_regions_part2 >dpgp_stops_dmel_part2.fa

## blast again D. simulans and D.sech on NCBI BLAST page, not flybase, get txt output, megablast
## combine output txt files for the two parts and extract query and subject alignments
cat JCNVN0KT013-Alignment.txt JCNV2M2M013-Alignment.txt >Dsech_alns
grep 'Query\|Sbjct' Dsech_alns | sed -e '/Scientific/d' -e 's/ID.*//'   -e 's/#//g' >Dsech_alns_sub

cat JCNN059Z016-Alignment.txt JCNSEMH5013-Alignment.txt >Dsim_alns
grep 'Query\|Sbjct' Dsim_alns | sed -e '/Scientific/d' -e 's/ID.*//'   -e 's/#//g' >Dsim_alns_sub

## this section gets the short intron variants for the PBS analyses
grep 'intron' dmel-all-r5.13.gff >dmel-all-r5.13_intron.gff  ## get short intron region and process and extract required region in dgn_lectin.R

## next part is to identify sites typed in at least 40% samplers per region and with MAF > 0.05% across all samples
cat america.txt europe_north_africa.txt southern_africa.txt >allthreeregions
vcftools --vcf coding_del.vcf --keep allthreeregions.txt --recode-INFO-all --out coding_del_3regions --recode
vcftools --vcf dpgp_masked.vcf --keep allthreeregions.txt --recode-INFO-all --out dpgp_masked_3regions --recode
bgzip coding_del_3regions.recode.vcf
tabix coding_del_3regions.recode.vcf.gz
bgzip dpgp_masked_3regions.recode.vcf
tabix dpgp_masked_3regions.recode.vcf.gz
bcftools concat  -o dpgp_masked_3regions_all.vcf -O v dpgp_masked_3regions.recode.vcf.gz coding_del_3regions.recode.vcf.gz
## only use SNPs with MAF of 0.05
vcftools --vcf dpgp_masked_3regions_all.vcf --recode --out dpgp_masked_3regions_0.05 --maf 0.05

## this part is a filter to removes sites with more than 40% missing data in each regions
vcftools --vcf dpgp_masked_3regions_0.05.recode.vcf --recode --out dpgp_masked_3regions_forindep_SA --max-missing 0.4  --keep southern_africa.txt
vcftools --vcf dpgp_masked_3regions_0.05.recode.vcf --recode --out dpgp_masked_3regions_forindep_NA --max-missing 0.4  --keep america.txt
vcftools --vcf dpgp_masked_3regions_0.05.recode.vcf --recode --out dpgp_masked_3regions_forindep_ENA --max-missing 0.4  --keep europe_north_africa.txt

## now get the filtered snps per regions
cut -f 1-2 dpgp_masked_3regions_forindep_ENA.recode.vcf >ENA_snp_list
cut -f 1-2 dpgp_masked_3regions_forindep_NA.recode.vcf >NA_snp_list
cut -f 1-2 dpgp_masked_3regions_forindep_SA.recode.vcf >SA_snp_list

## Rscript to only keep sites that intersect all 3 regions - generates file 'alllist'
Rscript combine_snps.R

## this part gets the regional frequency for short intron sites
vcftools --vcf dpgp_masked.vcf  --out dpgp_masked_introns --bed introns.txt --recode

## this part gets the regional frequency for short intron sites typed in at least 40% samplers per region and with MAF > 0.05% across all samples
vcftools --vcf dpgp_masked_introns.recode.vcf --freq2 --keep europe_north_africa.txt --out ENA_introns_freq  --positions alllist
vcftools --vcf dpgp_masked_introns.recode.vcf  --freq2 --keep america.txt --out NA_introns_freq --positions alllist
vcftools --vcf dpgp_masked_introns.recode.vcf  --freq2 --keep southern_africa.txt --out SA_introns_freq --positions alllist

## this part gets the pairwise fst for short intron sites typed in at least 40% samplers per region and with MAF > 0.05% across all samples
vcftools --vcf dpgp_masked_introns.recode.vcf  --positions alllist --weir-fst-pop europe_north_africa.txt --weir-fst-pop southern_africa.txt  --out ENA_SA_introns
vcftools --vcf dpgp_masked_introns.recode.vcf  --positions alllist --weir-fst-pop america.txt  --weir-fst-pop southern_africa.txt  --out NA_SA_introns 
vcftools --vcf dpgp_masked_introns.recode.vcf  --positions alllist --weir-fst-pop europe_north_africa.txt --weir-fst-pop america.txt   --out ENA_NA_introns 

## this part gets the overall pairwise fst for all short introns to generate the PBS tree.
vcftools --vcf dpgp_masked_introns.recode.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop southern_africa.txt  --out ENA_SA_introns_filtered  --positions alllist 
vcftools --vcf dpgp_masked_introns.recode.vcf --weir-fst-pop america.txt  --weir-fst-pop southern_africa.txt  --out NA_SA_introns_filtered  --positions alllist 
vcftools --vcf dpgp_masked_introns.recode.vcf --weir-fst-pop europe_north_africa.txt --weir-fst-pop america.txt   --out ENA_NA_introns_filtered  --positions alllist 
grep 'Weir and Cockerham mean Fst estimate' *_introns_filtered.log | sed 's/_introns_filtered.log:Weir and Cockerham mean Fst estimate://' >all_introns_fst.txt


## Lectin Ohana analyses

## subset vcf to only contain SNPs filtered for missingness, only sample 1 snp very 500bp for faster processing
vcftools --vcf dpgp_masked_3regions_0.05.recode.vcf --recode --out dpgp_masked_3regions_forindep2 --positions alllist --thin 500

## now convert vcf into file compatable for ohana
sed -e 's/2L/1/' -e 's/2R/2/' -e 's/3L/3/' -e 's/3R/4/' dpgp_masked_3regions_forindep2.recode.vcf  >dpgp_masked_3regions_forindep2.vcf 
## filter snps based on local LD
/scratch/Joel/software/plink2/v1.90b3.31/plink --vcf dpgp_masked_3regions_forindep2.vcf --indep 50 5 2 --recode 12 tab --out dpgp_masked_3regions_forindep2
/scratch/Joel/software/ohana/bin/convert ped2dgm dpgp_masked_3regions_forindep2.ped dpgp_masked_3regions_forindep2.dgm

## run structure analysis on intron snps from k=2 to k=6
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_forindep2.dgm -nt 3 -e 0.0001 -k 2 -qo dpgp_masked_3regions_forindep2_k2_Q.matrix -fo dpgp_masked_3regions_forindep2_k2_F.matrix &
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_forindep2.dgm -nt 3 -e 0.0001 -k 3 -qo dpgp_masked_3regions_forindep2_k3_Q.matrix -fo dpgp_masked_3regions_forindep2_k3_F.matrix &
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_forindep2.dgm -nt 3 -e 0.0001 -k 4 -qo dpgp_masked_3regions_forindep2_k4_Q.matrix -fo dpgp_masked_3regions_forindep2_k4_F.matrix &
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_forindep2.dgm -nt 3 -e 0.0001 -k 5 -qo dpgp_masked_3regions_forindep2_k5_Q.matrix -fo dpgp_masked_3regions_forindep2_k5_F.matrix &
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_forindep2.dgm -nt 3 -e 0.0001 -k 6 -qo dpgp_masked_3regions_forindep2_k6_Q.matrix -fo dpgp_masked_3regions_forindep2_k6_F.matrix &


## infer population variances and covariances 
/scratch/Joel/software/ohana/bin/nemeco dpgp_masked_3regions_forindep2.dgm dpgp_masked_3regions_forindep2_k2_F.matrix -e 0.0 -co dpgp_masked_3regions_forindep2_k2_C.matrix
/scratch/Joel/software/ohana/bin/nemeco dpgp_masked_3regions_forindep2.dgm dpgp_masked_3regions_forindep2_k3_F.matrix -e 0.0 -co dpgp_masked_3regions_forindep2_k3_C.matrix
/scratch/Joel/software/ohana/bin/nemeco dpgp_masked_3regions_forindep2.dgm dpgp_masked_3regions_forindep2_k4_F.matrix -e 0.0 -co dpgp_masked_3regions_forindep2_k4_C.matrix
/scratch/Joel/software/ohana/bin/nemeco dpgp_masked_3regions_forindep2.dgm dpgp_masked_3regions_forindep2_k5_F.matrix -e 0.0 -co dpgp_masked_3regions_forindep2_k5_C.matrix
/scratch/Joel/software/ohana/bin/nemeco dpgp_masked_3regions_forindep2.dgm dpgp_masked_3regions_forindep2_k6_F.matrix -e 0.0 -co dpgp_masked_3regions_forindep2_k6_C.matrix

## generate trees
/scratch/Joel/software/ohana/bin/convert cov2nwk dpgp_masked_3regions_forindep2_k2_C.matrix k2_forindep2_tree.nwk
/scratch/Joel/software/ohana/bin/convert cov2nwk dpgp_masked_3regions_forindep2_k3_C.matrix k3_forindep2_tree.nwk
/scratch/Joel/software/ohana/bin/convert cov2nwk dpgp_masked_3regions_forindep2_k4_C.matrix k4_forindep2_tree.nwk
/scratch/Joel/software/ohana/bin/convert cov2nwk dpgp_masked_3regions_forindep2_k5_C.matrix k5_forindep2_tree.nwk
/scratch/Joel/software/ohana/bin/convert cov2nwk dpgp_masked_3regions_forindep2_k6_C.matrix k6_forindep2_tree.nwk

## now working on main dataset
## get sites types in at least 40% samples per region
sed -e 's/2L/1/' -e 's/2R/2/' -e 's/3L/3/' -e 's/3R/4/'  alllist >alllist2
cat stop_codon.pos alllist2 > alllist3
Only analyses short intron variants and lectin stop codon
sed -e 's/2L/1/' -e 's/2R/2/' -e 's/3L/3/' -e 's/3R/4/'  -e 's/ /\t/g' introns.txt >introns2.txt
cat stop_codon.bed introns2.txt >introns3.txt
sort -k1 -k2 -V introns3.txt >introns4.txt
sed -e 's/2L/1/' -e 's/2R/2/' -e 's/3L/3/' -e 's/3R/4/' dpgp_masked_3regions_0.05.recode.vcf  >dpgp_masked_3regions_0.05.vcf 
## subset vcf to keep variants of interest
vcftools --vcf dpgp_masked_3regions_0.05.vcf --recode --out dpgp_masked_3regions_test2 --positions alllist3
vcftools --vcf dpgp_masked_3regions_test2.recode.vcf --recode --out dpgp_masked_3regions_test --bed introns4.txt 

## prepare dataset for ohana analyses
/scratch/Joel/software/plink2/v1.90b3.31/plink --vcf dpgp_masked_3regions_test.recode.vcf --recode 12 tab --out dpgp_masked_3regions_test
/scratch/Joel/software/ohana/bin/convert ped2dgm dpgp_masked_3regions_test.ped dpgp_masked_3regions_test.dgm 


## using estimated likelihood of the ancestral population proportions vector Q using k=2-6, infer the vector of admixture corrected allele frequencies F
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_test.dgm -nt 5 -qi dpgp_masked_3regions_forindep2_k2_Q.matrix -fo dpgp_masked_3regions_test_k2_F.matrix -e 0.0001 -fq 
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_test.dgm -nt 5 -qi dpgp_masked_3regions_forindep2_k3_Q.matrix -fo dpgp_masked_3regions_test_k3_F.matrix -e 0.0001 -fq 
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_test.dgm -nt 5 -qi dpgp_masked_3regions_forindep2_k4_Q.matrix -fo dpgp_masked_3regions_test_k4_F.matrix -e 0.0001 -fq 
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_test.dgm -nt 5 -qi dpgp_masked_3regions_forindep2_k5_Q.matrix -fo dpgp_masked_3regions_test_k5_F.matrix -e 0.0001 -fq 
/scratch/Joel/software/ohana/bin/qpas dpgp_masked_3regions_test.dgm -nt 5 -qi dpgp_masked_3regions_forindep2_k6_Q.matrix -fo dpgp_masked_3regions_test_k6_F.matrix -e 0.0001 -fq

## run global selection scan for covariance outliers by applying a likelihood model to each locus.
/scratch/Joel/software/ohana/bin/selscan dpgp_masked_3regions_test.dgm dpgp_masked_3regions_test_k2_F.matrix dpgp_masked_3regions_forindep2_k2_C.matrix  > scan_k2.txt
/scratch/Joel/software/ohana/bin/selscan dpgp_masked_3regions_test.dgm dpgp_masked_3regions_test_k3_F.matrix dpgp_masked_3regions_forindep2_k3_C.matrix  > scan_k3.txt
/scratch/Joel/software/ohana/bin/selscan dpgp_masked_3regions_test.dgm dpgp_masked_3regions_test_k4_F.matrix dpgp_masked_3regions_forindep2_k4_C.matrix  > scan_k4.txt
/scratch/Joel/software/ohana/bin/selscan dpgp_masked_3regions_test.dgm dpgp_masked_3regions_test_k5_F.matrix dpgp_masked_3regions_forindep2_k5_C.matrix  > scan_k5.txt
/scratch/Joel/software/ohana/bin/selscan dpgp_masked_3regions_test.dgm dpgp_masked_3regions_test_k6_F.matrix dpgp_masked_3regions_forindep2_k6_C.matrix  > scan_k6.txt

## generate c matrix for ancestry specific scans
cp dpgp_masked_3regions_forindep2_k3_C.matrix ans_specific_afr_k3
cp dpgp_masked_3regions_forindep2_k4_C.matrix ans_specific_afr_k4
cp dpgp_masked_3regions_forindep2_k5_C.matrix ans_specific_afr_k5
cp dpgp_masked_3regions_forindep2_k6_C.matrix ans_specific_afr_k6

## run african lineage selection scan for covariance outliers by applying a likelihood model to each locus.
/scratch/Joel/software/ohana/bin/selscan dpgp_masked_3regions_test.dgm dpgp_masked_3regions_test_k3_F.matrix dpgp_masked_3regions_forindep2_k3_C.matrix  -cs ans_specific_afr_k3 > scan_k3_afr.txt
/scratch/Joel/software/ohana/bin/selscan dpgp_masked_3regions_test.dgm dpgp_masked_3regions_test_k4_F.matrix dpgp_masked_3regions_forindep2_k4_C.matrix  -cs ans_specific_afr_k4 > scan_k4_afr.txt
/scratch/Joel/software/ohana/bin/selscan dpgp_masked_3regions_test.dgm dpgp_masked_3regions_test_k5_F.matrix dpgp_masked_3regions_forindep2_k5_C.matrix  -cs ans_specific_afr_k5 > scan_k5_afr.txt
/scratch/Joel/software/ohana/bin/selscan dpgp_masked_3regions_test.dgm dpgp_masked_3regions_test_k6_F.matrix dpgp_masked_3regions_forindep2_k6_C.matrix  -cs ans_specific_afr_k6 > scan_k6_afr.txt




