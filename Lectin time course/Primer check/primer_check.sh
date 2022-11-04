wget http://dgrp2.gnets.ncsu.edu/data/website/dgrp2.vcf
vcftools --vcf dgrp2.vcf --indv line_437 --recode-INFO-all --out dgrp2_437 --recode
vcftools --vcf dgrp2.vcf --indv line_892 --recode-INFO-all --out dgrp2_892 --recode
bgzip dgrp2_437.recode.vcf&
bgzip dgrp2_892.recode.vcf&
tabix dgrp2_892.recode.vcf.gz
tabix dgrp2_437.recode.vcf.gz
samtools faidx dmel-all-chromosome-r5.13.fasta 2L | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 2R | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 3L | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 3R | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta X | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta Uextra | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 2RHet | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 2LHet | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 3LHet | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 3RHet | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta U | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta XHet | ../bcftools/bcftools consensus -M N -s line_437 -p line_437_ dgrp2_437.recode.vcf.gz >>dgrp2_437.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 2L | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 2R | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 3L | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 3R | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta X | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta Uextra | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 2RHet | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 2LHet | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 3LHet | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta 3RHet | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta U | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
samtools faidx dmel-all-chromosome-r5.13.fasta XHet | ../bcftools/bcftools consensus -M N -s line_892 -p line_892_ dgrp2_892.recode.vcf.gz >>dgrp2_892.fasta
makeblastdb -in dgrp2_437.fasta -out dgrp2_437 -parse_seqids -dbtype nucl
makeblastdb -in dgrp2_892.fasta -out dgrp2_892 -parse_seqids -dbtype nucl
blastn -db dgrp2_892 -query lectin_primers.fa -out primers_892.out.csv  -task blastn -outfmt 10
blastn -db dgrp2_437 -query lectin_primers.fa -out primers_437.out.csv  -task blastn -outfmt 10
grep '_INS\|_DEL' dgrp2_437.recode.vcf | grep '^2L' >dgrp2_437_indel_2L.txt
grep '_INS\|_DEL' dgrp2_892.recode.vcf | grep '^2L' >dgrp2_892_indel_2L.txt
