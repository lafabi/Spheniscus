# Spheniscus
This repository include all scripts used in Spheniscus analyses



fastqc -o $PATH/ --noextract -t 2 -f fastq $PATH/${sample}*.fastq.gz

java -jar trimmomatic-0.39.jar PE -threads 10 -phred33 $PATH/${sample}*.fastq.gz $PATH/${sample}*.fastq.gz -baseout $PATH/${sample}*_trimmed.fq ILLUMINACLIP:$PATH/TruSeq2OVR-PE.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:25

bwa mem -t 20 -M -R '@RG\tID:sample\tLB:sample\tPL:ILLUMINA\tPU:A00887\tSM:sample' $PATH/index_Enovo_ref/GCA_010078495.1_BGI_Enov.V1_genomic.fna $PATH/${sample}*_trimmed_1P.fq $PATH/${sample}*_trimmed_2P.fq > $PATH/${sample}*.sam 

samtools view -q 20 -f 0x2 -bSh -@ 10 $PATH/${sample}*.sam > $PATH/${sample}*.bam
samtools sort -o $PATH/${sample}*_sorted.bam -@ 10 $PATH/${sample}*.bam
picard MarkDuplicates I=$PATH/${sample}*_sorted.bam O=$PATH/${sample}*_sorted_dedup.bam METRICS_FILE=$PATH/${sample}*_sorted_dedup.metrics.txt VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true CREATE_MD5_FILE=true TAGGING_POLICY=All ASSUME_SORT_ORDER=coordinate 


gatk3 -nt 6 -T RealignerTargetCreator -R $PATH/GCA_010078495.1_BGI_Enov.V1_genomic.fna -I $PATH/${sample}_sorted_dedup.bam -o $PATH/${sample}.realn.intervals -allowPotentiallyMisencodedQuals -S LENIENT

gatk3 -T IndelRealigner -R $PATH/index_Enovo_ref/GCA_010078495.1_BGI_Enov.V1_genomic.fna -I $PATH/${sample}_sorted_dedup.bam -targetIntervals $PATH/${sample}.realn.intervals -o $PATH/${sample}.realign.bam -allowPotentiallyMisencodedQuals -S LENIENT --generate_md5 

bcftools mpileup -Ou -f $PATH/index_Enovo_ref/GCA_010078495.1_BGI_Enov.V1_genomic.fna -b $PATH/bam.filelist -d 10000 -q 10 -Q 20 -a DP,SP  --threads 2 | bcftools call -vm  -Oz -f GQ  -o $PATH/scaffolds.vcf.gz --threads 2

BCFtools concat -f $PATH/concat.list --threads 10 -Ov -o $PATH/unfiltered.vcf

bcftools norm --check-ref w -f  $PATH/GCA_010078495.1_BGI_Enov.V1_genomic.fna  -o $PATH/unfiltered_norm.vcf.gz  -Oz --threads 8 $PATH/unfiltered.vcf

tabix -p vcf $PATH/unfiltered_norm.vcf.gz

#Data Set 1A

vcftools --gzvcf $PATH/unfiltered_norm.vcf.gz --minQ 30  --max-missing 1  --min-meanDP 2.5 --max-meanDP 7.5 --minDP 3 --recode --recode-INFO-all --out $PATH/filter.vcf


# Data Set 1B

plink --vcf $PATH/filter.vcf --recode --out neutral --allow-extra-chr

plink --file neutral --indep-pairwise 50 10 0.1 --allow-extra-chr

plink --file neutral --extract plink.prune.in

plink --vcf $PATH/filter.vcf --recode 12 --out $PATH/filter.ADMX --allow-extra-chr

admixture --cv $PATH/filter.ADMX.ped $K -j8 | tee log${K}.out

# Data Set 2A

vcftools --gzvcf $PATH/filter{banded_penguin_sp}.vcf --maf 0.05  --recode --recode-INFO-all --out $PATH/filter{banded_penguin_sp}_MAF.vcf

plink --vcf $PATH/filter{banded_penguin_sp}_MAF.vcf --recode --out neutral_{banded_penguin_sp} --allow-extra-chr

plink --file {banded_penguin_sp} --indep-pairwise 50 10 0.1 --allow-extra-chr

plink --file {banded_penguin_sp} --extract plink.prune.in

plink --vcf $PATH/filter{banded_penguin_sp}_MAF.vcf --recode 12 --out $PATH/filter_{banded_penguin_sp}.ADMX --allow-extra-chr

admixture --cv $PATH/filter_{banded_penguin_sp}.ped $K -j8 | tee log${K}.out


