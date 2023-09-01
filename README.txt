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

for m in {1..10}; do
    for i in {1..10}; do
        # Generate random seed
        s=$RANDOM
        treemix -i banded_penguin.txt.gz -o banded_penguin.${i}.${m} -global -m ${m} -k 500 -seed ${s}
    done
done

tar -zcvf /banded_peguin_dir.tar.gz

https://rfitak.shinyapps.io/OptM/

bed2diffs_v1 --bfile /banded_penguin_pop --nthreads 2 
runeems_snps --params banded_penguin_pop.ini

#Data Set 3A

samtools mpileup -C50 -uf $REF $PATH/${sample}.realign.bam  | bcftools view -c - | vcfutils.pl vcf2fq -d 15  -D 45 | gzip > $PATH/diploid_${sample}.fq.gz
fq2psmcfa -q20 $PATH/diploid_${sample}.fq.gz > $PATH/${sample}.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $PATH/${sample}.psmc  $PATH/${sample}.psmcfa

#Data Set 3B

samtools mpileup -C50 -uf $REF $PATH/${sample}.realign.bam  | bcftools view -c - | vcfutils.pl vcf2fq -d 3  -D 14 | gzip > $PATH/diploid_${sample}.fq.gz
fq2psmcfa -q20 $PATH/diploid_${sample}.fq.gz > $PATH/${sample}.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $PATH/${sample}.psmc  $PATH/${sample}.psmcfa

# Data Set4

vcftools --gzvcf $PATH/filter{banded_penguin_sp}.vcf --maf 0.05  --recode --recode-INFO-all --out $PATH/filter{banded_penguin_sp}_MAF.vcf

plink --vcf $PATH/filter{banded_penguin_sp}_MAF.vcf --recode --out neutral_{banded_penguin_sp} --allow-extra-chr

plink --file {banded_penguin_sp} --indep-pairwise 50 10 0.1 --allow-extra-chr

plink --file {banded_penguin_sp} --extract plink.prune.in

~/angsd/angsd -gl 1 \
-anc $PATH/GCA_010078495.1_BGI_Enov.V1_genomic.fna \
-ref $PATH/GCA_010078495.1_BGI_Enov.V1_genomic.fna \
-bam $PATH/banded_penguin_pop.filelist \
-rf $PATH/NoCDS.angsd.regions \
-out $PATH/POP \
-dosaf 1 \
-baq 1 \
-C 50 \
-minMapQ 30 \
-minQ 20 \
-P 20

java -cp stairway_plot_es Stairbuilder SP_blueprint
bash SP_blueprint.sh


bcftools consensus -s ${sample} -f $REF $VCF -o $FASTA/${sample}_consenso.fa

gffread Little_Blue_penguin.gff  -g $REF -x banded_penguin_sp -C  -V -H  -J

phyluce_probe_run_multiple_lastzs_sqlite --db test.sqlite --output $PATH/UCE/test-lastz --scaffoldlist list --genome-base-path . --probefile $PATH/UCE/uce-5k-probes.fasta --cores 8

phyluce_probe_slice_sequence_from_genomes --lastz $PATH/UCE/test-lastz --conf $PATH/UCE/genomes.conf --flank 750 --name-pattern "uce-5k-probes.fasta_v_{}.lastz.clean" --output $PATH/UCE/test_fasta

mafft banded_penguins_sequences.fasta > aligned_sequences.fasta

iqtree -s banded_penguins_UCE_alignment.fasta -m GTR+G4 -bb 1000 -nt 4 -o banded_penguin_UCE_tree.nwk

iqtree -s banded_penguins_CDS_alignment.fasta -m GTR+G4 -bb 1000 -nt 4 -o banded_penguin_UCE_tree.nwk

iqtree -s banded_penguins_exon_alignment.fasta -m GTR+G4 -bb 1000 -nt 4 -o banded_penguin_UCE_tree.nwk








