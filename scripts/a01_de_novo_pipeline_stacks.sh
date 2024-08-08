# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              De novo pipeline for building loci from ddRAD paired end experiments using STACKS + filtering of genotyping data to ensure high quality + genetic analyses including MAF, mean expected and observed heterozygosity, FST, IBD, PCA
# To add:		     
# Suggestions: 
# Fixes:  

spp=primula_palinuri

# create output folders
mkdir -p ${INPUT_DIR}/${spp}/de_novo_pipeline/populations_stacks

out1=${INPUT_DIR}/${spp}/de_novo_pipeline
out2=${INPUT_DIR}/${spp}/de_novo_pipeline/populations_stacks
# folder containing fastq sequences
raw_files=${INPUT_DIR}/${spp}/raw_reads

# Required population file inside raw_files folder 
# population_names.txt example:
# ID1 pop1
# ID2 pop1
# ID3 pop2
# ID4 pop2

###############################
# Run STACKS de novo pipeline #
###############################

perl ${STACKS_path}/denovo_map.pl -M 2 -T 24 -o ${out1} --popmap ${raw_files}/population_names.txt --samples ${raw_files} --paired &> ${out1}/de_novo_pipeline.log

# Run populations. At least 75% of the individuals with SNP info
${STACKS_path}/populations -P ${out1} -O ${out2} -M ${raw_files}/population_names.txt -t 24 --fstats --vcf -R 0.75 --max-obs-het 0.8 &> ${out2}/populations_stacks_pipeline.log

# Filter vcf file for individuals with at least 5 reads
vcftools --vcf ${out2}/populations.snps.vcf --minDP 5 --recode --out ${out2}/populations.snps.filtered &> ${out2}/vcf_filter.log

# Get allele count table from filtered vcf considering all the individuals
Rscript ${FUNC_DIR}/a01_allele_frequencies_individuals.r \
-S "ID1,ID2,ID3,ID4" \
-V ${out2}/populations.snps.filtered.recode.vcf \
-O ${out2}/allele_info.tbl &> ${out2}/allele_info.log

# Create column with minor allele considering all individuals together
sed 1d ${out2}/allele_info.tbl | \
awk 'OFS="\t" {if($NF<($(NF-2)+$(NF-1))) print $0}' | awk 'OFS="\t" {if($10<$11)print $0,$10;if($11<$10)print $0,$11;if($10==$11)print $0,$10}' | cat <(paste <(head -1 ${out2}/allele_info.tbl) <( echo -e "MinAll")) - > ${out2}/allele_info_MinAll.tbl

# Calculate minor allele frequency, expected heterozygosity and histogram with minor allele frequency
Rscript ${FUNC_DIR}/a02_HWE_analysis.r \
-F ${out2}/allele_info_MinAll.tbl \
-O ${out2}/allele_freq_and_exp_het.tbl \
-P ${out2}/histogram_MAF.jpeg


# Create column with pop name = ALL. For some plots we will include all the individuals together as a unique sample
sed 1d ${out2}/allele_freq_and_exp_het.tbl | \
awk 'OFS="\t" {print $0,"ALL"}' | \
cat <(paste <(head -1 ${out2}/allele_freq_and_exp_het.tbl) <(echo -e "population")) - > ${out2}/allele_freq_and_exp_het_All.tbl

# Filter positions prior to divide by populations. We only want positions that were genotyped with at least 5 reads in 50% of the individuals. file without vcf header "##"
join -a 1 -1 1 -2 1 <(awk 'OFS="\t" {print $1"|"$2}' ${out2}/allele_freq_and_exp_het.tbl | sed 1d | sort -k1,1 ) <(grep -v "##" ${out2}/populations.snps.filtered.recode.vcf | sed 1d | awk 'OFS="\t" {print $1"|"$2,$0}' | sort -k1,1) -t $'\t' | cut -f1 --complement | sort -n -k1,2 | \
cat <(grep -v "##" ${out2}/populations.snps.filtered.recode.vcf | head -1) - > ${out2}/populations.snps.filtered.recode_MIS_filt.vcf

# keep normal vcf format for other analyses
cat <(grep "##" ${out2}/populations.snps.filtered.recode.vcf) ${out2}/populations.snps.filtered.recode_MIS_filt.vcf > ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf

# Divide the vcf by populations. This will generate a separate file for each population + histogram with minor allele frequency for every population.
# output file for each population: the number of reference alleles (POP1_REF); the number of alternative alleles (POP1_ALT); the number of missing alleles in the population (POP1_MIS); the total number of alleles for that population (POP1_REF + POP1_ALT + POP1_MIS) (POP1_N); expected heterozygosity (POP1_He); count for the minor allele (POP1_minALL) and the minor allele frequency (POP1_minALLfreq); population ID (Popname)
# otput format example
# #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT          POP1_REF  POP1_ALT  POP1_MIS  POP1_N  POP1_He  POP1_minALL  POP1_minALLfreq     Popname
# 2       215  .   C    T    .     PASS    .     GT:DP:AD:GQ:GL  23        15        2         40      0.43125  15           0.394736842105263   POP1
# 2       30   .   G    T    .     PASS    .     GT:DP:AD:GQ:GL  38        0         2         40      0        0            0                   POP1
# 2       96   .   C    T    .     PASS    .     GT:DP:AD:GQ:GL  38        0         2         40      0        0            0                   POP1

cd ${out2}
Rscript ${FUNC_DIR}/a03_divide_populations_and_HWE_analysis.r \
-S "ID1,ID2,ID3,ID4" \
-V ${out2}/populations.snps.filtered.recode_MIS_filt.vcf \
-N 0 \
-F ${out2}/histogram_MAF_pops_ALL.pdf &> ${out2}/divide_populations_and_HWE.log &


# perform PCA analysis and IBD with SNPRelate
for step in pca_analysis ibd_analysis 
do
    Rscript ${FUNC_DIR}/a04_SNPRelate_PCA_IBD.r \
    -p ${spp} \
    -V ${out2}/populations.snps.filtered.recode_MIS_filt.vcf \
    -P ${raw_files}/sample_ID_and_populations.txt \
    -O ${out2}/SNPRelate/ \
    -s ${step} &> ${out2}/${step}.log
done

# FST heatmap. We need to adapt the file that was generated with stacks populations.fst_summary.tsv for the heatmap.
# # input table
# CIM  PPA        PC        ID        FIU       SGP        LAM
# CIM  0.192775   0.18101   0.182763  0.193849  0.0837179  0.181239
# PPA  0.240124   0.246545  0.251917  0.229188  0.0884559
# PC   0.15312    0.160193  0.219462  0.243646
# ID   0.0900836  0.225295  0.25797
# FIU  0.232218   0.260227
# SGP  0.230969

# # desired output table
# 0    CIM  PPA       PC        ID        FIU        SGP        LAM
# CIM  0    0.192775  0.18101   0.182763  0.193849   0.0837179  0.181239
# PPA  0    0         0.240124  0.246545  0.251917   0.229188   0.0884559
# PC   0    0         0         0.15312   0.160193   0.219462   0.243646
# ID   0    0         0         0         0.0900836  0.225295   0.25797
# FIU  0    0         0         0         0          0.232218   0.260227
# SGP  0    0         0         0         0          0          0.230969
# LAM  0    0         0         0         0          0          0

# prepare input table
awk -F'\t' -v OFS='\t' '{ for (i = 1; i <= NF; ++i) sub(/^$/, 0, $i) } 1' ${out2}/populations.fst_summary.tsv | cat - <(echo -e "LAM\t0\t0\t0\t0\t0\t0\t0") > ${out2}/populations.fst_summary_for_heatmap.tsv

# FST heatmap
Rscript ${FUNC_DIR}/a05_FST_heatmap.r \
-I ${out2}/populations.fst_summary_for_heatmap.tsv \
-O ${out2}/ &> ${out2}/FST.log

# plots
# boxplot with expected and observed heterozygosity and histograms with mean exp and obs heterozygosity
Rscript ${FUNC_DIR}/a06_plots_heterozygosity.r 

# Count heterozygous SNPs + scatterplot
Rscript ${FUNC_DIR}/a07_count_het_SNPs_per_ind.r \
-V ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf \
-O ${out2}/genotype_counts.tbl \
-P ${out2}/ratio_heterozygosity.jpeg

# Create phylogenetic trees
Rscript ${FUNC_DIR}/a08_phyl_tree_plot.r \
-V ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf \
-M ${raw_files}/sample_ID_and_populations.txt \
-O ${out2}/phyl_tree_NJ

# prepare files for procrustes analysis
# prepare genetic coordinates (using the PCA generated with SNPrealate)
# We will also need geographic coordinates
cut -f1-3 ${out2}/SNPRelate/PCA.all_snp.${spp}.txt | sed 's/sample\.id/ID/g' > ${out2}/genetic_coordinates.txt
# We will also need geographic coordinates (name of the file = coordinates_individuals.txt) as the following example:
# ID     Latitude   Longitude
# ID1    40.032339  15.272951
# ID2    40.032329  15.272896
# ID3    40.032261  15.272726
# ID4    40.032179  15.272605


Rscript ${FUNC_DIR}/a10_proclustes.r \
-C ${out2}/coordinates_individuals.txt \
-G ${out2}/genetic_coordinates.txt \
-M ${raw_files}/sample_ID_and_populations.txt

# "test" directory to verify that everything works properly
spp=primula_palinuri
INPUT_DIR=/projects/marroni/seedforce/primula_palinuri
FUNC_DIR=/projects/marroni/seedforce/primula_palinuri/github/functions
out1=/projects/marroni/seedforce/primula_palinuri/github/tests
out2=/projects/marroni/seedforce/primula_palinuri/github/tests
raw_files=${INPUT_DIR}/raw_reads
