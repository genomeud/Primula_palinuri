# Copyright:	    Gabriele Magris & Fabio Marroni & Paloma 2024
# Aim:              Perform admixture analysis + draw plot 
# To add:		     
# Suggestions: 
# Fixes:  

spp=primula_palinuri
raw_files=${spp}/raw_reads

INPUT_DIR=${spp}/de_novo_pipeline/populations_stacks
# Folder containing the functions
FUNC_DIR=functions
mkdir -p ${INPUT_DIR}/admixture/input
cd ${INPUT_DIR}/admixture

PLINK_path=/software/plink_linux_x86_64
ADMIXTURE_path=/software/admixture_linux-1.23

# ----------------- #
# Create input file #
# ----------------- #

# create vcf for plink. We need letters when no using chromosomes with low number!
awk 'OFS="\t" {print "s_"$1,$0}' ${INPUT_DIR}/populations.snps.filtered.recode_MIS_filt.vcf  | sed 1d | cut -f2 --complement | cat <(head -1 ${INPUT_DIR}/populations.snps.filtered.recode_MIS_filt.vcf) - > ${INPUT_DIR}/populations.snps.filtered.recode_MIS_filt_PLINK.vcf

# need to transform the ped file obtained from stacks - first create bed file, because of an error in plink 
${PLINK_path}/plink --vcf ${INPUT_DIR}/populations.snps.filtered.recode_MIS_filt_PLINK.vcf  --maf 0.000000001 --allow-extra-chr --out ${INPUT_DIR}/admixture/input/populations.snps.filtered.recode_MIS_filt_PLINK --make-bed &> ${INPUT_DIR}/admixture/input/plink_bed.log

# transform the bed file to ped /map 
${PLINK_path}/plink --bfile ${INPUT_DIR}/admixture/input/populations.snps.filtered.recode_MIS_filt_PLINK --allow-extra-chr --recode 12 --out ${INPUT_DIR}/admixture/input/populations.snps.filtered.recode_MIS_filt_PLINK &> ${INPUT_DIR}/admixture/input/plink_ped.log

# sort the ped file in order to have the data separately for each population 
# assign to each individual the population name, sort the individuals and run admixture <- first column can be replaced with the family ID 
awk 'FNR==NR{a[$1]=$2;next}{ print a[$1],$0}' <(sed 1d /projects/marroni/seedforce/${spp}/raw_reads/sample_ID_and_populations.txt) ${INPUT_DIR}/admixture/input/populations.snps.filtered.recode_MIS_filt_PLINK.ped | \
# remove second column
cut -f 2 --complement -d " " | \
# sort the dataframe by population name 
sort -V -k1,1 > ${INPUT_DIR}/admixture/input/populations.snps.filtered.recode_MIS_filt_PLINK_sorted.ped

# --------- #
# Admixture #
# --------- #

# 1. run Admixture 
mkdir -p ${INPUT_DIR}/admixture/output/
cd ${INPUT_DIR}/admixture/output/
# optimal number of clusters K from 1 to 15
for K in {1..15}
do
    ${ADMIXTURE_path}/admixture -j24 ${INPUT_DIR}/admixture/input/populations.snps.filtered.recode_MIS_filt_PLINK_sorted.ped ${K} 
done &> ${INPUT_DIR}/admixture/output/admixture_K.log

# the output with Q will give the probability for each ind to belong to pop 1,2,3 depending on K. This will be used to print the graph!!!
# the output with P, for each SNP, it will calculate the allele freq in the ancestral populations (one prob for each K)


# 2. perform Evanno analysis + Cross-validation
mkdir -p ${INPUT_DIR}/admixture/output/evanno

# from 1 to 20 because is to see how many times the test is repeated!! is not the number of k
for a in {1..20}
do
    cd ${INPUT_DIR}/admixture/output/evanno
    mkdir -p ${a}
    cd ${a}
    for K in {1..15}
    do
        # create random numbers for the test
        b=$((RANDOM*200000+10))
        # cv = crossvalidation, it will do 10 times the same
        ${ADMIXTURE_path}/admixture -s ${b} --cv=10 -j24 ${INPUT_DIR}/admixture/input/populations.snps.filtered.recode_MIS_filt_PLINK_sorted.ped $K | tee log${K}.${a}_MIS_filt_PLINK_sorted.out
    done
done

# 3. collect results 
cd ${INPUT_DIR}/admixture/output/evanno
for run in {1..20}
do
    rm ${INPUT_DIR}/admixture/output/evanno/${run}_cv-Loglikelihood_MIS_filt_PLINK_sorted.txt
    for k in {1..15}
    do
        likelihood=`grep "Loglikelihood" ${run}/log${k}.${run}_MIS_filt_PLINK_sorted.out | tail -n 1 | cut -f 2 -d " "`
        cv=`grep CV ${run}/log${k}.${run}_MIS_filt_PLINK_sorted.out | cut -f 4 -d " "`
        echo -e $k"\t"$cv"\t"$likelihood >> ${INPUT_DIR}/admixture/output/evanno/${run}_cv-Loglikelihood_MIS_filt_PLINK_sorted.txt
    done
done


cd ${INPUT_DIR}/admixture/output/evanno/
paste *_cv-Loglikelihood_MIS_filt_PLINK_sorted.txt | sed 's/\./,/g' > ${INPUT_DIR}/admixture/output/evanno/all_cv-Loglikelihood_MIS_filt_PLINK_sorted.txt

# calculate mean CV error per K and prepare summary table. test optimal number of clusters K from 1 to 15
awk 'OFS="\t" {print $1,$2,$5,$8,$11,$14,$17,$20,$23,$26,$29,$32,$35,$38,$41,$44,$47,$50,$53,$56,$59}' ${INPUT_DIR}/admixture/output/evanno/all_cv-Loglikelihood_MIS_filt_PLINK_sorted.txt | sed 's/,/\./g' | awk 'OFS="\t" {out=$1; for(i=2; i<=NF;i++) sum+=$i; print out,sum/(NF-1); sum=0 }' | \
cat <(echo -e "K\tCV_profile") - > ${INPUT_DIR}/admixture/output/populations.snps.filtered.recode_MIS_filt_evanno.txt


# 4. Draw CV plot - from the summary table. Check the cross validation output to choose the correct number of K where the cv error reach the "lowest" value. This step is subject to interpretation
Rscript ${FUNC_DIR}/a09_admixture.r \
-o ${INPUT_DIR}/admixture/output/ \
--prefix populations.snps.filtered.recode_MIS_filt \
--infile ${INPUT_DIR}/admixture/output/populations.snps.filtered.recode_MIS_filt_evanno.txt \
--step draw_cross-validation

# 5. draw admixture plot 
# set optimal number of clusters according to the results obtained with cross validation (k=optimal number of clusters=
k=3

Rscript ${FUNC_DIR}/a09_admixture.r \
-i ${INPUT_DIR}/admixture/input/ \
-o ${INPUT_DIR}/admixture/output/ \
--prefix populations.snps.filtered.recode_MIS_filt_PLINK_sorted \
-k ${k}.Q \
-P ${raw_files}/sample_ID_and_populations.txt \
--step draw_admixture
