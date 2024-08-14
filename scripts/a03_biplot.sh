# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              prepare table for biplot PCA. Run Biplot function for all SNPs present in all the individuals and subset of 100,80,60,50,40 and 20 top SNPs more contributing to population differentiation
# To add:		     
# Suggestions: 
# Fixes: 


spp=primula_palinuri

out1=${INPUT_DIR}/${spp}/de_novo_pipeline
out2=${INPUT_DIR}/${spp}/de_novo_pipeline/populations_stacks
# folder containing fastq sequences
raw_files=${INPUT_DIR}/${spp}/raw_reads
# Folder containing the functions
FUNC_DIR=functions

# create output folder
mkdir -p ${out2}/biplot
# Desired input --> each column represents a SNP, and each raw represents a single individual
grep  -v "##" ${out2}/populations.snps.filtered.recode_MIS_filt.vcf | \
# transform the allele info into 0,1 and 2 code. 0 = REF/REF; 1 = ALT/REF or REF/ALT; 2 = ALT/ALT
awk 'OFS="\t" { for (i=10; i<=NF; i++) { if ($i ~ /0\/0/) $i = 0 ; if ($i ~ /1\/0/) $i = 1 ; if ($i ~ /0\/1/) $i = 1; if ($i ~ /1\/1/) $i = 2; if ($i ~ /.\/./) $i = "NA"} } 1' | \
# transpose columns in raws
awk 'OFS="\t" {print $1"_"$2,$0}' | \
cut -f 2,3,4,5,6,7,8,9,10 --complement | \
awk 'OFS="\t" {
       for (f = 1; f <= NF; f++) { a[NR, f] = $f } 
     }
     NF > nf { nf = NF }
     END {
       for (f = 1; f <= nf; f++) {
           for (r = 1; r <= NR; r++) {
               printf a[r, f] (r==NR ? RS : FS)
           }
       }
    }' | sed 's/\s/\t/g' | sed 's/#CHROM_POS/CHROM_POS/g' > ${out2}/biplot/transposed_ind_SNPs.tbl

# table with only SNPs present in all the individuals
awk '
BEGIN { FS = "\t"; OFS = "\t" }
{
    for (i = 1; i <= NF; i++) {
        if ($i == "NA") {
            na_count[i]++
        }
        data[NR,i] = $i
    }
}
END {
    for (i = 1; i <= NF; i++) {
        if (i == NF) {
            printf("%d\n", na_count[i])
        } else {
            printf("%d\t", na_count[i])
        }
    }
    for (row = 1; row <= NR; row++) {
        for (col = 1; col <= NF; col++) {
            if (col == NF) {
                printf("%s\n", data[row,col])
            } else {
                printf("%s\t", data[row,col])
            }
        }
    }
}' ${out2}/biplot/transposed_ind_SNPs.tbl | awk 'BEGIN { FS = "\t" } NR == 1 { for (i = 1; i <= NF; i++) { if ($i == 0) { col[i] = 1 } } } { for (i = 1; i <= NF; i++) { if (col[i]) { printf "%s\t", $i } } printf "\n" }' | sed 1d > ${out2}/biplot/transposed_ind_SNPs_no_NA.tbl

# run biplot for all SNPs and subset
Rscript ${FUNC_DIR}/a11_biplot.r \
-I ${out2}/biplot/transposed_ind_SNPs_no_NA.tbl \
-M ${raw_files}/sample_ID_and_populations.txt \
-O ${out2}/biplot/
