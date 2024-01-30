
## remove sites from the reference genome that not exists in target population.


chr=$1 ## chromoson identifier
tag=$2 ## tag for reference panel
dir=$3 ## place to store temporary data, should be the same from the first step
sourcevcf=$4 ##  vcf input for the 1000 genome reference panel
samples=$5 ## ids for the ancestry populations, one id per line
cohortsnplist=$6 ## variants that exist in the target cohort, from the first step

#chr=19
#tag=1kg
#dir=/fh/scratch/delete90/kooperberg_c/yz/WHI-topmedimputed-data/test.data
#sourcevcf=../data/1kg.${chr}.vcf.gz
#samples=../data/YRI_CEU.txt
#cohortsnplist=${dir}/share.${chr}.snplist


outvcf=${dir}/${tag}.${chr}.vcf.gz
snplist=${dir}/${tag}.${chr}.snplist
tmpsnplist=${dir}/tmp.${tag}.${chr}.snplist


##join biallelic sites into multiallelic records
bcftools norm -m +any ${sourcevcf} | bcftools query -f '%CHROM:%POS:%REF:%ALT\n' > ${tmpsnplist}

##find out the shared sites with the targeted cohort
grep -f ${tmpsnplist} ${cohortsnplist} | cut -f1,2 -d":" | sed "s/:/\t/g"  > ${snplist}
#### alternative solution for finding the overlapping variants if the grep command is too slow
# cat ${tmpsnplist} ${cohortsnplist} | sort | uniq -D | cut -f1,2 -d":" | sed "s/:/\t/g" | sort -k2 -n | uniq > ${snplist}
rm ${tmpsnplist}

##filtration
bcftools view --types snps -S ${samples} -R ${snplist} ${sourcevcf} | bcftools annotate -x ^FORMAT/GT,INFO --force -Oz --output ${outvcf}
