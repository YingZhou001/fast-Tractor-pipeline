## process the target cohort
## select samples of interest
## keep high quality sites with R2>0.8 and maf>0.5%
## remove multiallelic sites
## for information quality information not inlcuded in the vcf, the processing could be different

chr=$1  ##chromoson identifier
tag=$2 ##tag for data
dir=$3 ##place to store temporary data
sourcevcf=$4 ## vcf input of the cohort of interest

#chr=19
#datatag=share
#dir=/fh/scratch/delete90/kooperberg_c/yz/WHI-topmedimputed-data/test.data
#sourcevcf=../data/share.${chr}.vcf.gz

mkdir -p ${dir}
outvcf=${dir}/${tag}.${chr}.vcf.gz
snplist=${dir}/${tag}.${chr}.snplist


bcftools norm -m +any ${sourcevcf} | bcftools view --type snps -i 'INFO/R2>0.8' --min-alleles 2 --max-alleles 2 --min-af 0.005:minor | bcftools annotate -x ^FORMAT/GT,^INFO/R2 | bcftools norm --rm-dup all -Oz --output ${outvcf}

bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${outvcf} > ${snplist}
