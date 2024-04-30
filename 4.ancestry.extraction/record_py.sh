
extractAncestry=../src.v1/vcf_anc_formatter.py

chr=$1 ##chromoson identifier
cohorttag=$2 ## tag for the cohort
dir=$3 ##  place to store temporary data



#chr=19
#reftag=1kg
#cohorttag=share
#dir=/fh/scratch/delete90/kooperberg_c/yz/WHI-topmedimputed-data/test.data

ancvcf=${dir}/${cohorttag}.${chr}.flare-lai.anc.vcf.gz
outpref=${dir}/${cohorttag}.${chr}.flare-lai


zcat ${ancvcf} | ${extractAncestry} ${outpref}
