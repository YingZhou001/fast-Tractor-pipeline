

extractAncestry=../src.v1/vcf_anc_formatter/vcf_anc_formatter



chr=$1 ##chromoson identifier
anc=$2 ##ancestry code, can be found in flare's output vcf
cohorttag=$3 ## tag for the cohort
dir=$4 ##  place to store temporary data



#chr=19
#reftag=1kg
#cohorttag=share
#dir=/fh/scratch/delete90/kooperberg_c/yz/WHI-topmedimputed-data/test.data

ancvcf=${dir}/${cohorttag}.${chr}.flare-lai.anc.vcf.gz
outpref=${dir}/${cohorttag}.${chr}.flare-lai

#rm ${out}.*

${extractAncestry} ${ancvcf} ${anc} ${outpref}
