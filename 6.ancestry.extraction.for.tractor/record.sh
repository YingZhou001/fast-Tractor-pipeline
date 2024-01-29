formatter=../src.v1/tractor_anc_formatter/target/release/tractor_anc_formatter

chr=$1 ##chromoson identifier
cohorttag=$2 ## tag for the cohort
dir=$3 ##  place to store temporary data

ancvcf=${dir}/${cohorttag}.${chr}.flare-lai.anc.vcf.gz
outpref=${dir}/${cohorttag}.${chr}.flare-lai

${formatter} ${ancvcf} ${outpref}
