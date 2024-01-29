
chr=$1 ##chromoson identifier
phe=$2 ##phenotype
cohorttag=$3 ## tag for the cohort
inpdir=$4 ##  place to store temporary ancestry information
outdir=$5 ## output tractor results

hapdosepref=${inpdir}/${cohorttag}.${chr}.flare-lai.anc
covar=../data/covars.for.tractor/${cohorttag}.${phe}.csv
out=${outdir}/${cohorttag}.${chr}.${phe}.tractor.out.txt

Rscript ../src.v1/run_tractor_yz.R \
--hapdose ${hapdosepref} \
--phe ${covar} \
--method logistic \
--out ${out}
