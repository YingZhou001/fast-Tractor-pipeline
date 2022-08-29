#tractor's gwas

chr=$1 ##chromoson identifier
anc=$2 ## 
cohorttag=$3 ## tag for the cohort
datadir=$4 ## place to store temporary data
pheno=$5 ##phenotype name
covarnames=$6 ## a string of covariates names, sepearated by comma
phenofile=$7 ## file for phenofiles in plink format
covarfile=$8 ## file for covariates in plink format
outdir=$9 ## output directory for plink
final=${10} ## output directory for secondary analysis


#chr=19
#datadir=/fh/scratch/delete90/kooperberg_c/yz/WHI-topmedimputed-data/test.data
#cohorttag=share
#pheno=MI
#covarnames=AGE,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10
#phenofile=../data/covars.for.plink/share.${pheno}.csv
#covarfile=../data/covars.for.plink/share.${pheno}.csv
#outdir=/fh/scratch/delete90/kooperberg_c/yz/WHI-topmedimputed-data/test.data.out


vcfpref=${datadir}/${cohorttag}.${chr}.flare-lai
outpref1=${outdir}/${cohorttag}.${pheno}.${chr}.allele
outpref2=${outdir}/${cohorttag}.${pheno}.${chr}.hapanc

## Tractor analysis

plink --vcf ${vcfpref}.anc${anc}.vcf.gz --vcf-half-call haploid --allow-no-sex --pheno ${phenofile} --pheno-name ${pheno} --covar ${covarfile} --covar-name  ${covarnames} --logistic beta hide-covar --out ${outpref1}.anc${anc} --ci 0.95 --freq case-control

paste -d' ' ${outpref1}.anc${anc}.assoc.logistic ${outpref1}.anc${anc}.frq.cc | sed "s/ \+/ /g" | sed "s/^ //g" | cut -f1-12,17-20 -d' ' | grep -v NA | gzip -c > ${outpref1}.anc${anc}.assoc.simple.gz

mv ${outpref1}.anc${anc}.assoc.simple.gz ${final}/

## admixture mapping

plink --vcf ${vcfpref}.hapanc${anc}.vcf.gz --allow-no-sex --pheno ${phenofile} --pheno-name ${pheno} --covar ${covarfile} --covar-name  ${covarnames} --logistic beta hide-covar --out ${outpref2}.anc${anc} --ci 0.95 --freq case-control

paste -d' ' ${outpref2}.anc${anc}.assoc.logistic ${outpref2}.anc${anc}.frq.cc | sed "s/ \+/ /g" | sed "s/^ //g" | cut -f1-12,17-20 -d' ' | grep -v NA | gzip -c > ${outpref2}.anc${anc}.admap.simple.gz

mv ${outpref2}.anc${anc}.admap.simple.gz ${final}/
