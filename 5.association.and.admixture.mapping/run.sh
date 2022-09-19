
covarnames=AGE,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10

for chr in 19 X
do
for cohorttag in cohort
do
dir=../tmp/${chr}
plinkoutdir=../out/${cohorttag}.ANC.plink.out
finaloutdir=../out/${cohorttag}.ANC.final.out
mkdir -p ${plinkoutdir} ${finaloutdir}

for pheno in MI CAD
do
phenofile=../data/covars.for.plink/${cohorttag}.${pheno}.csv
covarfile=../data/covars.for.plink/${cohorttag}.${pheno}.csv

for anc in 0 1
do

job=O.${chr}.${cohorttag}.${pheno}

sbatch --mem 40G -J ${job} -o ${job}.%j --wrap="
bash record.sh ${chr} ${anc} ${cohorttag} ${dir} ${pheno} ${covarnames} ${phenofile} ${covarfile} ${plinkoutdir} ${finaloutdir}
"
done
done
done
done
