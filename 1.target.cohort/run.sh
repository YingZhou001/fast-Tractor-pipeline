

for chr in 19 X
do
for cohorttag in cohort
do

dir=../tmp/${chr}
cohortvcf=../data/cohort/cohort.${chr}.vcf.gz

job=O.${chr}.${cohorttag}

sbatch --mem 40G -J ${job} -o ${job}.%j --wrap="
bash record.sh ${chr} ${cohorttag} ${dir} ${cohortvcf}
"


done
done
