

for chr in 19 X
do
for cohorttag in cohort
do

dir=../tmp/${chr}

for anc in 0 1
do

job=O.${chr}.${cohorttag}.anc${anc}

sbatch --mem 40G -J ${job} -o ${job}.%j --wrap="
bash record.sh ${chr} ${anc} ${cohorttag} ${dir}
"

done
done
done
