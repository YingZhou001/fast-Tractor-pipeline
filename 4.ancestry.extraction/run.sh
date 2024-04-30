

for chr in 19 X
do
for cohorttag in cohort
do

dir=../tmp/${chr}

if false
then
# for python formating
## slow!!!
job=O.${chr}.${cohorttag}.anc${anc}
sbatch --mem 40G -J ${job} -o ${job}.%j --wrap="
bash record_py.sh ${chr} ${cohorttag} ${dir}
"
fi

if true
then
# for C formating
for anc in 0 1
do

job=O.${chr}.${cohorttag}.anc${anc}

sbatch --mem 40G -J ${job} -o ${job}.%j --wrap="
bash record.sh ${chr} ${anc} ${cohorttag} ${dir}
"
fi

done
done
done
