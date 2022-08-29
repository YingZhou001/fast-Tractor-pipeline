#ml Java


reftag=1kg
chr=X
cohorttag=cohort

for chr in 19 X
do
for cohorttag in cohort
do

dir=../tmp/${chr}

map=../data/plink.map/plink.chr${chr}.GRCh38.map
nthreads=8

job=O.${chr}.${cohorttag}


sbatch --ntasks 1 --cpus-per-task ${nthreads} --mem 40G -J ${job} -o ${job}.%j --wrap="
bash record.sh ${chr} ${reftag} ${cohorttag} ${dir} ${map} ${nthreads}
"

done
done
