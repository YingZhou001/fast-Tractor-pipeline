
reftag=1kg

for chr in 19 X
do
for cohorttag in cohort
do

dir=../tmp/${chr}
refvcf=../data/1000genome/1000genome.${chr}.vcf.gz

samples=../data/1000genome/small.sample.txt
snplist=${dir}/${cohorttag}.${chr}.snplist

job=O.${chr}.${cohorttag}

sbatch --mem 40G -J ${job} -o ${job}.%j --wrap="
bash record.sh ${chr} ${reftag} ${dir} ${refvcf} ${samples} ${snplist}
"

done
done
