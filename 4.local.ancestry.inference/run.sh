#ml Java


reftag=1kg
refpanel=../data/sample.anc.map.txt
nthreads=8

for chr in 19 X
do
for cohorttag in cohort
do

dir=../tmp/${chr}
map=../data/plink.map/plink.chr${chr}.GRCh38.map

job=O.${chr}.${cohorttag}

sbatch --ntasks 1 --cpus-per-task ${nthreads} --mem 40G -J ${job} -o ${job}.%j --wrap="
bash record.sh ${chr} ${reftag} ${cohorttag} ${dir} ${map} ${refpanel} ${nthreads}
"
done
done
