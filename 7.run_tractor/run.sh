ml fhR

for chr in 19 X
do
for cohorttag in cohort
do

inpdir=../tmp/${chr}
outdir=../out/${cohorttag}.ANC.tractor.out

mkdir -p ${outdir}

for phe in CAD MI
do

job=O.${chr}.${cohorttag}.anc${phe}

#sbatch --mem 40G -J ${job} -o ${job}.%j --wrap="
bash record.sh ${chr} ${phe} ${cohorttag} ${inpdir} ${outdir}
#"

done
done
done
