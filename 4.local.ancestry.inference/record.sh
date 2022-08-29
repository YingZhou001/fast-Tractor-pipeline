#run local ancestry inference

flare=/home/yzhou3/fast/tools/pub/flare/flare.jar

chr=$1 ##chromoson identifier
reftag=$2 ##tag for reference panel, which is 1000 genome
cohorttag=$3 ## tag for the cohort
dir=$4 ##  place to store temporary data
map=$5 ## recombination map in plink format
refpanel=$6 ## ancestry information for each sample, in the format of "sample id"+"tab"+"ancestry"
nthreads=$7 ## number of threads to run local ancestry inference



#chr=19
#reftag=1kg
#cohorttag=share
#dir=/fh/scratch/delete90/kooperberg_c/yz/WHI-topmedimputed-data/test.data
#refpanel=../data/YRI_CEU.anc.map.txt
#map=../data/plink.chr${chr}.GRCh38.map
#nthread=5


refvcf=${dir}/${reftag}.${chr}.imputed.vcf.gz
gt=${dir}/${cohorttag}.${chr}.vcf.gz
out=${dir}/${cohorttag}.${chr}.flare-lai


java -jar ${flare} ref=${refvcf} ref-panel=${refpanel} gt=${gt} map=${map} out=${out} nthreads=${nthreads}
bcftools index -t -f ${out}.anc.vcf.gz
