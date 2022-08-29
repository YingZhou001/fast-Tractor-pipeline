##reference based imputation

beagle=/fh/fast/kooperberg_c/yz/tools/pub/beagle/beagle.22Jul22.46e.jar


chr=$1 ##chromoson identifier
reftag=$2 ##tag for reference panel, which is 1000 genome
cohorttag=$3 ## tag for the cohort
dir=$4 ## place to store temporary data
map=$5 ## recombination map in plink format
nthreads=$6 ## number of threads used for haplotype phasing

#chr=19
#dir=/fh/scratch/delete90/kooperberg_c/yz/WHI-topmedimputed-data/test.data
#map=../data/plink.chr${chr}.GRCh38.map


ref=${dir}/${cohorttag}.${chr}.vcf.gz
gt=${dir}/${reftag}.${chr}.vcf.gz
out=${dir}/${reftag}.${chr}.imputed

java -jar ${beagle} gt=${gt} ref=${ref} map=${map} out=${out} nthreads=${nthreads}
