# About

This is a pipeline aimed to conduct admixture mapping and simplify the Tractor analysis [link](https://github.com/Atkinson-Lab/Tractor/wiki). 
This pipeline can handle admixed population that has multiple ancestries.

This pipeline is developed under linux/unix environment, and it can be easily adapted to high-performance-computing server for massive computation. Any questions or complaints are welcomed to post in the issue section or emailed to yzhou3 at fredhutch.org .

Updated Date: 9/19/2022

# Requirement

To run this pipeline, you need to know the basic knowledge of linux/unix and some common bioinformatic tools list in the bellow section.

## Computational tools

* bcftools v1.11 or later [link](https://samtools.github.io/bcftools/)
* htslib v1.15.1 or later [link](https://github.com/samtools/htslib)
* plink 1.9 [link](https://zzz.bwh.harvard.edu/plink/)
* beagle 5.4 [link](https://faculty.washington.edu/browning/beagle/beagle.html)
* flare [link](https://github.com/browning-lab/flare)
* extractAncestry [included in this pipeline]

## Data

* 1000 genome high coverage WGS data [link](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
* Recombination map in plink format [link](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)

# How to use this pipeline

1. Download all folders to your local disk and install all required software.
2. Follow the procedures in the order of folders' names.
3. In each folder, the "record.sh" records function scripts, and the "run.sh" specifies input, output and also job submission command. Feel free to modify them based on your own understanding.

We included an example data set, so you can just run the example without changing any input/output parameters.

# Folders & Procedures

## 1.target.cohort

* requirement: bcftools

In this folder, we will process the cohort genotype data in vcf format. 
Assume the genotype data is well phased and imputed, we will process the genotype data and only keep binary sites with imputation $R^2$ larger than 0.8 and minimum MAF larger than 0.5%. Duplicated sites are removing by keeping the first occurrence in vcf record.

If the imputation score $R^2$ is not in the vcf file, you need to remove "-i 'INFO/R2>0.8'" and ",^INFO/R2" from the following command and use other tools to create filtered genotype data set.

```bash
bcftools norm -m +any ${sourcevcf} \
| bcftools view --type snps -i 'INFO/R2>0.8' --min-alleles 2 \
--max-alleles 2 --min-af 0.005:minor \
| bcftools annotate -x ^FORMAT/GT,^INFO/R2 \
| bcftools norm --rm-dup all -Oz --output ${outvcf}

bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${outvcf} > ${snplist}
```

We may also use '-S' option of "bcftools view" to subset samples for further analysis.

In this procedure, we will have high quality genotype data of the target cohort and the list of SNPs we are interested in.


## 2.reference.panel

* requirement: bcftools

In this folder, we will process the reference panel with diverse ancestries, which is the 1000 genome high coverage data set. 
We will only keep the sites that overlapped with the cohort genotypes.

```bash
##join biallelic sites into multiallelic records
bcftools norm -m +any ${sourcevcf} \
| bcftools query -f '%CHROM:%POS:%REF:%ALT\n' > ${tmpsnplist}

##find out the shared sites with the targeted cohort
grep -f ${tmpsnplist} ${cohortsnplist} | cut -f1,2 -d":" \
| sed "s/:/\t/g"  > ${snplist}
rm ${tmpsnplist}

##filtration
### 'samples' are the selected samples of ancestry panel
bcftools view --types snps -S ${samples} -R ${snplist} ${sourcevcf} \
| bcftools annotate -x ^FORMAT/GT,INFO --force -Oz --output ${outvcf}
```

Note: you need to re-header the X chromosome if you use the data downloaded from this [link](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/), the new header file is included in this pipeline.

```bash
bcftools reheader -h 1000g.new.header.txt ${old_vcf} \
| bcftools annotate -x INFO,^FORMAT/GT -Oz -o ${new_vcf}
bcftools index -t ${new_vcf}
```

## 3.local.ancestry.inference

* requirement: flare, bcftools

From the first two steps, we have the input genotypes for local ancestry inference, the other required inputs for flare are recombination map and ancestry map for each reference sample. See flare [link] for more about the input format.

ref: Fast, accurate local ancestry inference with FLARE [link](https://www.biorxiv.org/content/10.1101/2022.08.02.502540v1#:~:text=We%20present%20FLARE%20(Fast,techniques%20developed%20for%20genotype%20imputation.)

```bash
flare=flare.jar
java -jar ${flare} ref=${refvcf} ref-panel=${refpanel} gt=${gt} \
map=${map} out=${out} nthreads=${nthreads}
bcftools index -t -f ${out}.anc.vcf.gz

```
After this step, we will have ancestry information for each allele, the output is in compressed vcf format.


*Notes:* the chromosome identifier should be the same format between the vcf of genotypes and the recombination map.


## 4.ancestry.extraction 

* requirement: extractAncestry, htslib

"cd" into the script folder and compile the c program "extractAncestry". You need to install htslib first and modify the library link in the "makefile". Please go seek help from your IT department if you have any difficulties. 

```bash
## 'ancvcf' is output of flare
## 'anc' is the ancestry code for each ancestry, see the header of flare's output
## 'outpref' is the output prefix, two outputs will be generated: ${vcfpref}.anc${anc}.vcf.gz is haplotypes from genome tracts of the specific ancestry, defined by 'anc' variable; ${vcfpref}.hapanc${anc}.vcf.gz is ancestry copy number of every variant on each haplotype.
extractAncestry=../src.v0/extractAncestry
${extractAncestry} ${ancvcf} ${anc} ${outpref}
```
In this step, we will extract dosage and genotypes from each ancestry for admixture mapping and the Tractor analysis.


## 5.Tractor analysis and admixture mapping.

* requirement: plink1.9

With all data ready we can use plink to conduct association analysis on alleles from each ancestry tracts. Do feel free to modify the codes based on your own understanding of plink for your own research interest.

```bash
## Tractor analysis
plink --vcf ${vcfpref}.anc${anc}.vcf.gz --vcf-half-call haploid \
--allow-no-sex --pheno ${phenofile} --pheno-name ${pheno} \
--covar ${covarfile} --covar-name  ${covarnames} \
--logistic beta hide-covar --out ${outpref1}.anc${anc} \
--ci 0.95 --freq case-control

paste -d' ' ${outpref1}.anc${anc}.assoc.logistic ${outpref1}.anc${anc}.frq.cc \
| sed "s/ \+/ /g" | sed "s/^ //g" | cut -f1-12,17-20 -d' ' \
| grep -v NA | gzip -c > ${outpref1}.anc${anc}.assoc.simple.gz

mv ${outpref1}.anc${anc}.assoc.simple.gz ${final}/

## admixture mapping
plink --vcf ${vcfpref}.hapanc${anc}.vcf.gz --allow-no-sex \
--pheno ${phenofile} --pheno-name ${pheno} --covar ${covarfile} \
--covar-name  ${covarnames} --logistic beta hide-covar \
--out ${outpref2}.anc${anc} --ci 0.95 --freq case-control

paste -d' ' ${outpref2}.anc${anc}.assoc.logistic ${outpref2}.anc${anc}.frq.cc \
| sed "s/ \+/ /g" | sed "s/^ //g" | cut -f1-12,17-20 -d' ' \
| grep -v NA | gzip -c > ${outpref2}.anc${anc}.admap.simple.gz

mv ${outpref2}.anc${anc}.admap.simple.gz ${final}/

```

