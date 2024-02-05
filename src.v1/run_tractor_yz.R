#!/usr/bin/env Rscript

# version: 1.1.0
script_version <- "1.1.0 +"
cat(paste("Tractor Script Version:", script_version), "\n")

# Changelog:
# modified by Ying Zhou to 
### 1) support gzip input 
### 2) allow mismatches of samples between phenotypes and genetypes, only ovelapped samples will be used
### 3) mem-friendly to large files, stream the input line by line
# Python script to R script to better deal with overly-inflated values.
# This script will instead output NA's in such cases.

require("optparse")
library(stringr) 

option_list = list(
  make_option(c("--hapdose"), type="character", default=NULL,
              help="The prefix of hapcount and dosage files", metavar="character"
),
  make_option(c("--phe"), type="character", default=NULL,
              help="Phenotype and covariates file;
                1st column sample ID, 2nd column phenotype, other columns will be
 treated as covariates;
                The table should only contain numeric/integer values;
                Missing data is allowed", metavar="character"),
  make_option(c("--method"), type="character", default=NULL,
              help="Either <linear> or <logistic>", metavar="character"),
  make_option(c("--out"), type="character", default=NULL,
              help="Summary statistics file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



if (is.null(opt$hapdose)){
  print_help(opt_parser)
  stop("Missing hapdose files", call.=FALSE)
} else if (is.null(opt$phe)){
  print_help(opt_parser)
  stop("Missing phenotype/covariates file", call.=FALSE)
} else if (!(opt$method %in% c("linear", "logistic"))){
  print_help(opt_parser)
  stop("Method should be either <linear> or <logistic>", call.=FALSE)
} else if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Missing output file", call.=FALSE)
}

####### helper function
subset_mat_NA = function(rows, mat){
  t(sapply(rows, function(row, mat){if (row %in% row.names(mat)){return(mat[row,]
)}else{return(rep(NA, ncol(mat)))}},mat))
}



###########
RunTractor <- function(prefix, phefile, method, outfile){
#prefix="test"
#method = "logistic"
#phefile="test.pheno_covars.txt"
#outfile = paste0(prefix, ".out.test0.txt")

  hapFiles = Sys.glob(paste0(prefix, "*.hapcount.txt*"))
    doseFiles = Sys.glob(paste0(prefix, "*.dosage.txt*"))
    nAnc = length(doseFiles)
    inFiles =c(hapFiles, doseFiles)
    phe = read.csv(phefile, sep = "\t", header = T)
    ID = as.character(phe[,1])  ### --> yz

    inFileCons = lapply(inFiles,function(infile){
	file(infile, open = "rt")
	})

    #data = lapply(inFiles,function(file){read.table(file, header=T)})

# overlaped samples between genotypes and phenotypes
    hdrs = lapply(inFileCons, function(con){unlist(strsplit(readLines(con, n=1), '\t'))})
    ID = ID[ID %in% hdrs[[1]][-(1:5)]]

    if(length(ID) == 0) stop("[Exit] No overlapped samples are selected")

# updates all files
    rownames(phe) = phe[,1]
    phe = phe[ID, ]
    sel_colmn = c("CHROM", "POS", "ID", "REF","ALT", ID)

    ID = phe[,1]  ### --> yz
    y = phe[,2] ### --> yz

    if (ncol(phe) > 2){
      COV_ = as.matrix(phe[,3:ncol(phe)])
    } else {
      COV_ = NULL
    }

  LAcolnames = paste0("LA", 0:(nAnc-1))
    Gcolnames = paste0("G", 0:(nAnc-1))

    resDF = setNames(data.frame(matrix(data = NA, nrow = 0, ncol = (5 +  4 * nAnc + 2 * (nAnc - 1)))),
	c("CHR", "POS", "ID", "REF", "ALT",
	  paste0("AF_anc",0:(nAnc-1)),
	  paste0("LAprop_anc",0:(nAnc-1)),
	  paste0("LAeff_anc",0:(nAnc-2)),
	  paste0("LApval_anc",0:(nAnc-2)),
	  paste0("Geff_anc",0:(nAnc-1)),
	  paste0("Gpval_anc",0:(nAnc-1))))
    write.table(resDF, outfile,  quote = F, row.names = F, sep = "\t")

    i = 0    
    while(T){
      i = i + 1
# matrix of Local Ancestry and Genotype
	data = unlist(lapply(inFileCons, 
	      function(con){
	      row = unlist(strsplit(readLines(con, n=1), '\t'))
	      if (length(row) == 0) {quit(save="no")}
	      names(row) = hdrs[[1]]
	      row[sel_colmn] } ))
	if (is.null(data))  break 
	LAG_ = matrix(data, nrow = length(sel_colmn), ncol = length(inFileCons))
	snp_info = LAG_[1:5, 1]
	LAG_ = LAG_[-(1:5), ]
	class(LAG_) <- "numeric"
	colnames(LAG_) = c(LAcolnames, Gcolnames)

	AF = as.numeric(colSums(LAG_[,Gcolnames])/colSums(LAG_[,LAcolnames]))
	LAprop = as.numeric(colSums(LAG_[,LAcolnames])/sum(LAG_[,LAcolnames]))

	LAG_ = LAG_[,c(LAcolnames[1:(length(LAcolnames) - 1)], Gcolnames)]
	coef_rownames = paste0("LAG_", colnames(LAG_))
	LA_rownames = paste0("LAG_", LAcolnames[1:(length(LAcolnames) - 1)])
	G_rownames = paste0("LAG_", Gcolnames)

	if (method == "linear"){
	  if (!is.null(COV_)){
	    coefs = summary(lm(y ~ LAG_ + COV_))$coefficients
	  } else {
	    coefs = summary(lm(y ~ LAG_))$coefficients
	  }

	  reg_res = subset_mat_NA(coef_rownames, coefs[,c(1,4)])
	    LAeff = as.numeric(reg_res[LA_rownames,1])
	    LApval = as.numeric(reg_res[LA_rownames,2])
	    Geff = as.numeric(reg_res[G_rownames,1])
	    Gpval = as.numeric(reg_res[G_rownames,2])


	} else if (method == "logistic"){
	  if (!is.null(COV_)){
	    model = glm(y ~ LAG_ + COV_, family = binomial(link = "logit"))
	      coefs = summary(model)$coefficients
	  }else {
	    model = glm(y ~ LAG_, family = binomial(link = "logit"))
	      coefs = summary(model)$coefficients
	  }

	  reg_res = subset_mat_NA(coef_rownames, coefs[,c(1,4)])
	    if (model$converged == TRUE){
	      LAeff = as.numeric(reg_res[LA_rownames,1])
		LApval = as.numeric(reg_res[LA_rownames,2])
		Geff = as.numeric(reg_res[G_rownames,1])
		Gpval = as.numeric(reg_res[G_rownames,2])
	    } else {
# if glm doesn't converge, set effect size and P value as NA
	      LAeff = rep(NA, length(LA_rownames))
		LApval = rep(NA, length(LA_rownames))
		Geff = rep(NA, length(G_rownames))
		Gpval = rep(NA, length(G_rownames))
	    }
	}


      resDF[1, c("CHR", "POS", "ID", "REF", "ALT")] = snp_info
	resDF[1, paste0("AF_anc",0:(nAnc-1))] = round(AF,6)
	resDF[1, paste0("LAprop_anc",0:(nAnc-1))] = round(LAprop,6)
	resDF[1, paste0("LAeff_anc",0:(nAnc-2))] = round(LAeff, 6)
	resDF[1, paste0("LApval_anc",0:(nAnc-2))] = LApval
	resDF[1, paste0("Geff_anc",0:(nAnc-1))] = round(Geff,6)
	resDF[1, paste0("Gpval_anc",0:(nAnc-1))] = Gpval

	write.table(resDF,  outfile, quote = F, row.names = F, col.names = F, 
	    append= T, sep = "\t")
    }
}

RunTractor(prefix = opt$hapdose, phefile = opt$phe, 
    method = opt$method, outfile = opt$out)
