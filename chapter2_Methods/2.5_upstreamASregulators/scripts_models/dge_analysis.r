#!/usr/bin/env Rscript
library(optparse)
option_list <- list(
                    make_option(c("-c", "--cancer"), type="character",dest="ctype",default=NULL,
                                help="cancer", metavar="character"),
                    make_option(c("-i", "--inObject"), type="character",dest="inobj",default=NULL, 
                                help="inputObject", metavar="character"),
                    make_option(c("-o", "--outDir"), type="character",dest="outdir",default=NULL,
                                help="outputDir", metavar="character")
                    )

# Packages
#----------------------------------------------------------------

suppressPackageStartupMessages({
    library(getopt)
    library(dplyr)
    library(data.table)
    library(DESeq2)
    library(readr)
    library(tidyr)
})
options<-parse_args(OptionParser(option_list = option_list))


# Functions
#----------------------------------------------------------------

# Extract metadata and gene expression counts for one cancer type
process_cancer_data<-function(ctype,...){
    # Exract cancertype data 
    meta_cancer<-metadata %>% filter(cancer==ctype)
    ns<-nrow(meta_cancer)
    # Process covariates
    meta_cancer <- meta_cancer %>%
                   mutate(condition = factor(as.character(condition),levels=c("normal","tumor")),
                          sex = factor(ifelse(gender=="male",0,1),levels = c(0,1)),
                          age = age_at_index) %>%
                   filter(!(is.na(age)|is.na(sex)),
                          (!is.na(impurity) & condition=="tumor") | condition=="normal") %>%
                   select(submitter_id,condition,sex,age,impurity)
    rownames(meta_cancer)<-meta_cancer$submitter_id

    # Remove samples without impurity estimations
    nt<-sum(meta_cancer$condition=="tumor")
    nn<-sum(meta_cancer$condition=="normal")
    message("Removing ",ns-nt," tumor samples without purity estimations")

    # Filter events with low counts 
    cts<-round(as.matrix(gn_cts[,rownames(meta_cancer)]),0) 
    ngenes<-nrow(cts)
    cts<-cts[rowSums(cts >= min_counts) >= nn,]
    message("Filtering ",ngenes-nrow(cts)," low expressed genes")
    
    res_lis<-list(counts=cts,
                  meta=meta_cancer)
    return(res_lis)
}


# Run DESeq2
run_deseq<-function(deseq_inputs){
    
    message("Running DESeq2...")
    dds <- DESeqDataSetFromMatrix(countData = deseq_inputs$counts,
                                  colData = deseq_inputs$meta,
                                  design = ~ condition + sex + age + impurity)
    dds <- DESeq(dds, quiet = FALSE)
    
    return(dds)
} 

# Main
#----------------------------------------------------------
ctype<-options$ctype
inputs.file<-options$inobj
res_dir<-options$outdir
min_counts<-5

# Load data
message("Loading input data...")
load(inputs.file,verbose=F)

# Run analysis
message("Starting analysis...")
cancer_list<-process_cancer_data(ctype=ctype)
dds<-run_deseq(deseq_inputs=cancer_list)

# Save results
message("Saving dds object")
outfile<-paste0(res_dir,"/dds.RData")
save(dds,file=outfile)

contrasts<-resultsNames(dds)[c(-1)]
for(contrast in contrasts){
    
    res <- results(dds,name=contrast) %>% as.data.frame()
    write.table(res, quote=FALSE, sep='\t',
                file=paste(res_dir,"/res.",contrast,".tsv",sep=""))
    
    res <- lfcShrink(dds,coef=contrast,format = "DataFrame",quiet=TRUE,type = "normal") %>%
            as.data.frame() 
    write.table(res,quote=F,sep='\t',row.names = FALSE,
               file=paste(res_dir,"/res.lfcShrink.",contrast,".tsv",sep=""))
}
message("Finished successfully!")