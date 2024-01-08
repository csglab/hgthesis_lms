suppressPackageStartupMessages({
    library(getopt)
    library(optparse)
    library(dplyr)
    library(data.table)
    library(DESeq2)
})


# Input args

option_list <- list(
  make_option(c("-s", "--sample"), type='character',dest="sampleid",
              help="Sample identifier used to name the output files"),
  make_option(c("-e", "--event"), type='character',dest="evtype",
              help="Name of AS event type to quantify. Must be one of A3,A5,AF,AL,MX,RI,SE,ALL"),
  make_option(c("-d", "--directory"), type='character',dest="dir",
              help="Directory with evenCounts"),
  make_option(c("-r","--refmap"), type="character",dest = "isomapfile",
              help="Path to the ioe file from SUPPA2")
)

options<-parse_args(OptionParser(option_list = option_list))
ransample <- options$sampleid
data_dir <- options$dir
evtype<-options$evtype 

min_ecount<-10 # Min event counts

message("------------ Analyzing ",ransample,"------------")
message("---- Fitting models for ",evtype," events -----")

isomapfile<-options$isomapfile
outdir <- data_dir
ev_cts <- paste0(outdir,"/",evtype,".eventCounts.RData")

# Load data
map <- fread(isomapfile,data.table = F,stringsAsFactors = F) %>% 
        as.data.frame() %>%
        dplyr::filter(event_type==evtype)
load(ev_cts,verbose = F)

trex.meta<-trex.meta[colSums(assays(trex.meta)[["counts"]])!=0,] 
ns<-ncol(trex.meta)

# Analysis without covariates
#######################################

message("Fitting model without covariates")
analysis<-"cell_line"
m<-model.matrix(~trex_count + cell_line + cell_line:trex_count, trex.meta@colData)

# Create DESeq2 dataset
cts <- trex.meta[,rownames(m)]
assays(cts) <- assays(cts)[c("counts","length")]
dds <- suppressMessages(DESeqDataSet(cts, design=m))

# Filter low count events
keep <- rowSums(counts(dds)) >= min_ecount
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds, quiet = FALSE, minReplicatesForReplace=Inf)

# Save results
message("Saving dds object")
outfile<-paste0(outdir,"/",evtype,".dds.",analysis,".RData")
save(dds,file=outfile)

message("Saving res files")
contrasts<-resultsNames(dds)[-1]
for(contrast in contrasts){
    
    res <- results(dds,name=contrast) %>%
           as.data.frame() %>%
           tibble::rownames_to_column("event_id") %>%
           left_join(.,map %>% select(-total_transcripts),by = "event_id")
    outfile<-paste(outdir,"/",evtype,".res.",analysis,".",contrast,".tsv",sep="")
    write.table(as.data.frame(res),file=outfile,quote=F,sep='\t')
}



# Analysis with covariates
#######################################

message("Fitting model with covariates...")
analysis<-"cell_line.batch"
m<-model.matrix(~ trex_count + cell_line + batch + (cell_line+batch):trex_count, trex.meta@colData)

# Create DESeq2 dataset
cts <- trex.meta[,rownames(m)]
assays(cts) <- assays(cts)[c("counts","length")]
dds <- suppressMessages(DESeqDataSet(cts, design=m))

# Filter low count events
keep <- rowSums(counts(dds)) >= min_ecount
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds, quiet = FALSE, minReplicatesForReplace=Inf)

# Save results
message("Saving dds object")
outfile<-paste0(outdir,"/",evtype,".dds.",analysis,".RData")
save(dds,file=outfile)

message("Saving res files")
contrasts<-resultsNames(dds)[c(-1,-2)]
for(contrast in contrasts){
    
    res <- results(dds,name=contrast) %>%
           as.data.frame() %>%
           tibble::rownames_to_column("event_id") %>%
           left_join(.,map %>% select(-total_transcripts),by = "event_id")
    outfile<-paste(outdir,"/",evtype,".res.",analysis,".",contrast,".tsv",sep="")
    write.table(as.data.frame(res),file=outfile,quote=F,sep='\t')
}
message("Completed analysis!")

