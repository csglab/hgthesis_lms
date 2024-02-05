#!/usr/bin/env Rscript
library(optparse)
option_list = list(
      make_option(c("-w", "--windowSize"), type="character",dest="ws",default=NULL, help="Windowsize", metavar="character"),
      make_option(c("-e", "--eventType"), type="character",dest="event",default=NULL, help="event", metavar="character"),
      make_option(c("-o", "--outDir"), type="character",dest="outdir",default=NULL, help="outputDir", metavar="character"),
      make_option(c("-a", "--affDir"), type="character",dest="affdir",default=NULL, help="affimxDir", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
options <- parse_args(opt_parser)

###########################################
# Libraries 
###########################################

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggsci)
    library(matrixStats)
    library(ggpubr)
    library(patchwork)
    library(readr)
    library(data.table)
    library(tidyr)
    library(purrr)
    library(caret)
    library(stringr)
    library(parallel)
    library(glmnet)
    library(doParallel)
    library(furrr)
})

###########################################
# Functions
###########################################

plan(multisession, workers = 10)

## Process AffiMx output
#------------------------------------
process_file <- function(file_path) {
    
  data <- fread(file_path, header = FALSE, sep = "\t",data.table=FALSE) %>%
            mutate(motif_id=sub("_spliceSites.*","",basename(file_path)))
  colnames(data)[1:2]<-c("Name","affinity")
  data <- data %>%
          mutate(event_id = sub(".*ENS","ENS",Name) %>% sub("CHR","chr",.),
                 seq_id = sub("_ENS.*","",Name) %>% sub("^([A-Z0-9]+)_","",.) %>% sub("_([0-9]*|\\d+\\.\\d+E\\+\\d+)_","_",.))
  m<-mean(data$affinity)
  s<-sd(data$affinity)
  data<-data %>%
         mutate(zscore_affinity=(affinity-m)/s)
    
 return(data)
}

## Process all files 
#------------------------------------
process_all_motif_files <- function(directory_path, motif_ids, window_size, event_type) {
    
  files <- file.path(directory_path, paste0(motif_ids,"_spliceSites_windowSize_",window_size,"_up_and_down.affinity_",event_type,".txt"))  
  mdata <- furrr::future_map_dfr(files, process_file) %>%
                      mutate(window_size=window_size,
                             event_type=event_type)  
  return(mdata)
}

# Reformat TRex identifiers 
#------------------------------------
reformat_trex_id <- function(input_string) {
   
  # Define regular expressions for extracting relevant information
  pattern_gene <- "ENSG\\d+\\.\\d+"
  pattern_chromosome <- "chr\\w+"
  pattern_strand <- "[+-]$"

  # Extract gene, chromosome, and strand from the input string
  gene <- gsub(paste0(".*(", pattern_gene, ").*"), "\\1", input_string)
  chromosome <- gsub(paste0(".*(", pattern_chromosome, ").*"), "\\1", input_string)
  strand <- gsub(paste0(".*(", pattern_strand, ").*"), "\\1", input_string)

  # Extract and reformat the coordinates
  coordinates <- sub(".*chr\\w+-", "",input_string) %>% 
                 sub("\\-\\+$|\\-\\-$","",.) %>%
                 sub("^\\d+:","",.) %>%
                 sub(":\\d+$","",.)

  # Reformat the extracted information into the desired structure
  output_string <- paste(gene, ";SE:", chromosome, ":", coordinates, ":", strand, sep = "")  

  return(output_string)
}

# Estimate number of zero-affinity motifs
#------------------------------------
get_affinity_stats<-function(motif_data){
    stats<-motif_data %>%
          filter(motif_id!="M0.6") %>%
          group_by(motif_id) %>%
          summarize(num_zero = sum(round(affinity,3)==0),
                    num_sites = length(affinity),
                    median_aff=median(affinity,na.rm=T),
                    mean_aff=mean(affinity,na.rm=T),
                    sd_aff=sd(affinity,na.rm=T))  %>%
          mutate(p_per_event=(num_zero/num_sites)*length(unique(motif_data$seq_id))) %>%
          arrange(desc(p_per_event))
    return(stats)
}

###########################################
# Define inputs
###########################################

ws<-options$ws
event<-options$event
out_res_dir<-options$outdir
aff_dir<-options$affdir

rbp_info_file<-"/home/lmoral7/lmprojects/altsplicing-methods/upstream_reg/data/cisbp_Homo_sapiens_2023_06/RBP_Information_all_motifs.txt"
motif_ids_file<-"25112023_representative_motif_ids.txt"

###########################################
# Load data
###########################################

# Motif info
message("Loading Motif information...")
motifs<-read.table(motif_ids_file,header=FALSE) %>%
        rename("motif_id"="V1")

# AffiMx output
message("Loading AffiMx outpus of ",length(motifs$motif_id)," motifs...")
motif_data<-process_all_motif_files(directory_path = aff_dir,motif_ids = motifs$motif_id,window_size = ws, event_type = event) 
print(head(motif_data,3))

message("Storing as RDS object...")
saveRDS(motif_data,file=paste0(out_res_dir,"/objects/motif_data_",event,"_",ws,".RDS"))
message("Finished successfully!")