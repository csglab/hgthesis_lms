#!/usr/bin/env Rscript
library(optparse)
option_list = list(
      make_option(c("-w", "--windowSize"), type="character",dest="ws",default=NULL, help="Windowsize", metavar="character"),
      make_option(c("-c","--cancer"), type="character",dest="cancertype", default=NULL, help="cancerType",metavar="character"),
      make_option(c("-e", "--eventType"), type="character",dest="event",default=NULL, help="event", metavar="character"),
      make_option(c("-t", "--trexObject"), type="character",dest="trexobj",default=NULL, help="Trex res object", metavar="character"),
      make_option(c("-o", "--outDir"), type="character",dest="outdir",default=NULL, help="outputDir", metavar="character")
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
          summarize(num_zero = sum(round(affinity,6)==0),
                    num_sites = length(affinity),
                    median_aff=median(affinity,na.rm=T),
                    mean_aff=mean(affinity,na.rm=T),
                    sd_aff=sd(affinity,na.rm=T))  %>%
          mutate(p_zeros=num_zero/num_sites,
                 p_per_event=(num_zero/num_sites)*length(unique(motif_data$seq_id))) %>%
          arrange(desc(p_per_event))
    return(stats)
}

# Split dataset
#------------------------------------
split_data<-function(motif_data,p){
    
    train_ids<-sample(x = 1:nrow(motif_data),replace = FALSE,size = p*nrow(motif_data))
    train_data<-motif_data[train_ids,]
    test_data<-motif_data[-1*train_ids,]
    
    res<-list(train_data=train_data,
              test_data=test_data)
    return(res)
}

# Extract model coefficients from list
#------------------------------------

extract_model_coefficients<-function(model,fit_method){
    if(fit_method=="lm"){
        coefs <- summary(model)$coefficients %>%
                 as.data.frame()%>%
                 tibble::rownames_to_column("coefficient")
    }else if(fit_method=="glmnet"){
    
        coefs<-coef(model$finalModel,s=model$bestTune$lambda) %>%
                as.matrix() %>%
                as.data.frame() %>%
                tibble::rownames_to_column("coefficient")%>%
                rename("Estimate"="s1")
    }
    return(coefs)
}

# Run model training steps 
#------------------------------------

training_feature_selection<-function(train_data,fit_method="glmnet"){
    
    message("Running model training and feature selection...")
    train_ctrl <- trainControl(method = "cv", number = 10, p=0.7)
    train_model <- train(log2FoldChange ~ ., 
                   data = train_data, 
                   method = fit_method, 
                   trControl = train_ctrl)
    coefficients <- extract_model_coefficients(model = train_model,fit_method = fit_method) 
    res<-list(coefficients=coefficients,
              train_model=train_model)
    return(res)
}

# Run model testing steps
#------------------------------------
test_model<-function(test_data,fit_method="lm"){
    
    message("Testing model...")
    ctrl <- trainControl(method = "cv", number = 10, p=0.7)
    test_model <- train(log2FoldChange ~ ., 
                   data = test_data, 
                   method = fit_method, 
                   trControl = ctrl)
    coefs <- extract_model_coefficients(model = test_model,fit_method)
    
    res<-list(coefficients=coefs,
              test_model=test_model)
}

# Extract performance tables 
#------------------------------------

extract_performance_table<-function(training,testing){
    
    message("Extracting performance results...")
    
    train_perf<-training$train_model$results %>%
            filter(lambda==training$train_model$bestTune$lambda,
                   alpha==training$train_model$bestTune$alpha) %>%
            mutate(stage="training",
                   R = sqrt(Rsquared),
                   n = floor(median(unlist(lapply(training$train_model$control$index,length))))) %>%
            select(-alpha,-lambda)

    test_perf<-testing$test_model$results %>%
                mutate(stage="testing",
                       R = sqrt(Rsquared),
                       n = floor(median(unlist(lapply(testing$test_model$control$index,length))))) %>%
                select(-intercept)

    performance<-rbind(train_perf,test_perf)
    
    return(performance)
}

# Run model fit procedure
#------------------------------------

lm_model_fit_all<-function(motif_data_long){
        
    outfile<-unique(motif_data_long$outfile)
    
    motif_data <- motif_data_long %>%
                  select(-outfile) %>%
                  group_by(event_id,log2FoldChange) %>%
                  tidyr::pivot_wider(names_from = motif_seq_id,values_from = affinity,values_fn = max) %>% 
                  tibble::column_to_rownames("event_id")
    
    # Split data
    data_split<-split_data(motif_data,p=0.25)

    # Train model and select features
    training_fs<-training_feature_selection(train_data = data_split$train_data,fit_method = "glmnet")
    features<-training_fs$coefficients %>%
                filter(Estimate!=0,coefficient!="(Intercept)") %>%
                select(coefficient) %>%
                unlist() %>%
                as.character() 
    
    # Store coefficients table
    write.table(training_fs$coefficients,
                file=sub("coefficients","training_coefficients",outfile),
                row.names=FALSE,col.names=TRUE,quote=F,sep="\t")
    
    
    # Test model on selected features only
    #final_data<-motif_data[,c("log2FoldChange",features)]
    final_data<-motif_data
    final_model<-test_model(final_data,fit_method = "lm")
    
    # Store coefficients table
    write.table(final_model$coefficients,file=outfile,row.names=FALSE,col.names=TRUE,quote=F,sep="\t")
    
    # Extract training and testing performance
    model_res<-extract_performance_table(training = training_fs,testing = final_model)
    
    stopCluster(cl)
    
    message("Completed model fit!")
    return(model_res)
}

###########################################
# Define inputs
###########################################

ws<-options$ws
ctype<-options$cancertype
coef.sh.obj.fg.file<-options$trexobj
event<-options$event
out_res_dir<-options$outdir

aff_dir<-"../results/upstream_rbp_AffiMx_forward/by_event"
rbp_info_file<-"/home/lmoral7/lmprojects/altsplicing-methods/upstream_reg/data/cisbp_Homo_sapiens_2023_06/RBP_Information_all_motifs.txt"
motif_ids_file<-"../scripts_rbp_models/all_motifs.txt"

cl <- makePSOCKcluster(2)
registerDoParallel(cl)

###########################################
# Load data
###########################################

message("Analyzing ",event," events in ",ctype,"...")
outfile<-paste0(out_res_dir,"/performance_reports/models_summary_windowSize_",ctype,"_",event,"_",ws,".tsv")
#if(file.exists(outfile)){
#    message("Found results for this cancer type, exiting.")
#    exit()
#}

# Trex predictions for conditiontumor, valid events only
message("Loading TRex predictions...")
cond.res.fg<-readRDS(file=coef.sh.obj.fg.file) %>%
             filter(event_type==event) 

# Load pre-computed data object
motif_data<-readRDS(file=paste0(out_res_dir,"/objects/motif_data_",event,"_",ws,".RDS"))
message("Mofis scanned in ",nrow(motif_data)," sequences")

# Filter low affinity motifs
mstats<-get_affinity_stats(motif_data) %>%
        mutate(remove_motif=p_zeros>=0.75) 
print(head(mstats))
valid_motifs<-mstats$motif_id[!mstats$remove_motif]
motif_data<-motif_data %>%
            filter(motif_id %in% valid_motifs)
message("Discarded ",sum(mstats$remove_motif)," motifs due to high number of zero affinity sites")

###########################################
# Fit models
###########################################

message("Preparing data...")
motif_data_lis <- motif_data %>%
                  inner_join(.,cond.res.fg) %>%
                  mutate(motif_seq_id = paste(motif_id,seq_id,sep="_")) %>%
                  select(event_type,cancer,motif_seq_id,event_id,log2FoldChange,affinity) %>%
                  group_by(cancer) %>%
                  mutate(outfile=paste0(out_res_dir,"/coefficients/",cancer,"_",event_type,"_",ws,".tsv")) %>%
                  tidyr::nest(model_data=c(outfile,event_id,log2FoldChange,motif_seq_id,affinity)) 

message("Running model fits...")
options(future.rng.onMisuse="ignore")
motif_data_lis_fit <- motif_data_lis %>%
                      ungroup() %>%
                      mutate(res = furrr::future_map(model_data,lm_model_fit_all,.progress = FALSE)) 

# Write table of performance metrics

message("Writing performance table...")
perf<-motif_data_lis_fit %>%
      select(-model_data,-res) 
res.df<-cbind(perf,do.call(rbind,motif_data_lis_fit$res)) %>%
         mutate(windowSize=ws)
write.table(res.df,
            file = outfile,
            row.names = F,col.names = T,quote=F,sep="\t")

message("Finished fitting ",nrow(res.df)," models successfuly!")

