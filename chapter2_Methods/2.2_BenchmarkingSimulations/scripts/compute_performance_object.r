suppressPackageStartupMessages({
  library(dplyr)
  library(binr)
  library(RColorBrewer)
  library(cowplot) 
  library(data.table)
  library(ggpubr)
  library(SummarizedExperiment)
  library(ggplot2)
  library(viridis)
  library(grid)
  library(gridExtra)
  library(patchwork)
  library(ggpubr)
  library(ggsci)
  library(caret)
  library(precrec)
  library(caret)
  library(parallel)
  library(plotROC)
})

### FUNCTIONS

get_class_labels<-function(dpsi,dlog,...){
    ref <- ref.bm %>%
            mutate(class_labels = ifelse( abs(ref_dPSI)>=dpsi & abs(ref_logit_dPSI)>=dlog,1,0)) %>%
            mutate(num_pos = sum(class_labels==1),
                   num_neg = sum(class_labels==0),
                   dPSI_th = dpsi,
                   dLogitPSI_th = dlog)
    return(ref)
}

compute_auroc<-function(ref.data,preds.data){
    
    # Extract data frames

    preds.data <- inner_join(preds.data,
                             ref.data %>% select(event_id,class_labels,num_pos,num_neg,dPSI_th,dLogitPSI_th),
                             by="event_id") 
    
    pd.notna<-preds.data %>% na.omit()
    if(nrow(pd.notna)>=30){
        np<-sum(pd.notna$class_labels==1)
        nn<-sum(pd.notna$class_labels==0)
    }else{
        np<-0
        nn<-0
    }

    
    ## Compute AUROC
    
    if(np>=15 & nn>=15){
        auc<-evalmod(scores = preds.data$pred_score, 
                     labels = preds.data$class_labels,
                     modnames = 'auc',
                     na_worst = TRUE,
                     mode = "aucroc")  
    }else{
        auc<-list(uaucs=data.frame('modnames'= 'auc',
                                 'dsids'= NA,
                                 'aucs'= NA,
                                 'ustats'= NA))
    }
    
    ## Compute correlations
    
    if(sum(is.na(preds.data$dPSI))!=nrow(preds.data)){
        preds.data<-preds.data %>% 
                    ungroup() %>%
                    filter(!is.na(dPSI) & !is.na(ref_dPSI) & !is.na(ref_logit_dPSI)) %>%
                    filter(abs(ref_dPSI) > dPSI_th, 
                           abs(ref_logit_dPSI) > dLogitPSI_th) 
        
        cr1<-cor(x=preds.data$dPSI,y=preds.data$ref_dPSI)
        cr2<-cor(x=preds.data$dPSI,y=preds.data$ref_logit_dPSI)
    }else{
        c1<-NA
        c2<-NA
        
    }

    
    perf.df <- auc$uaucs %>%
               as_tibble() %>%
               mutate(num_pos = np,
                      num_neg = nn,
                      num_events = num_pos + num_neg,
                      cor_dPSI = cr1,
                      cor_dLogitPSI = cr2)
    return(perf.df)
}



auroc_wrapper<-function(row.data,...){
    
    row.perf <- ref.bm.labels %>%
                mutate(performace=purrr::map(.x=data,.f=compute_auroc,preds.data=row.data)) %>%      
                select(-data) %>%
                tidyr::unnest(c(performace)) 
    return(row.perf)

}

### MAIN

message("Loading performance objects")
load("tmp/perf.data.niso.RData",verbose=T)

message("Computing performace")
perf.data.niso <- perf.data.niso %>%
                  mutate(perf_list = mclapply(data,auroc_wrapper,ref.bm.labels,mc.cores = 10)) %>% 
                  select(-data)%>% 
                  tidyr::unnest(c(perf_list))


message("Saving performance results")
saveRDS(perf.data.niso,file = "../results/performance_metrics_numisos_v2.RDS")

message("Finished successfully!")