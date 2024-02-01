suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(gridExtra)
    library(matrixStats)
    library(data.table)
    source("trex_helper.r")
})


# IPUTS

args = commandArgs(trailingOnly=TRUE)

cancer<-args[1]
event<-args[2]
trex_dir <- args[3]
out_dir <- args[4]

min.psi<-0.01 # threshold on PSI
max.inv.samples<-0.95 # max proportion of samples with invalid PSI values allowed
max.na.samples<-0.95 #max numer of NA samples


#load event counts
message("Loading event counts")
cts_dir<-file.path(trex_dir,"event_counts")
load(file.path(cts_dir,paste0(event,".eventCounts.RData")),verbose = F)

analysis<-"condition"
res_dir<-file.path(trex_dir,"coefs_condition")
res.files <- list.files(file.path(res_dir),pattern = paste0(event,".res.lfcShrink.condition.trex_countA*"),recursive = F)
names(res.files)<-sub("\\.res.*","",res.files)%>%sub("/","_",.)

cond <- lapply(res.files,function(file,...){
                 res <- suppressWarnings(fread(file.path(res_dir,file),data.table=FALSE,stringsAsFactors = F,drop=1,nThread = 4))%>%
                        mutate(signif = padj<0.05) %>%                
                        mutate(cancer=cancer,
                                exp_var=sub(".*trex_countA\\.","",file) %>% sub("\\.tsv","",.)) %>%
                         select(-constitutive_transcripts,-alternative_transcripts,-seqname)    
        }) %>% do.call(rbind,.)%>%
        rename_at(vars(log2FoldChange:padj,signif),function(col){paste0(col,"_condition")} )


# split counts
message("splitting counts")
cts.tumor<-trex.meta@assays@data$counts[unique(cond$event_id),trex.meta@colData$condition=="tumor"]
cts.normal<-trex.meta@assays@data$counts[unique(cond$event_id),trex.meta@colData$condition=="normal"]


#estimate PSI
message("Estimating PSI")
psi.tumor<-cts.tumor[,grepl("A_",colnames(cts.tumor))]/(cts.tumor[,grepl("A_",colnames(cts.tumor))]+cts.tumor[,grepl("C_",colnames(cts.tumor))])
psi.normal<-cts.normal[,grepl("A_",colnames(cts.normal))]/(cts.normal[,grepl("A_",colnames(cts.normal))]+cts.normal[,grepl("C_",colnames(cts.normal))])

#Calculate PSI stats
message("Calculating stats")
psi.stats<-data.frame(mean_psi_tumor = rowMeans(psi.tumor,na.rm = T),
                      mean_psi_normal = rowMeans(psi.normal,na.rm = T), 
                      na_tumor = rowSums(is.na(psi.tumor))/ncol(psi.tumor),
                      na_normal = rowSums(is.na(psi.normal))/ncol(psi.normal),
                      nun_low_tumor = (rowSums(psi.tumor <= min.psi,na.rm=T))/ncol(psi.tumor),
                      nun_low_normal = (rowSums(psi.normal <= min.psi,na.rm=T))/ncol(psi.normal),
                      nun_high_tumor = (rowSums((1-psi.tumor) <= min.psi,na.rm=T))/ncol(psi.tumor),
                      nun_high_normal = (rowSums((1-psi.normal) <= min.psi,na.rm=T))/ncol(psi.normal)) %>%
                mutate(pi_tumor = nun_low_tumor+nun_high_tumor+na_tumor,
                       pi_normal = nun_low_normal+nun_high_normal+na_normal,
                       p_na = na_tumor+na_normal) %>%
                select(-nun_low_tumor,-nun_low_normal,-nun_high_tumor,-nun_high_normal)

psi.stats<-psi.stats %>% 
            mutate(psi_flag = (mean_psi_tumor <= min.psi & mean_psi_normal <= min.psi) | ((1-mean_psi_tumor) <= min.psi & 1-mean_psi_normal),
                   samples_flag = (pi_tumor >= max.inv.samples & pi_normal >= max.inv.samples),
                   na_flag = (na_tumor >= max.na.samples | na_normal >= max.na.samples)) %>%
            mutate(event_flag = ifelse((psi_flag&samples_flag)|na_flag,"invalid","valid")) %>%
            tibble::rownames_to_column("event_id")


message("Writting table")
write.table(psi.stats,file=paste0(out_dir,"/",event,"_PSI_",min.psi,"_flags.tsv"),row.names=F,col.names=T,sep="\t",quote=FALSE)


message("Finished successfully!")
