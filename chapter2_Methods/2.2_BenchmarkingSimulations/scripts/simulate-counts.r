library(readr)
library(data.table)
library(plyr)
library(VGAM)
library(drc)

# Input settings 

set.seed(7)

args = commandArgs(trailingOnly=TRUE)
dir_input<-args[1] # directory containing the files ground_truth_TPM.tsv and metadata.tsv
dir_results<-args[2] # directory to store the results
sample_file<-args[3] # file with sample ids 
dir_efsizes<-file.path(dir_results,"gt_effects") # sub directory to store ground truth

dir.create(dir_results,showWarnings = FALSE)
dir.create(dir_efsizes,showWarnings = FALSE)

# Load data
#####################################

message("Loading data from ",sample_file,"...")
iso_res_rsem<-read.table(sample_file,header=TRUE)
obs_TPM <- fread(file.path(dir_input,"ground_truth_TPM.tsv"),data.table=FALSE)
n_transcripts <- nrow(obs_TPM)
simulation_design <- fread(file.path(dir_input,"metadata.tsv"))
pseudocount <- 1e-3


# Fit mean-logFC model 
message("Fitting models on the ground truth...")
nonZero_TPMs <- rowSums( obs_TPM > 0 ) 
nonZero_indices <- which( nonZero_TPMs >= 4 ) 

log_obs_TPM_c1 <- ( log(obs_TPM[nonZero_indices,1]+pseudocount) + log(obs_TPM[nonZero_indices,2]+pseudocount) ) / 2
log_obs_TPM_c2 <- ( log(obs_TPM[nonZero_indices,3]+pseudocount) + log(obs_TPM[nonZero_indices,4]+pseudocount) ) / 2
obs_logFC <- log_obs_TPM_c2 - log_obs_TPM_c1

cv.data <- data.frame(mean_exp = (log_obs_TPM_c1+log_obs_TPM_c2)/2,logFC = obs_logFC )
cv.data <- cv.data[ order(cv.data$mean_exp), ]
cv.data$sd <- sapply( 1:nrow(cv.data), function(x) sd(cv.data$logFC[max(x-500,1):min(x+500,nrow(cv.data))]) )
                     
logFC_vs_mean_fit <- drm(sd ~ mean_exp, data = cv.data, fct = L.4(), type = "continuous")
                     
# Fit mean-variance model 
residual1 <- log(obs_TPM[nonZero_indices,1]+pseudocount) - log_obs_TPM_c1
residual2 <- log(obs_TPM[nonZero_indices,2]+pseudocount) - log_obs_TPM_c1
residual3 <- log(obs_TPM[nonZero_indices,3]+pseudocount) - log_obs_TPM_c2
residual4 <- log(obs_TPM[nonZero_indices,4]+pseudocount) - log_obs_TPM_c2
all_residuals <- c(residual1,residual2,residual3,residual4)

cv.data <- data.frame(mean_exp = c(log_obs_TPM_c1,log_obs_TPM_c1,log_obs_TPM_c2,log_obs_TPM_c2), residual = all_residuals )
cv.data <- cv.data[ order(cv.data$mean_exp), ]
cv.data$sd <- sapply( 1:nrow(cv.data), function(x) sd(cv.data$residual[max(x-500,1):min(x+500,nrow(cv.data))]) )
cv.data <- cv.data[cv.data$sd > 0 ,] # filter entries with sd = 0

# Fit a line that describes the logarithm of sd as a function of mean expression
res_vs_mean_fit <- drm(sd ~ mean_exp, data = cv.data, fct = L.4(), type = "continuous")   

# Generate TPM for c1 from a normal distribution
obs_TPM_NA <- obs_TPM
obs_TPM_NA[ obs_TPM_NA == 0 ] <- NA
obs_TPM_NA_mean <- apply(obs_TPM_NA,1,mean,na.rm=T)

c1_probs <- rank(obs_TPM_NA_mean,na.last="keep",ties.method="random")
c1_probs <- c1_probs / max(c1_probs+1,na.rm = T)

c1_probs[ is.na(c1_probs) ] <- min(c1_probs,na.rm=T)/1000

log_gt_TPM_c1 <- qnorm(c1_probs, 
                       mean=mean(log_obs_TPM_c1),
                       sd=sd(log_obs_TPM_c1))

design_mx <- model.matrix( ~ cell_line + batch, simulation_design )
                     
# Estimate condition effect
message("Estimating ground truth effects...")                     
gt_condition_effect <- rlaplace(n_transcripts, 
                                location=0,
                                scale=predict(logFC_vs_mean_fit,data.frame(mean_exp=log_gt_TPM_c1) )/sqrt(2))
gt_condition_effect[ nonZero_TPMs==0 ] <- 0 

# Estimate batch effect 
message("Running simulations...")                     
nran<-5
bes_list<-c(0.25,0.5,0.75,1)
for(n in 1:nran){
    
    message("------")
    message("Run = ",n)
    message("------")
    for(bes in bes_list){
        
        message("strenght = ",bes)
         
        gt_batch_effect <- rlaplace(n_transcripts, 
                                    location=0,
                                    scale=predict(logFC_vs_mean_fit,data.frame(mean_exp=log_gt_TPM_c1))/sqrt(2)*bes)
        gt_batch_effect[ nonZero_TPMs==0 ] <- 0 

        # Generate logTPM 
        gt_effects <- as.matrix(data.frame(intercept = log_gt_TPM_c1,
                                           condition_effect = gt_condition_effect,
                                           batch_effect = gt_batch_effect ))

        gt_log_TPM_mx <- gt_effects %*% t(design_mx)
        colnames(gt_log_TPM_mx)<-simulation_design$Sample

        noise <- rlaplace(length(gt_log_TPM_mx[nonZero_TPMs>0,]), 
                           location=0,
                           scale=predict(res_vs_mean_fit,data.frame(mean_exp=c(gt_log_TPM_mx[nonZero_TPMs>0,])))/sqrt(2))
        gt_log_TPM_mx[nonZero_TPMs>0,] <- gt_log_TPM_mx[nonZero_TPMs>0,] +
                                          noise

        # Convert logTMP values to TPM
        gt_TPM_mx <- exp(gt_log_TPM_mx)
        colnames(gt_TPM_mx) <- simulation_design$Sample
        gt_TPM_mx[nonZero_TPMs==0,] <- 0
        gt_TPM_mx <- scale(gt_TPM_mx, center = F, scale=colSums(gt_TPM_mx)/1000000)
        
        for(sample in colnames(gt_TPM_mx)){
            
            iso_res_rsem_batch<-iso_res_rsem
            iso_res_rsem_batch$TPM<-gt_TPM_mx[,sample]
            rownames(iso_res_rsem_batch)<-NULL

            dout<-file.path(dir_results,sample)
            dir.create(dout,showWarnings = FALSE)
            fout<-paste0(dout,"/",bes,"_rs.",n,"_rsem.isoforms.results")
            write.table(iso_res_rsem_batch,file = fout,sep="\t",quote=FALSE,row.names = F)
        }
        
        # Store ground truth effect sizes 
        fout<-paste0(dir_efsizes,"/",bes,"_rs.",n,"_ground_truth_effects.tsv")
        rownames(gt_effects)<-iso_res_rsem$transcript_id
        write.table(gt_effects,file=fout,sep="\t",quote=FALSE,row.names = T)
    }
    
    # Generate logTPM without batch effect
    
    gt_log_TPM_mx_noBatch <- gt_effects[,1:2] %*% t(design_mx)[1:2,]
    gt_log_TPM_mx_noBatch[nonZero_TPMs>0,] <- gt_log_TPM_mx_noBatch[nonZero_TPMs>0,] +
                                              rlaplace(length(gt_log_TPM_mx_noBatch[nonZero_TPMs>0,]), 
                                                       location=0,
                                                       scale=predict(res_vs_mean_fit,data.frame(mean_exp=c(gt_log_TPM_mx_noBatch[nonZero_TPMs>0,])))/sqrt(2))
    # Convert non-batch effect log counts to TPM  
    gt_TPM_mx_noBatch <- exp(gt_log_TPM_mx_noBatch)
    gt_TPM_mx_noBatch[nonZero_TPMs==0,] <- 0
    gt_TPM_mx_noBatch <- scale(gt_TPM_mx_noBatch, center = F, scale=colSums(gt_TPM_mx_noBatch)/1000000)
    colnames(gt_TPM_mx_noBatch) <- simulation_design$Sample

    # Write simulated TPM without batch effect
    for(sample in colnames(gt_TPM_mx_noBatch)){

        iso_res_rsem_noBatch <- iso_res_rsem                     
        iso_res_rsem_noBatch$TPM<-gt_TPM_mx_noBatch[,sample]
        rownames(iso_res_rsem_noBatch)<-NULL
        
        dout<-file.path(dir_results,sample)
        dir.create(dout, showWarnings = FALSE)
        fout<-paste0(dout,"/",0,"_rs.",n,"_rsem.isoforms.results")
        write.table(iso_res_rsem_noBatch,file = fout,sep="\t",quote=FALSE,row.names = F)
    }
    
}                     
                     
