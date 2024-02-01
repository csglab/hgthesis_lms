suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(gridExtra)
})


# IPUTS

args = commandArgs(trailingOnly=TRUE)

cancer<-args[1]
txifile<-args[2]
metfile<-args[3]
resfile<-args[4]

# MAIN 

message("Loading data of cancer ",cancer)
load(txifile,verbose = FALSE)
meta <- read.csv(metfile,row.names = 1)

message("Extracting count data")
X <- txi$abundance %>% 
     t() %>%
     janitor::remove_constant()
rownames(X) <- colnames(txi$counts)

# Compute PCA

message("Computing PCA")
pca<-prcomp(X,center = T,scale = T,rank=10)
U <- pca$x 
U.df <- U %>% as.data.frame() %>% cbind(.,meta[rownames(U),])


p1<-data.frame(PC=factor(paste0("PC",1:10),levels=paste0("PC",1:10)),
           var_explained=(pca$sdev[1:10])^2/sum((pca$sdev[1:10])^2)) %>%
    ggplot(.,aes(PC,var_explained,group=1))+
    geom_col()+
    theme_bw()+
    labs(y="Variance explained")+
    theme(text = element_text(size=12),plot.title=element_text(face="bold",size=16,hjust=0.5))+
    labs(title=cancer)
p2<-ggplot(U.df,aes(x=PC1, y=PC2,color=condition)) +
     geom_point()+ 
     coord_equal()+ 
     theme_bw()+    
     theme(text = element_text(size=12),legend.position="top",plot.title=element_text(face="bold",size=16,hjust=0.5))+
    labs(title=cancer)+
    facet_wrap(~gender)
p3<-ggplot(U.df,aes(x=PC1, y=PC3,color=condition)) +
     geom_point()+ 
     coord_equal()+ 
     theme_bw()+    
     theme(text = element_text(size=12),legend.position="none",plot.title=element_text(face="bold",size=16,hjust=0.5))+
    labs(title=cancer)+
    facet_wrap(~gender)

# Caculate distances 

message("Calculating mahalanobis distance")

U.dist <- U.df %>%
             select(file_id,PC1:PC10,condition) %>% 
             group_by(condition) %>%
             tidyr::nest(data=file_id:PC10) %>%
             mutate(data = lapply(data,function(d){
                                         out<-d %>% tibble::column_to_rownames('file_id') %>% as.matrix()
                                         return(out)
                                        }
                                 )
                   ) %>%
             mutate(nsamples = lapply(data,nrow) %>% unlist()) %>%
             filter(nsamples>11) %>%
             mutate(dist = lapply(data,function(U){
                                         dist<-bigutilsr::dist_ogk(U) %>%
                                               as.data.frame() %>%
                                               rename("dist"=".") %>%
                                               mutate(file_id = rownames(U),
                                                      PC1 = U[,1],
                                                      PC2 = U[,2],
                                                      PC3 = U[,3],
                                                      pval = pchisq(dist, df = 10, lower.tail = FALSE),
                                                      padj = p.adjust(pval,method = "fdr"),
                                                      is.outlier.dist = padj<0.0001)
                                         return(dist)
                                       }
                                 )
                   ) %>%
             mutate(llof = lapply(data,function(U){
                             llof.df <- bigutilsr::LOF(U,seq_k = seq(5,nrow(U),by=5)[1:3]) %>%
                                        as.data.frame() %>%
                                        rename("llof"=".") %>%
                                        mutate(file_id = rownames(U))
                             tk<-bigutilsr::tukey_mc_up(llof.df$llof)
                             llof.df <- llof.df %>%
                                        mutate(is.outlier.lof = llof>tk)
                             return(llof.df)
                            })
                   ) %>%
             select(-data,-nsamples) %>%
             tidyr::unnest(dist,llof) %>%
             select(-file_id1)

p4<-ggplot(U.dist,aes(PC1,PC2,color=is.outlier.dist))+
    geom_point()+
    facet_wrap(~condition)+
    theme_bw()+
    theme(text = element_text(size=14),plot.title=element_text(face="bold",size=16,hjust=0.5))+
    labs(title=cancer)
p5<-ggplot(U.dist,aes(PC1,PC3,color=is.outlier.dist))+
    geom_point()+
    facet_wrap(~condition)+
    theme_bw()+
    theme(text = element_text(size=14),plot.title=element_text(face="bold",size=16,hjust=0.5))+
    labs(title=cancer)
p6<-ggplot(U.dist,aes(PC1,PC2,color=is.outlier.lof))+
    geom_point()+
    facet_wrap(~condition)+
    theme_bw()+
    theme(text = element_text(size=14),plot.title=element_text(face="bold",size=16,hjust=0.5))+
    labs(title=cancer)
p7<-ggplot(U.dist,aes(PC1,PC3,color=is.outlier.lof))+
    geom_point()+
    facet_wrap(~condition)+
    theme_bw()+
    theme(text = element_text(size=14),plot.title=element_text(face="bold",size=16,hjust=0.5))+
    labs(title=cancer)

# Write results

message("Generating visualizations")

pdf(resfile, onefile = TRUE,width = 6,height = 4)
plots<-list(p1,p2,p3,p4,p5,p6,p7)
for (i in seq(length(plots))) {
  print(plots[[i]])  
}
par(mfrow=c(1,2))
hist(U.dist$pval,main = "Histogram of pvalues")
hist(U.dist$padj,main = "Histogram of FDR")
message(cancer," - Flagged ",sum(U.dist$is.outlier.dist)," distance outliers out of ",length(U.dist$is.outlier.dist)," data points (",round(sum(U.dist$is.outlier.dist)/length(U.dist$is.outlier.dist),4)*100,"%)")
message(cancer," - Flagged ",sum(U.dist$is.outlier.lof)," LOF outliers out of ",length(U.dist$is.outlier.lof)," data points (",round(sum(U.dist$is.outlier.lof)/length(U.dist$is.outlier.lof),4)*100,"%)")
dev.off()

message("Writing flagged metadata")

metfile.out<-sub("csv","OLflag.csv",metfile)
outliers.lof <- U.dist %>% ungroup() %>% filter(is.outlier.lof == T) %>% select(file_id) %>% unlist() %>% as.character()
outliers.mdist <- U.dist %>% ungroup() %>% filter(is.outlier.dist == T) %>% select(file_id) %>% unlist() %>% as.character()
meta <- meta %>%
        mutate(is_outlier.lof = ifelse(file_id %in% outliers.lof,TRUE,FALSE),
               is_outlier.mdist = ifelse(file_id %in% outliers.mdist,TRUE,FALSE))
write.csv(meta,file=metfile.out,quote=T,row.names=T)

message("Finished successfuly!")