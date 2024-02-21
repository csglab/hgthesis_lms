# TRex functions


load_tximport<-function(txifile,metadata){
  
  #############################################
  # This will load a tximport object and annotate the columns
  # with the information in the data.frame `metadata`
  # the order of the columns in the tximport object must be
  # the same as the rows in the metadata table. 
  #############################################
  
  load(txifile,verbose=F)
  trex.txi<-txi
  
  # Add names to columns from metadata
  colnames(trex.txi$counts)<-rownames(metadata)
  colnames(trex.txi$abundance)<-rownames(metadata)
  colnames(trex.txi$length)<-rownames(metadata)
  return(trex.txi)
}

## Convert transcript counts to AS event counts

trex_counts<-function(counts,trex.map,...){
  
  trex_sum<-aggregate(counts[trex.map$ensembl_transcript_id,],
                      list(event_id = trex.map$event_id),
                      function(x){sum(x,na.rm = F)}) %>% 
    tibble::column_to_rownames("event_id")
  return(trex_sum)
}

trex_means<-function(lengths,trex.map,...){
  
  trex_mean<-aggregate(lengths[trex.map$ensembl_transcript_id,],
                       list(event_id = trex.map$event_id),
                       function(x){mean(x,na.rm = F)}) %>% 
    tibble::column_to_rownames("event_id")
  return(trex_mean)
}

## Match alternative and constitutive isoforms for all events
  
match_events<-function(alt,const){
  
  colnames(alt)<-paste("alt",colnames(alt),sep='_')
  colnames(const)<-paste("const",colnames(const),sep='_')
  
  alt$event_id<-rownames(alt)
  const$event_id<-rownames(const)
  
  match<-join(alt,const,type = 'inner')%>%
         tibble::column_to_rownames('event_id')
  match[is.na(match)]<-0
  
  m_alt <- match[,grep('alt',colnames(match))]
  colnames(m_alt)<-sub('alt_','',colnames(m_alt))
  
  m_const <- match[,grep('const',colnames(match))]
  colnames(m_const)<-sub('const_','',colnames(m_const))
  
  mats<-list(m_alt,m_const)
  
  return(mats)
}

# Correct event counts for the average isoform length of all isoforms in the same count type group (A or C)

correct_event_counts<-function(trex){

    # and for the average library size of each count type (A or C). 
    # Returns tximprot like object with updated 'counts','abundance','length', and 'countsFromAbundance' slots
    
    mats<-list()
    for(tc in c("A","C")){
        cols <- grepl(paste0(tc,"_"),colnames(trex$counts))

        # Retrieve counts and length matrices
        countsMat <- trex$counts[,cols]
        lengthMat <- trex$length[,cols] %>% 
                     apply(.,1,function(row){
                         # replace missing lengths with average length across all samples
                         mr<-mean(row[row!=0],na.rm=T)
                         row[row==0|is.na(row)]<-mr
                         return(row)
                     }) %>% t()

        # Re-esimate abundances for the trex count type
        tpmMat <- countsMat/lengthMat
        tpmSum <- colSums(tpmMat)
        abundanceMat <- t(t(tpmMat)/tpmSum*1e6)

        # Estimate counts from abundances per trex count type
        eCountsMat <- abundanceMat * rowMeans(lengthMat) # Scale to average length over all samples
        scaleFactor <- colSums(countsMat)/colSums(eCountsMat)
        eCountsMat <- t(t(eCountsMat) * scaleFactor)

        mats[[tc]] <- list("counts"=eCountsMat,"abundance"=abundanceMat,"length"=lengthMat)
    }
    trex$counts<-cbind(mats[["A"]][["counts"]],mats[["C"]][["counts"]])
    trex$abundance<-cbind(mats[["A"]][["abundance"]],mats[["C"]][["abundance"]])
    trex$length<-cbind(mats[["A"]][["length"]],mats[["C"]][["length"]])
    trex$countsFromAbundance<-"lengthScaledTPM"
    
    return(trex)
}

## Compute delta PSI values across replicates
  
delta.psi<-function(psi,metadata,length=F){
  
  cond<-metadata[colnames(psi),"condition"] %>% levels()
  reps<-metadata[colnames(psi),"replicate"] %>% levels()
  
  delta.reps<-data.frame(row.names = 1:nrow(psi))
  for(rep in reps){
    delta.reps[,paste("dPSA_",rep,sep="")]<-psi[,paste(cond[2],rep,sep="_")]-psi[,paste(cond[1],rep,sep="_")]
  }
  delta.reps$mean_dPSI<-rowMeans(delta.reps)
  delta.reps$sdev_dPSI<-apply(delta.reps[,1:2],1,sd)
  rownames(delta.reps)<-rownames(psi)
  
  return(delta.reps)
}