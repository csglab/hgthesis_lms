suppressPackageStartupMessages({
  library(getopt)
  library(optparse)
})
source("trex_helper.r")

option_list <- list(
  make_option(c("-s", "--sample"), type='character',dest="sampleid",
              help="Sample identifier used to name the output files"),
  make_option(c("-e", "--event"), type='character',dest="evtype",
              help="Name of AS event type to quantify. Must be one of A3,A5,AF,AL,MX,RI,SE,ALL"),
  make_option(c("-t", "--txiformat"), type='character',dest="format",default="tximeta",
              help="Output format of trex object. When set to 'tximport' it will also save a tximport object containing only the event-level counts."),
  make_option(c("-c","--counts"), type="character",dest = "txifile",
              help="Path to the input tximport object stored as an .RData file"),
  make_option(c("-m","--metadata"), type="character",dest = "metafile",
              help="Path to a two-column tab-separated file containing sample metadata (condition and replicate). The order of the rows must be the same as the columns in the matrices of the input tximport object."),
  make_option(c("-r","--refmap"), type="character",dest = "isomapfile",
              help="Path to the ioe file from SUPPA2"),
  make_option(c("-o","--output"),type="character",dest = "outdir",
              help="Directory to store the output of TRex.")
)

options<-parse_args(OptionParser(option_list = option_list))

message("############### TRex setup ##############")

message(" 1/3 Loading required packages...")

# Import libraries
suppressPackageStartupMessages({
  library(data.table)
  library(plyr)
  library(dplyr)
  library(lme4)
  library(stringr)
  library(SummarizedExperiment)
  library(tximport)
  library(hms)
  library(DESeq2)
})

###############################################
# Define functions
###############################################

## Load object 

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


################## MAIN ###############
#######################################

start_time <- Sys.time()

# Input parameters 

message(" 2/3 Reading inputs...",appendLF = T)

sampleid<-options$sampleid
outputfmt<-options$format
txifile<-options$txifile
metafile<-options$metafile
isomapfile<-options$isomapfile
outdir<-options$outdir
evtype<-options$evtype 

# Create directory

resdir<-outdir
if(!dir.exists(resdir)){dir.create(resdir)}

# Load data 
##########################################

message(" 3/3 Loading files...",appendLF = T)

## Salmon output
trex.txi<-load_tximport(txifile,metadata)

## Metadata
metadata<-fread(metafile,data.table = F,stringsAsFactors = F)
rownames(metadata)<-metadata[,1]

## SUPPA2 output
isomap<-fread(isomapfile,data.table = F,stringsAsFactors = F) %>% as.data.frame() %>% unique()
rownames(isomap)<-isomap$event_id
if(evtype!="ALL"){
    isomap<-isomap%>%
            filter(event_type==evtype)
}

# TRex main 
########################################

## Alternative isoforms

message("############# TRex algorithm ############")
message("---- Analyzing sample ",sampleid," ----")
message("---- Calculating counts for ",evtype," events ----")
message(" 1/10 Mapping alternative isoforms...",appendLF = T)

alternative<-strsplit(isomap$alternative_transcripts,",")
niso<-lapply(alternative,length)%>%unlist()
alt_isos<-data.frame(ensembl_transcript_id=unlist(alternative),
                      event_id=rep(isomap$event_id,times=niso))
alt_isos<-alt_isos[alt_isos$ensembl_transcript_id%in%rownames(trex.txi$counts),]

## Constitutive isoforms
message(" 2/10 Mapping constitutive isoforms...",appendLF = T)
constitutive<-strsplit(isomap$constitutive_transcripts,",")
niso<-lapply(constitutive,length)%>%unlist()
const_isos<-data.frame(ensembl_transcript_id=unlist(constitutive),
                       event_id=rep(isomap$event_id,times=niso))
const_isos<-const_isos[const_isos$ensembl_transcript_id%in%rownames(trex.txi$counts),]

## Aggregate counts
message(" 3/10 Aggregating isoform counts...",appendLF = T)
cts<-list(trex.txi$counts,trex.txi$abundance)
names(cts)<-c("counts","abundance")
cts_alt<-lapply(cts,trex_counts,alt_isos)
cts_const<-lapply(cts,trex_counts,const_isos)

## Match events 
message(" 4/10 Finding events with annotated alternative and constitutive isoforms...",appendLF = T)
matrices<-list()
for(s in names(cts)){
  
  ml<-suppressMessages(match_events(cts_alt[[s]],cts_const[[s]]))
  im<-ml[[1]]
  colnames(im)<-paste("A_",colnames(im),sep="")
  tm<-ml[[2]]
  colnames(tm)<-paste("C_",colnames(tm),sep="")
  
  matrices[[paste("A_",s,sep="")]]<-im
  matrices[[paste("C_",s,sep="")]]<-tm
}

## Calculate mean effective length 
message(" 5/10 Calculating mean effective length per event...",appendLF = T)
colnames(trex.txi$length)<-rownames(metadata)

alt_isos<-alt_isos[alt_isos$event_id%in% rownames(matrices$A_counts),]
const_isos<-const_isos[const_isos$event_id%in% rownames(matrices$A_counts),]

alt_length<-suppressMessages(trex_means(trex.txi$length[alt_isos$ensembl_transcript_id,],alt_isos))
total_length<-suppressMessages(trex_means(trex.txi$length,unique(rbind(alt_isos,const_isos))))

ml<-suppressMessages(match_events(alt_length,total_length))
im<-ml[[1]]
colnames(im)<-paste("A_",colnames(im),sep="")
tm<-ml[[2]]
colnames(tm)<-paste("C_",colnames(tm),sep="")

# Create tximport object
message(" 6/10 Creating trex-tximport object with event counts...",appendLF = T)

trex.txi$counts<-cbind(matrices[["A_counts"]],matrices[["C_counts"]])
trex.txi$abundance<-cbind(matrices[["A_abundance"]],matrices[["C_abundance"]])
trex.txi$length<-cbind(im,tm)

message(" 7/10 Correcting count estimates...",appendLF = T)

trex.txi<-correct_event_counts(trex.txi)

rm(im,tm,total_length,alt_isos,const_isos,alt_length,ml,matrices,constitutive,alternative,cts,niso,cts_alt,cts_const)

# TRex main - tximeta object
########################################

# Re-format inferential replicates 

message(" 8/10 Building data frames with gene and sample-level annotations...",appendLF = T)

## Rows - AS events
rowinfo<-isomap[rownames(trex.txi$counts),]

## Columns - samples and count type

colinfo<-metadata
colinfo<-rbind(colinfo,metadata)

colinfo$sample<-rep(paste0("s",1:nrow(metadata)),2)
rownames(colinfo)<-c(paste("A",metadata$Sample,sep="_"),
                     paste("C",metadata$Sample,sep="_"))

colinfo$trex_count<-strsplit(colnames(trex.txi$counts),"_")%>%
                    lapply(.,function(n){n[[1]]})%>%
                    unlist()%>%
                    factor(.,levels=c("C","A"))

## Counts
mats<-trex.txi[c("abundance","counts","length")]
mats<-lapply(mats, function(m,...){
  n<-round(as.matrix(m),digits = 0)
  mode(n)<-"numeric"
  rownames(n)<-rownames(rowinfo)
  colnames(n)<-rownames(colinfo)
  return(n)
})


message(" 9/10 Creating baseline oject with event counts...",appendLF = T)
trex.meta<-SummarizedExperiment(assays = mats,
                                rowData = rowinfo,
                                colData = colinfo,
                                metadata = "lengthScaledTPM")

message(" 10/10 Writing trex output...",appendLF = T)

outfile<-paste(outdir,"/",evtype,".eventCounts.RData",sep="")
save(trex.meta,file=outfile)

outfile<-paste(outdir,"/",evtype,".eventCounts.tximport.RData",sep="")
save(trex.txi,file=outfile)

end_time <- Sys.time()
run_time <- difftime(end_time, start_time)
message(paste("Finished succesfully!\nEllapsed time:",as_hms(run_time)),appendLF = T)