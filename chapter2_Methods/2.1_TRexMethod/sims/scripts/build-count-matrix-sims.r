suppressPackageStartupMessages({
	library("tximport")
})


# IPUTS

args = commandArgs(trailingOnly=TRUE)

inpdir<-args[1]
outdir<-args[2]
metadat<-args[3]

# MAIN 

print("Reading sample metadata...")
samp.info<-read.csv(metadat,row.names=1,header=TRUE,sep='\t')

#Writing txi objects
################################################

print("=============Input settings============")
print(inpdir)

#Reading sample metadata
print("Reading sample files:")
files<-file.path(inpdir,rownames(samp.info),"quant.sf") 
print(files)

print("=============Output settings============")
print("tximport object will be created at")
objpath<-paste(outdir,"/txiobject.raw.RData",sep="")
outfile<-paste(outdir,"/counts.raw.csv",sep="")
print(objpath)

#Build tximport object
print("Reading files from input directory...")


txi<-tximport(files,type="salmon",txOut=TRUE,countsFromAbundance="no",dropInfReps=T)
colnames(txi$counts)<-rownames(samp.info)

#Save tximport object

save(txi,file=objpath)
print(paste("Writing object to:",objpath))

salmon_counts<-txi$counts
colnames(salmon_counts)<-rownames(samp.info)
print(colSums(salmon_counts))

print(paste("Writing count matrix to:",outfile))
write.csv(salmon_counts,file=outfile)


