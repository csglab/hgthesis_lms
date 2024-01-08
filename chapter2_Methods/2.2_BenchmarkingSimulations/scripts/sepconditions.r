#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print("Reading input files")
cts<-args[1]
psi<-args[2]
out.cst<-args[3]
out.psi<-args[4]


cts.dat<-read.table(cts,header=T,row.names=1)
psi.dat<-read.table(psi,header=T,row.names=1)

samples<-unique(unlist(lapply(strsplit(colnames(cts.dat),"_"),function(x){ paste(x[[1]],x[[2]],sep="_")})))

print(paste("Detected",length(samples),"samples:"))
print(samples)

# Write files

print("Writing individual sample files:")
for(sample in samples){

	print(sample)
	sname<-sub("\\.","_",sample)

	#Counts
	
	data<-cts.dat[,grep(sample,colnames(cts.dat))]
	colnames(data)<-paste(sname,1:ncol(data),sep='')

	#fout<-paste(out.cst,"/",sname,".counts.tsv",sep="")
	fout<-paste(out.cst,sname,".counts.tsv",sep="")
	write.table(data,fout,row.names=T,quote=F,sep='\t')

	#PSI 

	data<-psi.dat[,grep(sample,colnames(psi.dat))]
	colnames(data)<-paste(sname,1:ncol(data),sep='')

	fout<-paste(out.psi,sname,".SE.events.psi",sep="")
	write.table(data,fout,row.names=T,quote=F,sep='\t')

}