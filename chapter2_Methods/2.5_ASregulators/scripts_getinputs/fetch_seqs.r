library(dplyr)
library(GenomicRanges)
library(Biostrings)
library(parallel)

# Load data
message("Loading coordinate GRanges...")
ioe.gr.ev<-readRDS(file = "../data/objects/event_cords_gr.RDS")
message("Building human genome object...")
hg<-BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# Extract sequences
message("Fetching sequences...")
ioe.gr.ev <- sortSeqlevels(ioe.gr.ev)
ioe.gr.ev <- sort(ioe.gr.ev)

ws<-c(51,101,201)
for(w in ws){
    
    message("Window size ",w-1)
    wgr<-ioe.gr.ev[ioe.gr.ev@ranges@width==w]
    wgr.strands<-as.character(strand(wgr))
    
    # Extract DNA sequences
    message("extracting DNA sequences...")
    dna_seqs<-BSgenome::getSeq(hg, wgr)
    
    message("splitting sequences by strand...")
    dna_seqs_neg<-dna_seqs[wgr.strands=="-"]

    # Converte to mRNA sense 
    message("Converting to mRNA sense...")
    rna_seqs <- reverseComplement(dna_seqs_neg)
    
    # Combine seqs objects
    message("Merging back sequences from both strands...")
    seqs <- DNAStringSet(c(dna_seqs[wgr.strands=="+"],
                           rna_seqs))
    names(seqs)<-c(wgr$seq_id[wgr.strands=="+"],wgr$seq_id[wgr.strands=="-"])
    
    # Create fasta file
    message("Writting fasta file...")
    Biostrings::writeXStringSet(seqs, paste0('../data/fastas/spliceSites_windowSize_',w-1,"_up_and_down",".fasta"))
}
message("Finished successfully!")