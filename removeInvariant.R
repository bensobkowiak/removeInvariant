#' Remove invariant sites and excess ambiguity from FASTA file
#' @param fastafile Sequence file in FASTA format
#' @param indexfile Positions of SNPs in .txt file, single column with header (optional)
#' @param prefix Output file prefix
#' @param missing Code for missing sites (default = ? and -)
#' @param keepAmb If TRUE, include ambiguity as variant code, otherwise just consider A, C, G, T (default = FALSE)
#' @param excessAmbPerc Maximum percentage of ambiguous calls allowed to retain site, set to 0 to keep all sites (default = 10)
#' @return new FASTA file and position index file after processing
#' @export


removeInvariant<-function(fastafile,indexfile=NULL,prefix="output",missing=c("?","-"),keepAmb=FALSE,excessAmbPerc=10){
  require(seqinr)
  options(stringsAsFactors = F)
  fasta<-read.fasta(fastafile,forceDNAtolower = F)
  fastnames<-names(fasta)
  if (is.null(indexfile)){
    pos<-data.frame(Position = 1:length(fasta[[1]]))
  } else {
    pos<-read.table(indexfile,header = T)
  }
  nucs<-c("A","C","G","T","a","c","g","t")
  invariant<-numeric()
  for (i in 1:length(fasta[[1]])){
    nucls<-unique(as.character(sapply(fasta, "[[", i)))
    nucls<-nucls[!nucls %in% missing]
    if(!keepAmb){
      nucls<-nucls[nucls %in% nucs]
    }
    if(length(nucls)<2){
      invariant<-c(invariant,i)
    }
  }
  if (length(invariant)>0){
    fasta<-lapply(1:length(fasta), function(x){fasta[[x]][-invariant]})
    pos<-pos[-invariant,]
  }
  if (excessAmbPerc>0){
    remove<-numeric()
    for (i in 1:length(fasta[[1]])){
      sites<-as.character(unlist(lapply(fasta, "[[", i)))
      if (length(which(sites == "N" | sites == "?"))>length(fasta)*(excessAmbPerc/100)){
        remove<-c(remove,i)
      }
    }
    if (length(remove)>0){
      fasta<-lapply(1:length(fasta), function(x){fasta[[x]][-remove]})
      pos<-pos[-remove,]
    }
  }
  write.fasta(newfasta,fastnames,paste0(prefix,".fasta"),open = "w")
  write.table(pos,paste0(prefix,"_index.txt"),quote=F,row.names=F)
}

