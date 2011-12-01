library(Biostrings)
library(GenomicRanges)
simulate.large.indels <- function(genome.seq, n.ins, n.dels=n.ins,
                                 # pos.ins = NULL, pos.del = NULL,
                                  length.ins = rpois(n.ins, 1e5),
                                  length.dels = rpois(n.dels, 1e5),
                                  insert.seq = NULL,
                                  strand.prob=c("+"=0.5,"-"=0.5),
                                  uniform.by.base=TRUE
                                
                                  ins.source.sequence=NULL,
                                  ins.source.regions=NULL,
                                  ins.masked.regions= NULL,
                                  del.masked.regions=ins.masked.regions,
                                  make.diploid=TRUE,
                                  homozygous.loci.prob = 0,
                                  sex.chromosomes=c("X","Y")
                                  
                                       
                                  ) {
 border.delete.policy <- match.arg(border.delete.policy) 
 stopifnot(inherits(genome.seq, "XStringSet"))
 if (! is.null(insert.seq))
   stopifnot(inherits(insert.seq, "XStringSet"))
 stopifnot(is.integer(n.ins) && n.ins >= 0)
 stopifnot(length(strand.prob==2) && names(strand.prob) %in% c("+","-"))

 }

insert.at <- function(sequence, insert, pos, strand="+") {
  insert <- DNAString(insert)
  if (strand=="-") insert <- reverseComplement(insert)
  subseq(sequence, start=pos, width=0) <- insert
  return (sequence)
}

rand.chrom <- function(n, genome, uniform.by.base=TRUE){
 prob <- if (uniform.by.base) width(genome)/sum(width(genome))
 else NULL  
  sample (n,x=1:length(genome), prob=prob, replace=TRUE )
}

rand.delete.loci <- function(n, genome, lengths, masks, uniform) {

  chrs <- rand.chrom(n, genome, uniform)
  starts <- sapply(width(genome[chrs]), function(x){sample(1:x, size=1)})
  gr <- GRanges(seqnames=Rle(names(genome[chrs])), ranges=IRanges(start=starts, width=lengths), mutation.type=rep("deletion", n)) 
  if (! is.null(masks)) {
    gr <- redraw.on.overlap(gr, masks, genome, chrs)
  }
  gr
  
}

rand.insert.loci <- function(n, genome, lengths, masks, strand.prob, uniform) {
  chrs <- rand.chrom(n, genome, uniform)
  strands <- sample(c("+","-"), size=n, replace=T, prob=strand.prob[c("+","-")])
  starts <- sapply(width(genome[chrs]), function(x){sample(1:x, size=1)})
  gr <- GRanges(seqnames=Rle(names(genome[chrs])), ranges=IRanges(start=starts, width=lengths),
                strand=Rle(strands),  mutation.type=rep("insertion", n))

  ## check if any range is in overlap with masks
  ## draw a new start position for them
  ## this is not a very good implementation, if there are large masks this will take a while
  if (! is.null(masks)) {
    gr <- redraw.on.overlap(gr, masks, genome, chrs)
  }
  gr
  
}

redraw.on.overlap <- function(gr, masks, genome, chrs) {
  while(any (maskovl <- (countOverlaps(gr,masks) > 0) ) ) {
    new.starts <-  sapply(width(genome[chrs[maskovl]]), function(x){sample(1:x, size=1)})
    start(gr[maskovl]) <- new.starts
  }
  gr
}
