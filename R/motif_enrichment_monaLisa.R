#' @importFrom Biostrings alphabetFrequency oligonucleotideFrequency
NULL

#' 
#' Functions folowing/inspired by HOMER to do motif enrichment analysis 
#' 
#'  
#' 
#' 
#' 
#' 


#' @title Filter Bad Sequences
#'
#' @description We filter sequences similarly to how HOMER does it. Namely, sequences with more
#'   than 70% (default) N are removed.
#'   
#' @param seqs object of class \code{DNAStringSet} containing all the sequences (foreGround and backGround)
#' @param frac the fraction of Ns allowed in a sequence (default set to 0.7)
#'   
#' @return object of class \code{DNAStringSet} containing the sequences that passed the filter.
#'   
#' @importFrom Biostrings alphabetFrequency
#' 
#' @export
filterSeqs <- function(seqs=NULL, frac=0.7) {
  
  # checks
  if(!class(seqs)=="DNAStringSet"){stop("'seqs' mus be of class DNAStringSet")}
  if(is.null(seqs)){stop("'seqs' in NULL")}
  
  # fraction of Ns per sequence
  fracN <- Biostrings::alphabetFrequency(seqs)[, "N"] / lengths(seqs)
  
  # remove sequences with fraction > frac
  w <- which(fracN>frac)
  if(!isEmpty(w)){
    seqs <- seqs[-w]
  }
  
  # return seqs passing filter
  seqs
  
}


#' @title Get GC weights for background regions 
#'
#' @description In this implementation we still follow HOMER. But this can be re-implemented in 
#'   a more continuous fashion. Each sequence is put in a specific GC bin depending on its GC content. 
#'   Then, for each GC bin, the number of foreGround and backGround sequences in that bin is calculated.
#'   Weights are calculated for the backGround sequences in bin i as follows:
#'   weight_i = (number_fg_seqs_i/number_bg_seqs_i)*(number_bg_seqs_total/number_fg_seqs_total)
#'   
#' @param foreGround_seq object of class \code{DNAStringSet} containing all foreGround sequences
#' @param backGround_seq object of class \code{DNAStringSet} containing all backGround sequences
#'   
#' @return a list containing:
#'   \itemize{
#'     \item{sequenceWeights}{: a \code{dataframe} comtaining the sequence GC content, GC bins
#'       they were assigned to, and the weight to correct for GC differences between foreGround
#'       and backGround sequences}
#'     \item{sequenceNucleotides}{: a \code{DNAStringSet} object containing the raw sequences}
#'   }
#' 
#' @importFrom Biostrings oligonucleotideFrequency  
#'   
#' @export
getGCweight <- function(foreGround_seq=NULL, backGround_seq=NULL) {
  
  # checks
  if(is.null(foreGround_seq)){stop("'foreGround_seq' is NULL")}
  if(is.null(backGround_seq)){stop("'backGround_seq' is NULL")}
  if(!class(foreGround_seq)=="DNAStringSet"){stop("'foreGround_seq' must be a DNAStringSet")}
  if(!class(backGround_seq)=="DNAStringSet"){stop("'backGround_seq' must be a DNAStringSet")}

  # calculate GC fraction for each sequence
  # ... foreGround
  fMono <- Biostrings::oligonucleotideFrequency(foreGround_seq, width=1, as.prob=TRUE)
  fDi <- Biostrings::oligonucleotideFrequency(foreGround_seq, width=2, as.prob=TRUE)
  fg_frac <- fMono[,"G"]+fMono[,"C"]
  # ... backGround
  fMono <- Biostrings::oligonucleotideFrequency(backGround_seq, width=1, as.prob=TRUE)
  fDi <- Biostrings::oligonucleotideFrequency(backGround_seq, width=2, as.prob=TRUE)
  bg_frac <- fMono[,"G"]+fMono[,"C"]
  
  # HOMER's GC breaks/bins
  gc_breaks <- c(0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8)
  
  # assign each sequence to a GC bin (separately for foreGround and backGround)
  GC_bin_fg <- findInterval(x = fg_frac, vec = gc_breaks)
  GC_bin_bg <- findInterval(x = bg_frac, vec = gc_breaks)
  
  # keep bins that have at least 1 foreGround and 1 backGround sequence
  bins <- unique(c(GC_bin_fg, GC_bin_bg))
  keep <- bins%in%unique(GC_bin_fg) & bins%in%unique(GC_bin_bg)
  bins <- bins[keep]
  
  # keep sequences belonging to these bins
  # ... foreGround
  keep_fg <- GC_bin_fg%in%bins
  GC_bin_fg <- GC_bin_fg[keep_fg]
  foreGround_seq <- foreGround_seq[keep_fg]
  fg_frac <- fg_frac[keep_fg]
  # ... backGround
  keep_bg <- GC_bin_bg%in%bins
  GC_bin_bg <- GC_bin_bg[keep_bg]
  backGround_seq <- backGround_seq[keep_bg]
  bg_frac <- bg_frac[keep_bg]
  
  # total number of foreGround and backGround sequences
  total_fg <- length(GC_bin_fg)
  total_bg <- length(GC_bin_bg)
  
  # calculate GC weight per bin
  weight_per_bin <- sapply(bins, function(b){
    n_fg_b <- sum(GC_bin_fg%in%b) # number of fg seqs in b
    n_bg_b <- sum(GC_bin_bg%in%b) # number of bg seqs in b
    (n_fg_b/n_bg_b)*(total_bg/total_fg)
  })
  
  # assign calculated GC weight to each backGround sequence (foreGround get a weight of 1)
  weights_fg <- rep(1, length(GC_bin_fg))
  weights_bg <- rep(1, length(GC_bin_bg))
  for(i in 1:length(bins)){
    b <- bins[i]
    w <- weight_per_bin[i]
    weights_bg[GC_bin_bg%in%b] <- w
  }
  
  # return list of sequences (that made it) and their weights
  df <- data.frame(foreGround=c(rep(1, length(weights_fg)), rep(0, length(weights_bg))), 
                   GCfraction=c(fg_frac, bg_frac), 
                   GCbin=c(GC_bin_fg, GC_bin_bg), 
                   GCweight=c(weights_fg, weights_bg))
  seqs <- c(foreGround_seq, backGround_seq)
  
  list(sequenceWeights=df, sequenceNucleotides=seqs)
  
}
























