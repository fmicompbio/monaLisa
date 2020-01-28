#' @title Sample random regions.
#'
#' @description Sample random regions from the mappable parts of the genome with
#'      a given fraction from CpG islands.
#'
#' @param N Number of regions to sample.
#' @param regWidth Region width.
#' @param fractionCGI fraction of regions to be sampled from CpG islands. Not implemented yet.
#' @param mappableRegions A unstranded GRanges object of the mappable parts of the genome. This can be a GRanges of
#'                        of the entire genome with one range from start to end of each chromosome if mappability
#'                        is not an issue.
#' @param seed Seed for seeding the random number generator.
#'                  
#' 
#' @details to be filled in. 
#'
#' @return A GRanges object with randomly sampled mappable regions of width \code{regWidth}
#'         with \code{fractionCGI} coming from CpG islands.
#'
#' @export
sample_random_regions <- function(N = 100L, regWidth = 200L, fractionCGI = 0.5,
																	mappableRegions = NULL, seed = 123){
	
	stopifnot(exprs = {
		length(N) == 1
		is.numeric(N)
		length(regWidth) == 1
		is.numeric(regWidth)
		length(fractionCGI) == 1
		is.numeric(fractionCGI)
		(fractionCGI >= 0) & (fractionCGI <= 1)
		class(mappableRegions) == "GRanges"
		length(mappableRegions) > 0
		sum(width(mappableRegions) >= regWidth) > 0
		unique(strand(mappableRegions)) == "*"
		length(seed) == 1
		is.numeric(seed)
	})
	
	#in case these parameters are not integers, round to integer
	N <- round(N)
	regWidth <- round(regWidth)
	seed <- round(seed)
	
	#reduce to ranges that are larger than or equal to width and adjust
	#end position such that the object only contains
	#the start positions that when extended by width are 
	#entirely within the mappable parts of the genome
	gr <- mappableRegions[width(mappableRegions) >= regWidth]
	end(gr) <- end(gr) - regWidth + 1
	
	#first sample regions according to the number of positions
	#then sample position uniformly within each region
	
	#probability of sampling each region
	widths <- width(gr)
	ps <- widths/sum(widths)
	
	#sample regions
	set.seed(seed)
	indx <- sample(1:length(gr), size = N, prob = ps, replace=TRUE)
	
	#for each range, sample 1 position uniformly
	gr.sampled <- gr[indx]
	st.sampled <- start(gr.sampled)
	end.sampled <- end(gr.sampled)
	chr.sampled <- seqnames(gr.sampled)
	
	pos.sampled <- unlist(lapply(1:length(gr.sampled), function(i){
		sample(st.sampled[i]:end.sampled[i], 1)
	}))
	gr.sampled <- GRanges(seqnames = chr.sampled,
												ranges=IRanges(start = pos.sampled, width = regWidth),
												seqlengths = seqlengths(gr))
	gr.sampled[order(as.character(seqnames(gr.sampled)), start(gr.sampled))]
	
}
