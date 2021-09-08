.binomEnrichmentTest <- function(matchCountBg, totalWeightBg, matchCountFg,
                                 totalWeightFg, verbose) {
    if (verbose) {
        message("using binomial test to calculate ",
                "log(p-values) for enrichments")
    }
    
    prob <- matchCountBg / totalWeightBg
    minProb <- 1 / totalWeightBg
    maxProb <- (totalWeightBg - 1) / totalWeightBg
    if (any(i <- (prob < minProb))) {
        # warning("some background match probabilities are below ",
        #         "minProb (for example when there were zero hits) ",
        #         "and will be given a value of minProb=1/totalWeightBg")
        prob[i] <- minProb
    }
    if (any(i <- (prob > maxProb))) {
        # warning("some match probabilities a above",
        #         "maxProb (for example when all sequences had hits) ",
        #         "and will be given a value of ",
        #         "maxProb=(totalWeightBg-1)/totalWeightBg")
        prob[i] <- maxProb
    }
    
    pbinom(q = matchCountFg - 1,
           size = totalWeightFg,
           prob = prob, lower.tail = FALSE, log.p = TRUE)
}

.fisherEnrichmentTest <- function(matchCountBg, totalWeightBg, matchCountFg,
                                  totalWeightFg, verbose) {
    if (verbose) {
        message("using Fisher's exact test (one-sided) to calculate ",
                "log(p-values) for enrichments")
    }
    
    # contingency table per sequence for Fisher's exact test (rounded to integer):
    #              withHit  noHit
    #   foreground    x       y
    #   background    z       w
    #
    logP <- log(vapply(structure(seq_along(matchCountFg),
                                 names = names(matchCountFg)),
                       function(i) {
                           ctab <- rbind(c(matchCountFg[i],
                                           totalWeightFg - matchCountFg[i]),
                                         c(matchCountBg[i],
                                           totalWeightBg - matchCountBg[i]))
                           ctab <- round(ctab)
                           fisher.test(x = ctab, alternative = "greater")$p.value
                       }, FUN.VALUE = numeric(1)))
}

.calcPearsonResiduals <- function(matchCountBg, totalWeightBg, matchCountFg,
                                  totalWeightFg) {
    obsTF <- matchCountFg
    expTF <- totalWeightFg * (matchCountFg + matchCountBg) /
        (totalWeightFg + totalWeightBg)
    N <- totalWeightFg + totalWeightBg
    enr <- (obsTF - expTF) / sqrt(expTF * (1 - totalWeightFg / N) *
                                      (1 - (matchCountFg +
                                                matchCountBg) / N))
    enr[is.na(enr)] <- 0
    enr
}

.calcExpFg <- function(matchCountBg, totalWeightBg, matchCountFg, 
                       totalWeightFg) {
    totalWeightFg * (matchCountFg + matchCountBg) / 
        (totalWeightFg + totalWeightBg)
}

.calcLog2Enr <- function(matchCountBg, totalWeightBg, matchCountFg, 
                         totalWeightFg, pseudocount) {
    minTot <- min(totalWeightFg, totalWeightBg)
    normFg <- log2(matchCountFg/totalWeightFg * minTot + pseudocount) 
    normBg <- log2(matchCountBg/totalWeightBg * minTot + pseudocount)
    log2enr <- normFg - normBg
    log2enr
    
    # D <- enrich1[, c("sumForegroundWgtWithHits", "sumBackgroundWgtWithHits")]
    # nTot <- unlist(enrich1[1, c("totalWgtForeground", "totalWgtBackground")])
    # D.norm <- t(t(D) / nTot * min(nTot))
    # DL <- log2(D.norm + pseudocount.log2enr)
    # log2enr <- DL[, 1] - DL[, 2]
    # log2enr
}

# enr <- .calcPearsonResiduals(matchCountBg = enrich1[, "sumBackgroundWgtWithHits"],
#                              totalWeightBg = enrich1[, "totalWgtBackground"],
#                              matchCountFg = enrich1[, "sumForegroundWgtWithHits"],
#                              totalWeightFg = enrich1[, "totalWgtForeground"])
# names(enr) <- enrich1[, "motifName"]



