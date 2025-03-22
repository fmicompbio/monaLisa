library(monaLisa)
lmrfile <- system.file("extdata", "LMRsESNPmerged.gr.rds", 
                       package = "monaLisa")
lmr <- readRDS(lmrfile)
set.seed(1)
lmrsel <- lmr[ sample(x = length(lmr), size = 10000, replace = FALSE) ]
bins <- bin(x = lmrsel$deltaMeth, binmode = "equalN", nElement = 800, 
            minAbsX = 0.3)
(gg <- plotBinDensity(lmrsel$deltaMeth, bins, legendPosition = "none") + 
        theme(axis.text = element_text(size = 14), 
              axis.title = element_text(size = 14)) + 
        labs(x = "atacPeaksChange"))
ggsave(plot = gg, filename = "man/figures/monaLisa_binning_small_ggplot.png")
