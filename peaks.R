library(ChIPpeakAnno)

ZF19 <- toGRanges("macs2/8406-ZF-19_GCTACGCT-ATAGAGAG_S26_summits.bed", format="BED", header =FALSE)
ZF22 <- toGRanges("macs2/8406-ZF-22_GTAGAGGA-ATAGAGAG_S27_summits.bed", format="BED", header =FALSE)
ol <- findOverlapsOfPeaks(ZF19, ZF22)
ol <- addMetadata(ol, colNames="score", FUN=mean)
ol$peaklist[["ZF19///ZF22"]][1:2]

pdf(file="ZF_19_22.pdf", width =12) 
makeVennDiagram(ol, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"))
dev.off() 

