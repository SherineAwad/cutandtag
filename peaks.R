library(ChIPpeakAnno)

ZF19 <- toGRanges("macs2/8406-ZF-19_GCTACGCT-ATAGAGAG_S26_summits.bed", format="BED", header =FALSE)
ZF22 <- toGRanges("macs2/8406-ZF-22_GTAGAGGA-ATAGAGAG_S27_summits.bed", format="BED", header =FALSE)
ZF7 <- toGRanges("macs2/8406-ZF-7_CAGAGAGG-GCGATCTA_S17_summits.bed", format="BED", header =FALSE)
ZF6 <- toGRanges("macs2/8406-ZF-6_CTCTCTAC-GCGATCTA_S16_summits.bed", format="BED", header =FALSE)
ZF1 <- toGRanges("macs2/8406-ZF-1_TAAGGCGA-GCGATCTA_S12_summits.bed", format="BED", header =FALSE)
ZF22 <- toGRanges("macs2/8406-ZF-22_GTAGAGGA-ATAGAGAG_S27_summits.bed", format="BED", header =FALSE)
ZF19 <- toGRanges("macs2/8406-ZF-19_GCTACGCT-ATAGAGAG_S26_summits.bed", format="BED", header =FALSE)
ZF29 <- toGRanges("macs2/8406-ZF-29_CTCTCTAC-AGAGGATA_S33_summits.bed", format="BED", header =FALSE)
ZF4 <- toGRanges("macs2/8406-ZF-4_TCCTGAGC-GCGATCTA_S14_summits.bed", format="BED", header =FALSE)
ZF23 <- toGRanges("macs2/8406-ZF-23_TAAGGCGA-AGAGGATA_S28_summits.bed", format="BED", header =FALSE)
ZF16 <- toGRanges("macs2/8406-ZF-16_TAGGCATG-ATAGAGAG_S23_summits.bed", format="BED", header =FALSE)
#ZF14 <- toGRanges("macs2/8406-ZF-14_TCCTGAGC-ATAGAGAG_S21_summits.bed", format="BED", header =FALSE)
#ZF17 <- toGRanges("macs2/8406-ZF-17_CTCTCTAC-ATAGAGAG_S24_summits.bed", format="BED", header =FALSE)
#ZF25 <- toGRanges("macs2/8406-ZF-25_AGGCAGAA-AGAGGATA_S30_summits.bed", format="BED", header =FALSE)



ol1 <- findOverlapsOfPeaks(ZF19, ZF22)
ol1 <- addMetadata(ol1, colNames="score", FUN=mean)
ol1$peaklist[["ZF19///ZF22"]][1:2]

ol2 <- findOverlapsOfPeaks(ZF16, ZF22)
ol2 <- addMetadata(ol2, colNames="score", FUN=mean)
ol2$peaklist[["ZF16///ZF22"]][1:2]

ol3 <- findOverlapsOfPeaks(ZF16, ZF19)
ol3 <- addMetadata(ol3, colNames="score", FUN=mean)
ol3$peaklist[["ZF16///ZF19"]][1:2]

ol4 <- findOverlapsOfPeaks(ZF4, ZF7)
ol4 <- addMetadata(ol4, colNames="score", FUN=mean)
ol4$peaklist[["ZF4///ZF7"]][1:2]

ol5 <- findOverlapsOfPeaks(ZF23, ZF29)
ol5 <- addMetadata(ol5, colNames="score", FUN=mean)
ol5$peaklist[["ZF23///ZF29"]][1:2]

ol6 <- findOverlapsOfPeaks(ZF1, ZF6)
ol6 <- addMetadata(ol6, colNames="score", FUN=mean)
ol6$peaklist[["ZF1///ZF6"]][1:2]

pdf(file="peaks_intersection.pdf", width =12) 
makeVennDiagram(ol1, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"))

makeVennDiagram(ol2, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"))


makeVennDiagram(ol3, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"))

makeVennDiagram(ol4, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"))

makeVennDiagram(ol5, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"))

makeVennDiagram(ol6, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"))

dev.off()




