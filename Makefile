ZF-25.out:
	python count_genes.py 8406-ZF-7_CAGAGAGG-GCGATCTA_S17.annotatednarrowpeaks ZF-7.out 11
	python count_genes.py 8406-ZF-6_CTCTCTAC-GCGATCTA_S16.annotatednarrowpeaks ZF-6.out 11
	python count_genes.py 8406-ZF-1_TAAGGCGA-GCGATCTA_S12.annotatednarrowpeaks ZF-1.out 11
	python count_genes.py 8406-ZF-22_GTAGAGGA-ATAGAGAG_S27.annotatednarrowpeaks ZF-22.out 11
	python count_genes.py 8406-ZF-19_GCTACGCT-ATAGAGAG_S26.annotatednarrowpeaks ZF-19.out 11
	python count_genes.py 8406-ZF-29_CTCTCTAC-AGAGGATA_S33.annotatednarrowpeaks ZF-29.out 11
	python count_genes.py 8406-ZF-4_TCCTGAGC-GCGATCTA_S14.annotatednarrowpeaks ZF-4.out 11
	python count_genes.py 8406-ZF-23_TAAGGCGA-AGAGGATA_S28.annotatednarrowpeaks ZF-23.out 11
	python count_genes.py 8406-ZF-16_TAGGCATG-ATAGAGAG_S23.annotatednarrowpeaks ZF-16.out 11
	python count_genes.py 8406-ZF-14_TCCTGAGC-ATAGAGAG_S21.annotatednarrowpeaks ZF-14.out 11
	python count_genes.py 8406-ZF-17_CTCTCTAC-ATAGAGAG_S24.annotatednarrowpeaks ZF-17.out 11 
	python count_genes.py 8406-ZF-25_AGGCAGAA-AGAGGATA_S30.annotatednarrowpeaks ZF-25.out 11


counts: 
	python common_genes.py ZF-4.out ZF-7.out  
	python common_genes.py ZF-16.out ZF-19.out 
	python common_genes.py ZF-16.out ZF-22.out 
	python common_genes.py ZF-19.out ZF-22.out 
	python common_genes.py ZF-23.out ZF-29.out 
	python common_genes.py ZF-1.out ZF-6.out 



ZF-25.bed: 
	cut -f 1-6 8406-ZF-7_CAGAGAGG-GCGATCTA_S17.annotatednarrowpeaks > ZF-7.bed
	cut -f 1-6 8406-ZF-6_CTCTCTAC-GCGATCTA_S16.annotatednarrowpeaks > ZF-6.bed
	cut -f 1-6 8406-ZF-1_TAAGGCGA-GCGATCTA_S12.annotatednarrowpeaks > ZF-1.bed
	cut -f 1-6 8406-ZF-22_GTAGAGGA-ATAGAGAG_S27.annotatednarrowpeaks > ZF-22.bed
	cut -f 1-6 8406-ZF-19_GCTACGCT-ATAGAGAG_S26.annotatednarrowpeaks > ZF-19.bed
	cut -f 1-6 8406-ZF-29_CTCTCTAC-AGAGGATA_S33.annotatednarrowpeaks > ZF-29.bed
	cut -f 1-6 8406-ZF-4_TCCTGAGC-GCGATCTA_S14.annotatednarrowpeaks > ZF-4.bed
	cut -f 1-6 8406-ZF-23_TAAGGCGA-AGAGGATA_S28.annotatednarrowpeaks > ZF-23.bed
	cut -f 1-6 8406-ZF-16_TAGGCATG-ATAGAGAG_S23.annotatednarrowpeaks > ZF-16.bed
	cut -f 1-6 8406-ZF-14_TCCTGAGC-ATAGAGAG_S21.annotatednarrowpeaks > ZF-14.bed
	cut -f 1-6 8406-ZF-17_CTCTCTAC-ATAGAGAG_S24.annotatednarrowpeaks > ZF-17.bed
	cut -f 1-6 8406-ZF-25_AGGCAGAA-AGAGGATA_S30.annotatednarrowpeaks > ZF-25.bed


