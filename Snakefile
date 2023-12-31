with open(config['SAMPLES']) as fp:
    SAMPLES = fp.read().splitlines()

print(SAMPLES) 

rule all:
         input:
            expand("galore/{sample}_R1_001_val_1.fq.gz", sample = SAMPLES),
            expand("galore/{sample}_R2_001_val_2.fq.gz", sample = SAMPLES),
            expand("{sample}.sam", sample = SAMPLES),
            expand("{sample}.bam", sample = SAMPLES), 
            expand("{sample}.sorted.bam", sample =SAMPLES),
            expand("{sample}.sorted.rmDup.bam", sample =SAMPLES),
            expand("{sample}.rmDup.unique.bam", sample =SAMPLES), 
            expand("{sample}.unique.counts", sample = SAMPLES),
            expand("{sample}.bigwig", sample = SAMPLES),
            expand("{sample}_tagdir/tagCountDistribution.txt", sample =SAMPLES), 
            expand("{sample}_tagdir/peaks.txt", sample = SAMPLES), 
            expand("macs2/{sample}_summits.bed", sample = SAMPLES)
rule trim: 
       input: 
           r1 = "{sample}_R1_001.fastq.gz",
           r2 = "{sample}_R2_001.fastq.gz"
       output: 
           "galore/{sample}_R1_001_val_1.fq.gz",
           "galore/{sample}_R2_001_val_2.fq.gz"
       conda: 'env/env-trim.yaml'
       shell: 
           """
           mkdir -p galore 
           mkdir -p fastqc 
           trim_galore --gzip --retain_unpaired --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """ 
rule align:
              input:               
                   "galore/{sample}_R1_001_val_1.fq.gz",
                   "galore/{sample}_R2_001_val_2.fq.gz"
              params:
                   index=config['INDEX'],
                   mem = config['MEMORY'],
                   cores = config['CORES']
              conda: 'env/env-align.yaml'
              output:
                   "{sample}.sam",
                   "{sample}_hist.txt" 
              shell:
                   """
                   bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 700 -p {params.cores} -x {params.index} -1 {input[0]} -2 {input[1]} -S {output[0]}  &> {output[1]}
                   """
rule samTobam:
             input: 
                 "{sample}.sam"
             output: 
                 "{sample}.bam"
             conda: 'env/env-align.yaml'
             shell: 
                   """
                   samtools view -bS {input} > {output}
                   """

rule sort: 
             input: 
                  "{sample}.bam"
             output: 
                    "{sample}.sorted.bam"
             conda: 'env/env-picard.yaml'
             shell: 
                   """ 
                   picard SortSam I={input}  O={output} SORT_ORDER=coordinate 
                   """ 

rule remove_duplicates:
       input: 
        "{sample}.sorted.bam"
       output: 
         "{sample}.sorted.rmDup.bam",
         "{sample}.rmDup.txt"
       conda: 'env/env-picard.yaml' 
       shell: 
           """
            picard MarkDuplicates I={input} O={output[0]} REMOVE_DUPLICATES=true METRICS_FILE={output[1]} 
           """

rule index: 
      input: 
         "{sample}.sorted.rmDup.bam"
      output: 
         "{sample}.sorted.rmDup.bam.bai"
      conda: 'env/env-align.yaml'  
      shell: 
          """
          samtools index {input} 
          """ 
rule bamCoverage: 
       input: 
        "{sample}.sorted.rmDup.bam",
        "{sample}.sorted.rmDup.bam.bai" 
       output: 
        "{sample}.bigwig" 
       params: 
         genome_size = config['Genome_Size'], 
         binsize = config['BINSIZE'], 
         num_processors = config['Num_Processors'] 
       conda: 'env/env-coverage'
       shell: 
          """ 
          bamCoverage -b {input[0]} -p {params.num_processors}  --normalizeUsing RPGC --effectiveGenomeSize {params.genome_size} --binSize {params.binsize} -o {output} 
          """ 
rule tag_dir: 
      input: 
       "{sample}.sorted.rmDup.bam"
      output: 
         "{sample}_tagdir/tagCountDistribution.txt" 
      params: 
        "{sample}_tagdir" 
      conda:'env/env-peaks.yaml' 
      shell: 
        """ 
        makeTagDirectory {params} -single {input} 
        """ 
rule findPeaks: 
      input: 
        "{sample}_tagdir/tagCountDistribution.txt" 
      params: 
        "{sample}_tagdir"
      output: 
        "{sample}_tagdir/peaks.txt",
      conda: 'env/env-peaks.yaml'  
      shell: 
         """ 
         findPeaks {params} -o auto 
         """


rule macs_bed: 
      input: 
         "{sample}.sorted.rmDup.bam"
      params: 
         "{sample}", 
         genome_size = config['Genome_Size'] 
      output: 
          "macs2/{sample}_summits.bed" 
      conda: 'env/env-peaks.yaml' 
      shell: 
           """
           macs2 callpeak -t {input} 
 	-f BAM -g {params[1]} \
	-n {params[0]} \ 
        --nomodel 
	--outdir macs2

         """
######	For visualisation and debugging 
rule get_unique:  
     input: 
       "{sample}.sorted.rmDup.bam"
     output: 
       #"{sample}.rmDup.unique.bam"
        "{sample}.unique.counts"
     conda: 'env/env-align.yaml' 
     shell: 
       """
       #samtools view -c -b -q 10 {input} > {output}
       samtools view -c -q 10 {input} > {output}
       """  
