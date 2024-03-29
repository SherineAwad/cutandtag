with open(config['CONTROL']) as fp:
    CONTROL = fp.read().splitlines()

print(CONTROL) 

with open(config['TREAT']) as fp:
    TREAT = fp.read().splitlines()

print(TREAT)

rule all:
         input:
            #Run on Control
            expand("galore/{sample}_R1_001_val_1.fq.gz", sample = CONTROL),
            expand("galore/{sample}_R2_001_val_2.fq.gz", sample = CONTROL),
            expand("{sample}.sam", sample = CONTROL),
            expand("{sample}.bam", sample = CONTROL), 
            expand("{sample}.sorted.bam", sample =CONTROL),
            expand("{sample}.sorted.rmDup.bam", sample =CONTROL),
            expand("{sample}.rmDup.unique.bam", sample =CONTROL), 
            expand("{sample}.unique.counts", sample = CONTROL),
            expand("{sample}.bigwig", sample = CONTROL),
            #Run on Treat   
            expand("galore/{sample}_R1_001_val_1.fq.gz", sample = TREAT),
            expand("galore/{sample}_R2_001_val_2.fq.gz", sample = TREAT),
            expand("{sample}.sam", sample = TREAT),
            expand("{sample}.bam", sample = TREAT),
            expand("{sample}.sorted.bam", sample =TREAT),
            expand("{sample}.sorted.rmDup.bam", sample =TREAT),
            expand("{sample}.rmDup.unique.bam", sample =TREAT),
            expand("{sample}.unique.counts", sample = TREAT),
            expand("{sample}.bigwig", sample = TREAT),
            
            expand("macs/{sample}_summits.bed", sample = TREAT),
            expand("macs/{sample}_peaks.narrowPeak", sample = TREAT),
                       
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

rule macs_bed: 
      input: 
         "{sample}.sorted.rmDup.bam"
      params: 
         "{sample}", 
         genome_size = config['Genome_Size'] 
      output: 
          "macs/{sample}_summits.bed",
          "macs/{sample}_peaks.narrowPeak" 
      conda: 'env/env-peaks.yaml' 
      shell: 
           """
           bash macs2.sh 
           """


rule annotateNarrowPeaks: 
      input: 
         "macs/{sample}_peaks.narrowPeak" 
      params: 
           genome= config['GENOME'], 
           gtf = config['GTF']  
      output: 
         "{sample}.annotatednarrowpeaks", 
         "{sample}.annotatednarrowpeaks.stats"
      shell: 
          """
          annotatePeaks.pl {input} {params.genome} -gtf {params.gtf}   -annStats {output[1]}  > {output[0]}    
          """ 

rule annotateSummitPeaks: 
      input:
         "macs/{sample}_summits.bed"
      params:
           genome= config['GENOME'],
           gtf = config['GTF']
      output:
         "{sample}.annotatedSummitpeaks",
         "{sample}.annotatedSummitpeaks.stats"
      shell:
          """
          annotatePeaks.pl {input} {params.genome} -gtf {params.gtf}   -annStats {output[1]}  > {output[0]}
          """

