with open(config['SAMPLES']) as fp:
    SAMPLES = fp.read().splitlines()

print(SAMPLES) 

rule all:
         input:
            expand("galore/{sample}.R1_001_val_1.fastq.gz", sample = SAMPLES),
            expand("galore/{sample}.R2_001_val_2.fastq.gz", sample = SAMPLES),
            expand("{sample}.bam", sample = SAMPLES),
            expand("{sample}.sorted.bam", sample =SAMPLES),
            expand("{sample}.dupMark.bam", sample =SAMPLES),
            expand("{sample}.sorted.rmDup.bam", sample =SAMPLES),
            expand("{sample}.rmDup.unique.bam", sample =SAMPLES)
rule trim: 
       input: 
           r1 = "{sample}_R1_001.fastq.gz",
           r2 = "{sample}_R2_001.fastq.gz"
       output: 
           "galore/{sample}.R1_001_val_1.fastq.gz",
           "galore/{sample}.R2_001_val_2.fastq.gz"
       conda: 'env/env-trim.yaml'
       shell: 
           """
           mkdir -p galore 
           mkdir -p fastqc 
           trim_galore --gzip --retain_unpaired --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """ 
rule align:
              input:
                   "galore/{sample}.R1_001_val_1.fastq.gz",
                   "galore/{sample}.R2_001_val_2.fastq.gz"
              params:
                   genome=config['GENOME'],
                   mem = config['MEMORY'],
                   cores = config['CORES']
              output:
                   "{sample}.bam",
                   "{sample}_hist.txt" 
              conda: 'env/env-align.yaml'
              shell:
                   """
                   bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p {params.cores} -x {params.genome} -1 {input[0]} -2 {input[1]} -b {output[0]}  &> {output[1]}
                   """

rule sort: 
             input: 
                  "{sample}.bam"
             output: 
                    "{sample}.sorted.sam"
             shell: 
                   """ 
                   picard SortSam I={input}  O={output} SORT_ORDER=coordinate 
                   """ 

rule mark_duplicates: 
             input: 
               "{sample}.sorted.bam"
             output: 
               "{sample}.dupMark.bam",
               "{sample}.dupMark.txt"
             shell: 
                 """ 
                  picard MarkDuplicates I= {input} O ={output[0]} METRICS_FILE={output[1]}
                 """ 

rule remove_duplicates:
       input: 
        "{sample}.sorted.bam"
       output: 
         "{sample}.sorted.rmDup.bam",
         "{sample}.rmDup.txt" 
       shell: 
           """
            picard MarkDuplicates I={input} O={output[0]} REMOVE_DUPLICATES=true METRICS_FILE={output[1]} 
           """


rule get_unique:  
     input: 
       "{sample}.sorted.rmDup.bam"
     output: 
       "{sample}.rmDup.unique.bam"
     shell: 
       """
       samtools view -b -q 10 {input} > {output}
       """  
