
Morning lecture 9-12

  - Long reads for metagenomics
  - Introduction to nanopore sequencing [20]
     - Portable sequencing applications [20]
     - Nanopore metagenomics
        - Challenges
        - Zymo stuff 
     - Bit of stuff on fermentome
     - Running the MinION
     - Break [20]
     - Understanding nanopore data
     - Basic data handling
        - Basecalling
        - QC  
     - Perspective on ultra-long reads [20]

Afternoon practical

  - MinION QC [10]
  - Taxonomic identification [30]

```kraken2 -db kraken2-microbial-fatfree --threads 12 ebame18.fastq --report ebame18.report```

  Questions:
    
  - What are the dominant species?
  - How many species are present in the sample?


  - Long read assembly [30]

```minimap2 -x ava-ont -t 12  Kefir_RBK.fastq Kefir_RBK.fastq | gzip -1 > Kefir_RBK.paf.gz```

```miniasm -f Kefir_RBK.fastq Kefir_RBK.paf.gz > Kefir_RBK.contigs.gfa```

```awk '/^S/{print ">"$2"\n"$3}' Kefir_RBK.contigs.gfa | fold > Kefir_RBK.contigs.fa```

  - How big is the assembly? How many contigs?
  - What is the largest contig?


## Bandage plots

  - What is everything?
 
## Back to Kraken

  - How many taxa are present now at species level?

  
  - Identification of contigs 

```kraken --db minikraken_20171013_4GB -t 4 Kefir_RBK.contigs.fa  > Kefir_RBK.contigs.kraken.txt```

  - Polishing [20]
     - Racon
      
  - Functional assignments
  	  - 16S with Anvi'o 

  - Phylogenetics
     - Core genome analysis


#Flye/bin/flye --genome-size 18000000 --out-dir Kefir_RBK-metaflye --threads 24 --nano-raw Kefir_RBK.fastq```


