
Morning lecture 9-12

  - Introduction to nanopore sequencing [20]
  - Introduction to fieldable experiments [20]
  - Ultra-long reads [20]
  - Introduction to long read metagenomics [20]
  - Break [20]
  - Interactive experiment:
      - Running the MinION [30]
      - Basic data handling [30]
 
Afternoon practical

  - MinION QC [10]
  - Taxonomic identification [30]

```kraken2 -db kraken2-microbial-fatfree --threads 12 ebame18.fastq --report ebame18.report```


  - Long read assembly [30]

```minimap2 -x ava-ont -t 12  Kefir_RBK.fastq Kefir_RBK.fastq | gzip -1 > Kefir_RBK.paf.gz```

```miniasm -f Kefir_RBK.fastq Kefir_RBK.paf.gz > Kefir_RBK.contigs.gfa```

```awk '/^S/{print ">"$2"\n"$3}' Kefir_RBK.contigs.gfa | fold > Kefir_RBK.contigs.fa```

  - Bandage plots

  - Identification of contigs 

  - Polishing [20]
  - Functional assignments
  	  - 16S with Anvi'o 


#Flye/bin/flye --genome-size 18000000 --out-dir Kefir_RBK-metaflye --threads 24 --nano-raw Kefir_RBK.fastq```


