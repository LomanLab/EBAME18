# Afternoon practical
Orgnanise yourself into groups.

### MinION QC [10m]
As you saw in our morning demonstration, we handled the problem of basecalling the nanopore *squiggles* for you, using `guppy`: Oxford Nanopore's GPU-accelerated basecaller.
We provide the sequenced reads as `FASTQ`.
Don't forget, the `FASTQ` contains both sequence and associated per-base quality scores -- so your yield is about half of that of the file size. You may notice the reads are pretty long! They're also of variable length.
Before we do anything, we'll want to get an idea of how good our sequence data is (especially since it was prepared by two bioinformaticians).
For nanopore data, there's a handy package for generating statistics and plots, aptly named `NanoPlot`. Other packages are available.

If you've worked with Illumina short-read data, you may be familiar with `FastQC`, which doesn't work with long-read data and panics about the quality.

```
NanoPlot -p ebame18 -o ebame18_nanoplot/ -t 12 --fastq_rich ebame18.fastq
```

`--fastq_rich` instructs `NanoPlot` to generate summary information from the `FASTQ` itself.
This is possible as each read header is encoded with some metadata concerning the read.
Use `scp` to copy the new directory to your local machine to inspect the plots.

```
scp -r <username>@<host>:ebame_nanoplot/ .
```

#### Questions
  - Use `NanoPlot` to generate a report from the `FASTQ`
  - What is the distribution of read length?
  - What is the distribution of read quality?
  - What is the average error rate? You can convert `phred` scores to a probability with `10^(-Q/10)`

### Taxonomic identification of reads [30m]
As a means to broadly and quickly inspect the contents of the sample, we can assign taxons to the sequenced reads using `kraken2`. Given a library of sequenced genomes with known taxonomy (such as all of `RefSeq`), `kraken2` breaks up the sequences into *k-mers* (DNA substrings of length `k`), and keeps an index of all taxons associated with that k-mer.
By default `kraken2` will index all observed 35-mers.

We've already built a database beforehand, and you downloaded it just before our mid-morning break. Our database contains all of the `archaea`, `bacteria`, `fungi`, `protozoa`, `viral` sequences from NCBI RefSeq. As well as `UniVec_Core` - a set of sequences that are useful for screening for contaminants such as common sequencing library adaptors and cloning vectors.

It's worth noting that the database we've provided is reasonably large (30 GB), and `kraken2` will load this database into RAM to conduct taxonomic assignments.
If you're short of resources, you could split up your database and run `kraken2` multiple times; alternatively, you can reduce the number of k-mers in your index.
The latter option will drastically reduce the size of your database, representing only a small percentage of the original observed k-mers.
This will come at the expense of some accuracy, though the default 8 GB "`minikraken`" database can provide a good first-pass (it doesn't index any k-mers from eukaryotes, though).

```
kraken2 --db kraken2-microbial-fatfree --threads 12 ebame18.fastq --report ebame18.report.txt > ebame18.kraken.txt
```

The output that we've redirected to `ebame18.kraken.txt` will show the classification of each input sequence, as well as some metadata on the k-mers that were used to make that decision.
The `ebame18.report.txt` is a more general overview that aggregates the classifications.
We can use some more command line magic to sort and filter the `kraken2` report.

```
sort -k2 -nr ebame18.report.txt | less

# Sort by abundance and filter by genus, count with grep -c
sort -k2 -nr ebame18.report.txt | grep "\sG\s"

# Filter by species
sort -k2 -nr ebame18.report.txt | grep "\sS\s"
```

It's important to note that `kraken2` will do its best to assign taxonomy to your sequences, even from just a few k-mers.
We discussed the presence of *Campylobacter jejuni* (NCBI ID `197`) in our mystery sample, which might make one think twice about consuming it. Yet only one read was given this classification.
We can find more about that read with some `awk` and the non-report `kraken2` file:

```
awk '$3 == 197 {print $0}' ebame18.kraken.txt
# C    f0c49900-01ed-4912-98b7-3c8b50437f39    197    3158    0:918 197:5 0:2201
```
The final column in this output describes how k-mers across the read were assigned.
Here, we see that `kraken2` failed to find a taxon for the first 918 k-mers (`0:918`, recalling that ID `0` is used for the `unassigned` taxon), followed by just 5 k-mers for *C. jejuni* and 2201 more unassigned k-mers.
So, the read was assigned to *C. jejuni*, despite just 0.16% of the k-mers being assigned to it  -- or indeed at all!
I wrote a [small script](https://gist.github.com/SamStudio8/138fd1df9a215e87da3b917e6e564fe8) to drop any result for a read that had fewer than 25% of its k-mers assigned, and found we could reduce the number of unique taxa for the mystery sample from 398, to 53.
Be aware of this when making inference based on the `kraken2` report alone.

For some further reading, somebody recommended `krakenuniq` [https://github.com/fbreitwieser/krakenuniq](https://github.com/fbreitwieser/krakenuniq) to help avoid this pitfall.

#### Questions
  - How many species are present in the sample?
  - What are the dominant species?
  - What is the likely source of our mystery sample?


### Long read *de novo* assembly [30m]
Classifying the reads provides some insight to the underlying community, but for biological inference we need more detail, and so we will attempt to reconstruct the genomic structures.
We can reconstruct larger structures from the sequenced reads by leveraging overlaps in the DNA sequences to **assemble** reads into larger fragments. Through `minimap2` and `miniasm` this can be achieved in reasonable time.
`minimap2` and `miniasm` were written by Heng Li (who also started the `samtools` project) and are designed for fast mapping and assembly of noisy long reads, overcoming some of the assumptions made by tools designed for very-high quality short read sequences.

`minimap2` is an exhaustive alignment approach, performing pairwise alignment of all possible sequenced read pairs, as a means to enumerate all possible overlaps.
`minimap2` is easy to use and its parameters are set automatically based on the `-x` preset. Here we specify `ava-ont` (*all-versus-all*: maps each read against every other read, *Oxford Nanopore* reads).
We pipe the `stdout` through `gzip` (that `-1` is a one -- specifying a fast compression), as the `PAF` is a plaintext format and can become quite large.
Note the output file is `PAF` (**pairwise alignment format**), not `SAM`.

```
minimap2 -x ava-ont -t 12 Kefir_RBK.fastq Kefir_RBK.fastq | gzip -1 > Kefir_RBK.paf.gz
```

Given the reads and their overlaps, we can now use `miniasm` to perform assembly.
The tools are fast as they avoid error correction of the reads.

```
miniasm -f Kefir_RBK.fastq Kefir_RBK.paf.gz > Kefir_RBK.contigs.gfa
```

`miniasm` leverages the overlaps found by `minimap2` to assemble reads into **contigs** (*contigious DNA fragments*).
The output from `miniasm` is a `GFA` (**graphical fragment assembly** format).
That is, we not only have the assembled contigs, but information on how those contigs overlap with other contigs -- a graph.
It's of note that representing DNA sequences and their overlaps in the form of the a graph is a new area of research for assemblers, and for storing multiple reference genomes.
As a result, many tools are not yet compatible with `GFA` files, so we'll use a little command line magic to extract the contigs as a `FASTA`.

```
awk '/^S/{print ">"$2"\n"$3}' Kefir_RBK.contigs.gfa | fold > Kefir_RBK.contigs.fa
```
It's important to remember that this assembly represents a consensus, often with species specificity, but not necessarily.
Closely related sequences are often "squashed" and the true population variation (*i.e.* the *haplotypes*) is lost.

#### Questions
  - How big is the assembly? How many contigs?
  - What is the largest contig?


### Bandage plots [5m]

You probably want to see your graph.
Luckily for you, Ryan Wick has written a nifty tool called `Bandage` that will read a `GFA` and draw you a pretty picture.
Install `Bandage` ([https://rrwick.github.io/Bandage/](https://rrwick.github.io/Bandage/)) to your local machine and copy across your contigs `FASTA`.

#### Questions
  - What does the graph look like?
  - What structures do we have, what has been assembled?
  - How good do you think our consensus is?

### Polishing with `racon` [20m]
As mentioned, `miniasm` is fast as it avoids any form of error correction, so you'll expect the consensus sequence error to be as least as high as the underlying read error.
All is not lost, as we can improve the quality of our consensus assembly through **polishing**.
The intuition behind this is to take your assembled contigs, align the sequenced reads back to them and use the evidence from the reads to *correct* the contigs.

First we map the reads to our assembled consensus.
Note the `minimap2` preset is now `-x map-ont` (not `-x ava-ont`), as we want to map our reads to our assembled consensus.

```
minimap2 -t 12 -x map-ont Kefir_RBK.contigs.fa Kefir_RBK.fastq | gzip -1 > Kefir_RBK.reads-assembly.paf.gz
```

Now polish:
```
racon -t 12 Kefir_RBK.fastq Kefir_RBK.reads-assembly.paf.gz Kefir_RBK.contigs.fa > Kefir_RBK.contigs.racon.fa
```


### Taxonomic identification of contigs [15m]

Earlier, we assigned taxons to our reads.
We can assign taxons to our polished contigs in the same way.

```
kraken --db kraken2-microbial-fatfree --threads 12 Kefir_RBK.contigs.racon.fa > Kefir_RBK.contigs.kraken
```

#### Questions
  - How many taxa are present now at the species level?
  - Identify some contigs of interest using your `Bandage` plot, what taxon has been assigned to some contigs of interest?




### Functional assignment [20m] ?

  - Functional assignments
  	  - 16S with Anvi'o

  - Phylogenetics
     - Core genome analysis


#Flye/bin/flye --genome-size 18000000 --out-dir Kefir_RBK-metaflye --threads 24 --nano-raw Kefir_RBK.fastq```
