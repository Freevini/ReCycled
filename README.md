# Circles
Michael Schmid & Vincent Somerville


Here, we are creating a pipeline that checks the circularity of genome assemblies and circluarises them according to the location of the dnaA.

In order to to this we are working on the following parts:

## 1. check query genome assembly

- check contigs: how many? How large?
- try to identify contig origines (bacteria, plasmids, phages)
- identify  if only 1 bacterial contigs

  ==> Take home: contig originates from?

## 2. check circularity

- map ends (10kb) onto eachother (minimap2) --> overlap or not
- merge contig ends and map raw reads+ look for:
  - even coverage (output: X coverage--> is this in the range)
  - spanning reads (number of reads)
    - PE (are pairs mapping)
    - ONT/Pacbio (is read spanning)

    ==> Take home: contig is circler or not!

## 3. circleries

  - take dnaA diamond mapping (+prodigal)
  - check if no genes are on reverse strand (prodigal)
  - check for presence of DNAa
  - check for DNAa orientation
  - cut and shift fasta start

## 4. Polish the newly stitched region (optional)

  - QC how many locations have change!
  - advise if additional polishing is required!

## 5 put Illumina contigs into synteny

Most genomes however are not finished (complete and circular), therefore it would be good to have an option for these too.

- check dnaA location of the REF sequence
- check dnaA lcoation of the query sequence
- minimap2 query sequence to REF sequence
- split first contigs at dnaA (+)
- arrange contigs according to minimap2 alignment
- (optional: change headers)

## 6. check contaminations

- PhiX
- humman contaminations
