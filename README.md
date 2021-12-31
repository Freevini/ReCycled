# Circles
Michael Schmid & Vincent Somerville

*Sensu stricto* circular genomes do not have a start and end. However bacterial chromosomes have a single, unique origin of replication (oriC). DNA replication starts at the first position of the replication initiator protein, also called *dnaA* This results in a skewed read coverage distribution when sequencing growing bacteria. They generally show higher coverages at the origin of replication (i.e. *dnaA*) then at the terminus (Fig. 1).

<p align = "center">
<img src = "07_figures/PTR_figure.png" width="400">
</p>
<p align = "center">
Fig.1 - Peak-to-trough ratio (PTR) for a growing or non-growing bacteria. Taken from Korem et al. 2015 (DOI: 10.1126/science.aac4812)
</p>

In order to compare genomes with each other, it has been decided to start align assembled genomes at the beginning of the *dnaA* gene, which should be located on the positive strand (REF). Previous tools, e.g. circlator (REF) attempt to identify completely assembled genomes and start align them. However it is limited in it's utility of using different sequencing methods and the high throughput of many sequencing projects today. Very few genome assembly tools incorporate start alignment into the workflow. Most notably Unicycler (REF) and Tricylcler (REF) identifies circular, bacterial contigs and start aligns them accordingly. Nevertheless the large majority of genomes on NCBI remain not start aligned (Fig. 2). This hinders a streamlined comparative genomic approach, e.g. genome synteny cannot easily be inferred.


<p align = "center">
<img src = "07_figures/Plot_chromosome_startAligned_NCBI.png" width="400">
</p>
<p align = "center">
Fig.2 -Number of start aligned and non-start aligned genomes on NCBI assgined to chromosome level assemblies (done on 30.12.2021))
</p>


Here, we are creating a pipeline that checks the circularity of genome assemblies and circularises them according to the location of the dnaA.

Aims:

1. Create a tool that identifies and circularises complete bacterial contigs
2. Create a fast and scalable approach
3. use a few dependencies as possible

In order to to this we are working on the following parts:

## 1. Check and prepare query genome assembly

- check contigs: how many? How large?
- try to identify contig origines (bacteria, plasmids, phages)
- identify if only 1 bacterial contigs

  ==> Take home: contig originates from?

## 2. Check circularity

- map ends (1kb) onto eachother (minimap2) --> overlap or not
- merge contig ends and map raw reads to look for:
  - spanning reads (number of reads)
    - PE (are pairs mapping)
    - ONT/Pacbio (is read spanning)

    ==> Take home: contig is circler or not!

## 3. Circleries

  - take dnaA minimap2 mapping
  - check for presence of DNAa
  - check for DNAa orientation
  - cut and shift fasta start
  - advise if additional polishing is required!
  - merge all non-start aligned contigs and start aligned contig back together


## 4. put Illumina contigs into synteny (not implemented)

Most genomes however are not finished (complete and circular), therefore it would be good to have an option for these too.

- check dnaA location of the REF sequence
- check dnaA lcoation of the query sequence
- minimap2 query sequence to REF sequence
- split first contigs at dnaA (+)
- arrange contigs according to minimap2 alignment
- (optional: change headers)

## 5. check contaminations (not implemented)

- PhiX
- humman contaminations
