# **scAllele**

<a href="https://test.pypi.org/project/scAllele/"><img alt="PyPI"></a>

## About


scAllele is a versatile tool to detect and analyze nucleotide variants in scRNA-seq. 
scAllele makes use local reassembly via de-Bruijn graph to identify sequence differences and infers nucleotide variants at the read level. \
The read level variant-call allows for the analysis of the role variants in the context of splicing. 
Using mutual information, scAllele identifis allelic linkage between nucleotide variants and splicing
isoforms.

## Table of contents

[Download](#Download)

[Usage](#Usage)

[Output](#Output)


## Download

scAllele is available through PyPi. To download simply type:

```
$ pip install -i https://test.pypi.org/simple/ scAllele
```
The download was tested with PiPy>=20.0.1

If succesful, the program is ready to use. The intallation incorporates console scripts entrypoints to directly call scAllele:
```
$ scAllele
```
```
[Wed, 17 Nov 2021, 16:59:19] 0.157   COMMAND = '/usr/bin/scAllele'
A variant caller and variant analysis for scRNA-seq data
Usage:
     scAllele -b <file.bam> -g <genome.fa> -o <output prefix>
```

You can also view all the options by typing:

```
$ scAllele --help
```

```
Usage: 
	scAllele -b <file.bam> -g <genome.fa> -o <output prefix>

A variant caller and variant analysis for scRNA-seq data

Options:
  -h, --help            show this help message and exit
  -b INPUT_BAM, --input-bam=INPUT_BAM
                        [Required] Input bam file, sorted and indexed
  -o OUTPUT_VCF, --output-vcf=OUTPUT_VCF
                        [Required] Prefix of the output files
  -g GENOME, --genome-file=GENOME
                        [Required] Reference Genome (Fasta format)

  Filtering Options:
    --AB=MINRATIOVAR    Minimum allelic ratio for the variant allele. Default:
                        0.01
    --AC=MINCOUNTVAR    Minimum read depth supporting the variant allele.
                        Default: 2
    --DP=MINCOVERAGE    Minimum read depth at the position of the variant.
                        Default: 10
    --min-base_position=MINREADPOS
                        Minimum mean distance for the variant from the read
                        ends. Default = 7
    --min-base_quality=MINQUAL
                        Minimum mean base quality score to consider SNPs.
                        Default = 20

  Run Mode Options:
    --run_mode=RUN_MODE
                        Select <Variant_Caller> or <Training> mode. Default:
                        Variant_Caller
    --glm_clf_name=GLM_CLASSIFIER_NAME
                        Prefix of the GLM pickle objects with the GLM models

  Advanced Options:
    -c SEARCH_REGION, --region=SEARCH_REGION
                        Limit search to this region (chrom:start-end or
                        chrom). Default: All
    -n NODES, --nodes=NODES
                        Number of threads for parallel execution. Default = 16
    --min-map_quality=MINMAPQ
                        Minimum mapping quality to keep a read. Default = 40
    --max_soft_clip=MAXSOFTCLIP
                        Maximum length of soft-clipping allow at each end of
                        the read. Default = 5
    --kmer-size=KMER    k-mer size for de-Bruijn graph assembly. Default: 15
    --strandedness=STRANDEDNESS
                        Select from ['fr-firststrand', 'fr-secondstrand', 'fr-
                        unstrand']. Default: 'fr-unstrand'
    --maxSeqErrorRate=MAXSEQERRORRATE
                        Maximum estimate of sequencing error rate. Default:
                        0.01
    --Ploidy=PLOIDY     Maximum ploidy to be considered. Default: 2.
 ```
 
## Usage 

### Basic usage

The minimum requirements to run scAllele are:
1. A bam file (sorted and indexed) `samtools sort file.bam file.sorted ; samtools index file.sorted.bam` 
2. A reference genome fasta file (indexed) `samtools faidx genome.fa`
3. A prefix for the output files.

```
$ scAllele -b testdata/gm12878.chr1.bam -g testdata/hg38.chr21.fa -o path/to/output_prefix
```

### Stranded data
If your scRNA-seq is strand-specific, then you can specify the strandedness of your data (default: non-strand specific). /
Strand-specific data helps resolve ambiguous alignments on overlapping genes. It also helps detect more accurate ASAS events. /
Most strand-specific libraries in RNA-Seq are `fr-firststrand` (second read pair is sense to the RNA). You can specify this in your command:
```
$ scAllele -b testdata/gm12878.chr1.bam -g testdata/hg38.chr21.fa -o path/to/output_prefix --strandedness='fr-firststrand'
```
Alternatively, you can use the option `--strandedness=fr-secondstrand` if the first read pair is sense to the RNA.  

### Local search
By default, scAllele searches for variants in all the regions of the transcriptome covered by reads. If you wish to search for variants in a custom genomic interval, you can do so with the `-c` option. 
```
## Only search chromosome 1
$ scAllele -b testdata/gm12878.chr1.bam -g testdata/hg38.chr21.fa -o path/to/output_prefix -c chr1

## Only search within these coordinates
$ scAllele -b testdata/gm12878.chr1.bam -g testdata/hg38.chr21.fa -o path/to/output_prefix -c chr1:154582111-154628004

## Only search the chromosomes, or coordinates specified in this file
$ scAllele -b testdata/gm12878.chr1.bam -g testdata/hg38.chr21.fa -o path/to/output_prefix -c my_regions_of_interest.txt
```

The latter command reads a file with one region per line e.g.
```
$ cat my_regions_of_interest.txt

chr1:154582111-154628004
chr2
chr3
chr21:4589110-4595910
```
scAllele will search for read clusters within these regions only. Bare in mind that it's possible to find no read clusters in the spcified region, and that, if a specified region does not contain the entirety of a gene, it may miss some ASAS events. 


### Filtering variants 
Although it is recommended to filter variants downstream of your analysis (via bcftools or others), it's possible to filter variants from the start. If you wish, for example, to only report variants with 3 reads supporting the alternative allele (AC) and 5 reads overall, then you can run the following command:
```
$ scAllele -b testdata/gm12878.chr1.bam -g testdata/hg38.chr21.fa -o path/to/output_prefix --AC=3 --DP=5
```
The default is `AC=2 and DP=2`. 

## Output
