# **scAllele**
_______________________________________
[![](https://img.shields.io/badge/scAllele-v0.0.9.2-blue)](https://test.pypi.org/project/scAllele/)

[Github](https://github.com/gxiaolab/scAllele/)

## **About**


scAllele is a versatile tool to detect and analyze nucleotide variants in scRNA-seq [(Preprint)](https://www.biorxiv.org/content/10.1101/2022.03.29.486330v1). \
scAllele makes use local reassembly via de-Bruijn graph to identify sequence differences and infers nucleotide variants at the read level. \
The read level variant-call allows for the analysis of the role variants in the context of splicing. 
Using mutual information, scAllele identifis allelic linkage between nucleotide variants and splicing
isoforms.

## **Table of contents**
- [Outline](#Outline)
- [Download](#Download)
- [Usage](#Usage)
	- [Basic usage](#Basic-usage)
	- [Preprocessing](#Preprocessing)
	- [Stranded data](#Stranded-data)
	- [Local variant call](#Local-variant-call)
	- [Filtering variants](#Filtering-variants)
	- [Training a new classifier](#Training-a-new-classifier)
- [Output](#Output)
- [Debug](#Debug)
 
_______________________________________

## **Outline**

| ![alt text](img/screenshot_Fig1.png) |
|:--:|
| *a. Illustration of the main algorithm of scAllele for variant calling. The reads and the reference genomic sequence overlapping a read cluster (RC) are decomposed into k-mers and are reasembled into a de Bruijn graph. b. Variants (green box in a) identified from the graph are then scored using a generalized linear model (GLM). The GLM was trained with different features (green box) to assign a confidence score to the variants. c. To identify allele-specific splicing (i.e., variant linkage), scAllele performs a mutual information calculation between nucleotide variants (SNVs, microindels) and intronic parts (where the ‘alleles’ are the different overlapping introns), to calculate allelic linkage of splicing isoforms.* |  

_______________________________________

## **Download**

scAllele is available through test.PyPi. To download simply type:

```
$ pip install -i https://test.pypi.org/simple/ scAllele
```

The download was tested with PyPi version >= 20.0.1

If succesful, the program is ready to use. The intallation incorporates console scripts entrypoints to directly call scAllele:
```
$ scAllele

Usage: 
	scAllele -b <file.bam> -g <genome.fa> -o <output prefix>

A variant caller and variant analysis tool for scRNA-seq data.

Options:
  -h, --help            show this help message and exit
  -b INPUT_BAM, --input-bam=INPUT_BAM
                        [Required] Input bam file, (or comma-seprated list of
                        bam files) sorted and indexed
  -o OUTPUT_PREFIX, --output-vcf=OUTPUT_PREFIX
                        [Required] Prefix of the output files
  -g GENOME, --genome-file=GENOME
                        [Required] Reference genome file (fasta format)

  Filtering Options:
    --AB=MINRATIOVAR    Minimum allelic ratio for the variant allele. Default:
                        0.01
    --AC=MINCOUNTVAR    Minimum read depth supporting the variant allele.
                        Default: 2
    --DP=MINCOVERAGE    Minimum read depth at the position of the variant.
                        Default: 5
    --min-base_position=MINREADPOS
                        Minimum mean distance for the variant from the read
                        ends. Default = 7
    --min-base_quality=MINBASEQUAL
                        Minimum mean base quality score to consider SNPs.
                        Default = 20

  Run Mode Options:
    --run_mode=RUN_MODE
                        Select <Variant_Caller> <Full> or <Training> mode.
                        Default: Full
    --glm_clf_name=GLM_CLASSIFIER_NAME
                        Prefix of the GLM pickle objects with the GLM models

  Linkage Options:
    --link_min_count=LINK_MIN_COUNT
                        Minimum number of common reads for linkage analysis.
                        Default = 10
    --link_min_mi=LINK_MIN_MI
                        Minimum mutual information for linkage analysis.
                        Default = 0.52

  Advanced Options:
    -c SEARCH_REGION, --region=SEARCH_REGION
                        Limit search to this region (chrom:start-end or
                        chrom). Default: All
    -n NODES, --nodes=NODES
                        Number of threads for parallel execution. Default = 64
    --min-map_quality=MINMAPQ
                        Minimum mapping quality to keep a read. Default = 40
    --max_soft_clip=MAXSOFTCLIP
                        Maximum length of soft-clipping allow at each end of
                        the read. Default = 5
    --kmer-size=KMER    k-mer size for the de-Bruijn graph assembly. Default:
                        15
    --strandedness=STRANDEDNESS
                        Select from ['fr-firststrand', 'fr-secondstrand', 'fr-
                        unstrand']. Default: 'fr-unstrand'
    --maxSeqErrorRate=MAXSEQERRORRATE
                        Maximum estimate of sequencing error rate. Default:
                        0.01
    --Ploidy=PLOIDY     Maximum ploidy to be considered. Default: 2.
 ```
_______________________________________

## **Usage**

### *Basic usage*

The minimum requirements to run scAllele are:
1. A bam file (sorted and indexed) `samtools sort file.bam file.sorted ; samtools index file.sorted.bam` 
2. A reference genome fasta file (indexed) `samtools faidx genome.fa`
3. A prefix for the output files.

```
$ scAllele -b testdata/gm12878.chr1.bam -g testdata/hg38.chr1.fa -o path/to/output_prefix
```
### *Preprocessing* 

scAllele only requires a bam file and a reference genome fasta file, however, in order to get optimal results it is recommended to pre-process the data:

| ![alt text](img/Fig_S1.png) |
|:--:| 
| *Recommended pipeline* |

### *Stranded data*
If your scRNA-seq is strand-specific, then you can specify the strandedness of your data (default: non-strand specific). \
Strand-specific data helps resolve ambiguous alignments on overlapping genes. It also helps detect more accurate ASAS events. \
Most strand-specific libraries in RNA-Seq are `fr-firststrand` (second read pair is sense to the RNA). You can specify this in your command:
```
$ scAllele \
    -b testdata/gm12878.chr1.bam \
    -g testdata/hg38.chr1.fa \
    -o path/to/output_prefix \
    --strandedness='fr-firststrand'
```
Alternatively, you can use the option `--strandedness=fr-secondstrand` if the first read pair is sense to the RNA.  

### *Local variant call*
By default, scAllele searches for variants in all the regions of the transcriptome covered by reads. If you wish to search for variants in a custom genomic interval, you can do so with the `-c` option. 
```
## Only search chromosome 1
$ scAllele \
    -b testdata/gm12878.chr1.bam \
    -g testdata/hg38.chr1.fa \
    -o path/to/output_prefix \
    -c chr1

## Only search within these coordinates
$ scAllele \
    -b testdata/gm12878.chr1.bam \
    -g testdata/hg38.chr1.fa \
    -o path/to/output_prefix \
    -c chr1:154582111-154628004
```

scAllele will search for read clusters within these regions only. Bare in mind that it's possible to find no read clusters in the spcified region, and that, if a specified region does not contain the entirety of a gene, it may miss some ASAS events. 


### *Filtering variants* 
Although it is recommended to filter variants downstream of your analysis (via bcftools or others), it's possible to filter variants from the start. If you wish, for example, to only report variants with 3 reads supporting the alternative allele (AC) and 5 reads overall, then you can run the following command:
```
$ scAllele \
    -b testdata/gm12878.chr21.bam \
    -g testdata/hg38.chr1.fa \
    -o path/to/output_prefix \
    --AC=3 \
    --DP=5
```
The default is `AC=2 and DP=2`. 

### *Training a new classifier* 
scAllele offers the option to retrain the variant classifier. Sequencing data from different platforms or resulting from different library preparation protocols may have different error profiles. If you wish to retrain scAllele's classifier run it in training mode: 

```
$ scAllele \
    -b testdata/gm12878.chr1.bam \
    -g testdata/hg38.chr1.fa \
    -o path/to/new_clf \
    --run_mode='Training' 
```

This will return a feature file (`.feature_matrix.tab`) containing the variant calls and all the features used for the training of the classifier.  
Then, run scAllele's training function. The supervised classifier will require a set of ground-truth variants to fit the model.  

```
$ scAllele_train \
    -i path/to/new_clf.feature_matrix.tab \
    -v truth.vcf \
    -g testdata/hg38.chr1.fa
```

This will return 3 pickle objects:
- path/to/new_clf.feature_matrix.tab.DELETION.glm.pickle
- path/to/new_clf.feature_matrix.tab.INSERTION.glm.pickle
- path/to/new_clf.feature_matrix.tab.SNP.glm.pickle

Finally, to use these new classifiers to call variants run:

```
$ scAllele \
    -b testdata/gm12878.chr1.bam \
    -g testdata/hg38.chr1.fa \
    -o new_path/to/output_prefix \
    --glm_clf_name path/to/new_clf.feature_matrix.tab  
```

_____________________________________


## **Output**

scAllele generates 4 files as output:
```
path/to/output_prefix.vcf
path/to/output_prefix.mi_summary.tab
path/to/output_prefix.read_cluster_info.tab
path/to/output_prefix.intronic_parts.bed
```
The first file `.vcf` is a standard vcf file reporting all the nucleotide variants found. The description of the tags and values are specified in the header.\
The second file `.mi_summary.tab` reports all the linkage events between variants or between variant and intronic part (ASAS). The mututal information and number of common reads between pairs of variants are reported. **NOTE**: All the testable linkages are presented. It is recommended to filter linkage events based on the mutual information and number of common reads as explained in the main publication.\
The third file `.read_cluster_info.tab` reports all the read clusters identified in the file. \
The fourth file `.intronic_parts.bed` reports the intronic parts identified together with the introns that form them. 

## **Debug**

If the installed version of scAllele is not the latest one, please try:

```
pip install -i https://test.pypi.org/simple/ scAllele==0.0.9.1
```

If one or more dependencies fail to install, make sure you have the latest version of pip:

```
pip install --upgrade pip
```

If the error persists, download the `requires.txt` file from this repository and install the dependencies prior to scAllele installation:

```
pip install -r requires.txt
pip install -i https://test.pypi.org/simple/ scAllele==0.0.9.1
``` 
