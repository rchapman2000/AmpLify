# AmpLify - an automated amplicon sequencing analysis pipeline
The purpose of ampLify is to automate the process of assembling a consensus genome from amplicon sequencing reads while providing users control into how their consensus is generated. Existing tools are built specifically for certain amplicon panels. While ampLify has been built with SWIFT Biosciences products in mind, it could be adapted for other types of amplicon sequencing. Additionally, other tools have been "black-box" where coverage and base-quality cutoff have been determined by the developer. AmpLify allows users to specify cutoff parameters which can help tailor the results based on the quality of their sequences.

## Installation
AmpLify can be installed using conda:
```
conda create -n ampLify-env
conda install -c rchapman2000 amplify
```

## Usage and Description
ampLify consists of three modules which can be viewed by entering:
```
ampLify -h
```
which displays
```
AmpLify Pipeline
Usage: ampLify <module> [module options]

Modules:

  - preProcess
  - calcDepth
  - generateConsensus
```

---

### Pre-Process
This module takes in a directory of paired-end fastq files, performs quality trimming, generates quality metrics, aligns the reads to a provided reference, and removes amplicon primers using the tool PrimerClip, developed by swift biosciences.

**Trimming:**

AmpLify performs quality/adapter trimming using [trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf). Default, it performs adapter trimming with a maximum of 2 seed mismatches, a palindrome clip threshold of 30, a simple clip threshold of 10, a minimum adapter length of 1, and retains the reverse read after adapter trimming. The tool also performs a sliding window trim, with a window size of 4 bases and a minimum average quality of 20. It performs moving leading and trailing clips from reads if qualities are less than 15. Finally, it keeps reads with a minimum length of 75 bp.
Optionally, users can specify a number of bases to hard-clip off of the head and tail ends of the reads, to improve read retention. We advising looking at the FastQC files to determine if the data should be rerun with these parameters.

See the [Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for more details

**Alignment:**

AmpLify performs local alignment to the provided reference using [Bowtie2](https://github.com/BenLangmead/bowtie2). Bowtie2's inherent --very-fast-local alignment scoring is used to ensure that indels are still allowed despite local alignemtn being used. When we tested the pipeline using a gapped alignment strategy rather than local, chimeric and singleton reads were found to be aligned in several locations, skewing the consensus. Additionally, AmpLify does not allow for disconcordant or singleton alignments to be included to help reduce the presence of false positive alignments.

**PrimerClip:**

AmpLify utilizes the [PrimerClip](https://github.com/swiftbiosciences/primerclip) software developed by SwiftBioSciences for their products. It takes an alignment file, and using a provided master file, trims the primer sequences. Masterfiles for swift products can be found on their websites, however, masterfiles could also be created for custom primer sets (but this has not been tested).

Master file format (Tab Delimited):
```diff
sequence_name target_start  target_end  amplicon_name fprimer_start fprimer_end fprimer_name  rprimer_start rprimer_end rprimer_name  forward_sequence
```

**Output:**
The aligned and clipped bam files will be found in the output directory along with a statistics file containing information about trimming and alignment. Additionally, a 'Working' directory will be found in the output directory containing intermediate files for every sample.

#### Usage:
```
ampLify preProcess -i DIR -r FILE -m FILE [options]
```

#### Options:
**Required Arguments**
| Argument: | Description |
|-----------|-------------|
| -i, --input | Path to a directory containing fastq input files (in format *_R1.fastq and *_R2.fastq) |
| -r, --reference | Path to a reference genome in .fasta or .fa format|
| -m, --master-file | Path to a swift primer masterfile for primer clipping |

**Optional Arguments**
| Argument: | Description |
|-----------|-------------|
| -o, --output | Name and path of the directory to output files to [Default = './results'] |
| -t, --theads | Number of CPU threads to use |
| --head-crop | Number of bases to trim off of hte front of reads [Default = 0] |
| --end-crop | Number of cases to trim off of the end of reads [Default = 0] |
| -h | Displays usage instructions |

---

### Calc Depth
AmpLify takes in a directory of clipped .bam alignment files and calculates the depth of coverage for each sample. We suggest running this to determine how much coverage a no template control might have, which should be used as a minimum coverage cutoff for the general consensus module. It should be noted that this method is not foolproof, if the NTC has 1 region with high amplification, it will have a large maximum and average coverage. If these values seem too large, the .bam file can be visualized to determine this as well.

**Output:**
A CSV file containing the maximum and average coverage for each sample along with a 'Working' directory containing intermediate files. Files will be located in the output directory.

#### Usage:
```
ampLify calcDepth -i DIR -r FILE [options]
```

#### Options:
**Required Arguments**
| Argument: | Description |
|-----------|-------------|
| -i, --input | Path to a directory containing clipped bam input files |
| -r, --reference | Path to a reference genome in .fasta or .fa format|


**Optional Arguments**
| Argument: | Description |
|-----------|-------------|
| -o, --output | Name and path of the directory to output files to [Default = './results'] |
| -h | Displays usage instructions |

---

### Generate Consensus
AmpLify takes in a directory of clipped .bam alignment files, performs variant calling, applies these variants to the reference, and masks the reference to reflect the areas that were sequenced with sufficient coverage.

**Variant Calling and Masking:**

AmpLify uses a [bcftools](http://samtools.github.io/bcftools/bcftools.html) and [pyVCF](https://pyvcf.readthedocs.io/en/latest/) to call/parse variants. First, a pileup is generated based on the bam file with all bases above the minimum base quality threshold. Then, the entries with variants are filtered to remove any entries where the depth is below the minimum coverage threshold and where the percent abundance of the alternative allele is below the minimum percentage cutoff. 

Next, AmpLify separates out the SNPs and Deletions from the insertions. the SNPs and deletions are then applied to the reference maintaining a '-' for any deletion. Our reason for doing this is that when bcftools reports nucleotides, the position is not relative to the deletion. Thus, this would affect what positions are masked given our masking strategy.

To determine what positions to mask, AmpLify utilizes the pileup generated during variant calling with positions above the minimum base quality threshold. Then, every position in the pileup that fits the depth of coverage threshold (after filtering low-quality bases) is written to a bed file. This bed file is compared to a bed file containing every nucleotide in the reference genome using [bedtools](https://bedtools.readthedocs.io/en/latest/). Positions that are found in the reference genome and not in the alignment positions are then placed into a new bed file to be masked. These positions are then applied to the consensus using [bedtools](https://bedtools.readthedocs.io/en/latest/).

Finally, these insertions that pass the thresholds are applied to the consensus.

**Output:**
The output directory will contain consensus fasta files for each sample as well as Variant calling statistics for each sample. The variant calling statistics will contain various metrics to view whether the variant calling is performing as you would like. It is a tsv with the format:
```
Chromosome    Position    Ref-Base    Alt-base    Ref-base-count-no-quality-filtering    Alt-base-count-no-quality-filtering    Alt-base-percent-no-quality-filtering    Ref-base-count-post-filtering    Alt-base-count-post-filtering    Alt-base-percent-post-filtering
```
Additionally, there will be a sample statistics file containing the number of variants applied for each sample and the percent coverage of the consensus compared to the reference. Finally, there is a 'Working' directory where intermediate files will be stored.

#### Usage:
```
ampLify generateConsensus -i DIR -r FILE [options]
```

#### Options:
**Required Arguments**
| Argument: | Description |
|-----------|-------------|
| -i, --input | Path to a directory containing clipped bam input files |
| -r, --reference | Path to a reference genome in .fasta or .fa format|

**Optional Arguments**
| Argument: | Description |
|-----------|-------------|
| -o, --output | Name and path of the directory to output files to [Default = './results'] |
| -t, --threads | Number of CPU threads to use |
| --min-depth | Minimum depth of coverage for a variant/region to be included in the final consensus [Default = 10] |
| --min-qual | Minimum base quality required for a base to be included in the coverage at that position [Default = 20 ] |
| --percentCutoff | Minimum Percent abundance of a variant for it to be applied as opposed to the reference [Default = 50.00] |
| --indel-gap | Minimum distance between two indels required for them to both be included. Below this threshold, one of the indels will be filtered out [Default = 10] |
| -h | Displays usage instructions |

## Useful Tips
### Automating
The pre-process and generate consensus modules can be strung together in bash to create a singular pipeline:
```
ampLify preProcess -i FILE -r REFERENCE -m MASTERFILE -o OUTDIR/preProcess && ampLify generateConsensus -i 
```

## Citating AmpLify:
If you have used ampLify in your reserach, please include a citation in your work:
- Authors: Ryan C. Chapman, Michael R. Wiley
- Date: August, 2021
- Url: https://github.com/rchapman2000/AmpLify
- Institution: Univeristy of Nebraska Medical Center

Please also include citations for following dependencies of ampLify:

- bcftools & samtools - https://academic.oup.com/gigascience/article/10/2/giab008/6137722
- bedtools - https://academic.oup.com/bioinformatics/article/26/6/841/244688
- biopython - dx.doi.org/10.1093/bioinformatics/btp163
- bowtie2 - https://www.nature.com/articles/nmeth.1923
- fastqc - https://github.com/s-andrews/FastQC
- pyVCF - https://github.com/jamescasbon/PyVCF
- trimmomatic - https://academic.oup.com/bioinformatics/article/30/15/2114/2390096
- tqdm - https://zenodo.org/record/5149124#.YQwftDpOlpg

## Citations:

Andrews, S. (2010). FASTQC. https://github.com/s-andrews/FastQC

Bolger, A. M., Lohse, M., &amp; Usadel, B. (2014). Trimmomatic: A FLEXIBLE trimmer FOR Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170 

Casbon J. (2012). PyVCF. https://github.com/jamescasbon/PyVCF

Casper da Costa-Luis, Stephen Karl Larroque, Kyle Altendorf, Hadrien Mary, richardsheridan, Mikhail Korobov, Noam Yorav-Raphael, Ivan Ivanov, Marcel Bargull, Nishant Rodrigues, Guangshuo CHEN, Antony Lee, Charles Newey, James, Joshua Coales, Martin Zugnoni, Matthew D. Pagel, mjstevens777, Mikhail Dektyarev, … Nikolay Nechaev. (2021). tqdm: A fast, Extensible Progress Bar for Python and CLI (v4.62.0). Zenodo. https://doi.org/10.5281/zenodo.5149124

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., &amp; de Hoon, M. J. (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163 

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., &amp; Li, H. (2021). Twelve years of samtools and BCFtools. GigaScience, 10(2). https://doi.org/10.1093/gigascience/giab008 

Langmead, B., &amp; Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923 

Quinlan, A. R., &amp; Hall, I. M. (2010). BEDTools: A Flexible suite of utilities for comparing Genomic features. Bioinformatics, 26(6), 841–842. https://doi.org/10.1093/bioinformatics/btq033 
