# Amplicon Sequencing Pipeline
The purpose of this tool is automate the process of assembling a consensus genome from amplicon sequencing reads while providing users control into how thier consensus is generated. Existing tools are built specifically for certain amplicon panels. While this tool has been built with SWIFT Biosciences products in mind, it could be adapted for other types of amplicon sequencing. Additionally, other tools have been "black-box" where coverage and base-quality cutoff have been determined by the developer. This tool allows users to specify cutoff parameters which can help tailor the results based on the quality of their sequences.

## Installation
To be Added

## Usage and Description
This tool consists of three modules which can be viewed by entering:
```
ampSeq -h
```
which displays
```
AmpSeq Pipeline
Usage: ampseq <module> [module options]

Modules:

  - preProcess
  - calcDepth
  - generateConsensus
```

---

### Pre-Process
This module takes in a directory of paried end fastq files, performs quality trimming using trimmomatic, generates quality metrics using FastQC, aligns the reads to a provided reference using bowtie2, and removes amplicon primers using the tool PrimerClip, developed by swift biosciences. 

**Trimming:**

This tool performs quality/adapter trimming using trimmomatic. Defaultly, it performs adapter trimming with a maximum of 2 seed mismatches, a palindrome clip threshold of 30, a simple clip threshold of 10, a minimum adapter length of 1, and retains the reverse read after adapter trimming. The tool also performs a sliding window trim, with a window size of 4 bases and a minimum average quality of 20. It performs moving leading and trailing clips from reads if qualities are less than 15. Finally, it keeps reads with a minimum length of 75 bp.
Optionally, users can specify a number of bases to hard-clip off of the head and tail ends of the reads, to improve read retention. We advising looking at the FastQC files to determine if the data should be rerun with these parameters.

See the [Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for more details

**Alignment:**

This tool performs local alignment to the provided reference using [Bowtie2](https://github.com/BenLangmead/bowtie2)

**PrimerClip:**

This tool utilizes the [PrimerClip](https://github.com/swiftbiosciences/primerclip) software developed by SwiftBioSciences for their products. It takes a an alignment file, and using a provided master file, trims the primer sequences. Masterfiles for swift products can be found on thier websites, however masterfiles could also be created for custom primer sets (but this has not been tested).

Master file format (Tab Delimited):
```js
sequence_name\ttarget_start\ttarget_end\tamplicon_name\tfprimer_start\tfprimer_end\tfprimer_name\trprimer_start\trprimer_end\trprimer_name\tforward_sequence
```

**Output:**
The aligned and clipped bam files will be found in the output directory along with a statistics file containing information about trimming and alignment. Additionally, a 'Working' diretory will be found in the output directory containing intermediate files for every sample. 

#### Usage:
```
ampSeq.sh preProcess -i DIR -r FILE -m FILE [options]
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
This module takes in a directory of clipped .bam alignment files and calculates the depth of coverage for each sample. We suggest running this to determine how much coverage a no template control might have, which should be used as a minimum coverage cutoff for the generate consensus module. It should be noted that this method is not fool proof, if the NTC has 1 region with high amplification, it will have a large maximum and average coverage. If these values seem too large, the .bam file can be visualized to determine this as well.

**Output:**
A CSV file containing the maximum and average coverage for each sample along with a 'Working' directory containing intermediate files. Files will be located in the output directory.

#### Usage:
```
ampSeq calcDepth -i DIR -r FILE [options]
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
This module takes in a directory of clipped .bam alignment files, performs variant calling, applies these variants to the reference, and masks the reference to reflect the areas that were sequenced with sufficient coverage.

**Variant Calling:**

This tool uses a [bcftools](http://samtools.github.io/bcftools/bcftools.html) and [pyVCF](https://pyvcf.readthedocs.io/en/latest/) to call/parse variants. First, a pileup is generated based on the bam file with all bases above the minimum base quality threshold. Then, the entries with variants are filtered to remove any entries where the depth is below the minimum coverage threshold and where the percent abundance of the alternative allele is below the minimum percentage cutoff. Finally, these variants are applied to the consensus using [bedtools](https://bedtools.readthedocs.io/en/latest/). 

**Masking:**

To determine what positions to mask, the tool utilizes the pileup generated during variant calling with positions above the minimum base quality threshold. Then, every position in the pileup that fits the depth of coverage threshold (after filtering low quality bases) is written to a bed file. This bed file is compared to a a bed file containing every nucleotide in the reference genome using [bedtools](https://bedtools.readthedocs.io/en/latest/). Positions that are found in the referene genome and not in the alignment positions are then placed into a new bed file to be masked. These positions are then applied to the conensus using [bedtools](https://bedtools.readthedocs.io/en/latest/).

**Output:**
The output directory will contain consensus fasta files for each sample as well as Variant calling statistics for each sample. The variant calling statistics will contain various metrics to view whether the variant calling is performing as you would like. It is a tsv with the format:
```
Chromosome	Position	Ref-Base	Alt-base	Ref-base-count-no-quality-filtering	Alt-base-count-no-quality-filtering	Alt-base-percent-no-quality-filtering	Ref-base-count-post-filtering	Alt-base-count-post-filtering	Alt-base-percent-post-filtering
```
Additionally, there will be a sample statistics file containing the number of variants applied for each sample and the percent coverage of the consensus comared to the reference. Finally, there is a 'Working' directory where intermediate files will be stored.

#### Usage:
```
ampSeq.sh generateConsensus -i DIR -r FILE [options]
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
| -t, --threads | Number of CPU threas to use |
| --min-depth | Minimum depth of coverag for a variant/region to be included in the final consensus [Default = 10] |
| --min-qual | Minimum base quality required for a base to be included in the coverage at that position [Default = 20 ] |
| --percentCutoff | Mimum Percent abundance of a variant for it to be applied as opposed to the referecne [Default = 50.00] |
| -h | Displays usage instructions |

