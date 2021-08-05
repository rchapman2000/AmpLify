#!/usr/bin/env python
import argparse
import os
import re
import shutil
import subprocess
import sys
import vcf
from tqdm import tqdm
from Bio import SeqIO
from utilities import runCommand
from utilities import getSamples
from utilities import parseCommonOptions

consoleSize = shutil.get_terminal_size().columns


# createVarStats - generates a file with statistics regarding variant positions
# to help provide context for why a variant was considered. The stats file contains
# the columns [Chromosome, position, reference_base, alternate_base, 
# pre_filtering_ref_base_count, pre_filtering_alternate_base_count, percent_alt_base,
# post_filtering_ref_base_count, post_filtering_alternate_base_count, percent_alt_base_post_filtering.
# 
# Parameters:
# sample - name of sample being run
# ref - path to the reference fasta file
# minCov - the minimum coverage threshold to filter variants
# minQual - the minimum base quality used in filtering
# minPerc - the minimum alternate base percentage to be 
#           accepted as a variant
def createVarStats(threads, inDir, sample, ref, minCov, minQual, minPerc):
    
    # Creats the file to contain the variable statistics
    filename = sample + '-varStats.tsv'
    statFile = open(filename, 'w+')
    
    # Writes two header lines to the file
    statFile.write('Variant Stats for Sample: {0} Coverage Cutoff: {1} Base Quality Cutoff: {2} Percent Cutoff: {3}\n'.format( \
        sample, minCov, minQual, minPerc))
    statFile.write('Chromosome\tPosition\tRef-Base\tAlt-base\tRef-base-count-no-quality-filtering\tAlt-base-count-no-quality-filtering\tAlt-base-percent-no-quality-filtering\tRef-base-count-post-filtering\tAlt-base-count-post-filtering\tAlt-base-percent-post-filtering\n')
    
    # Creats a pileup with no base-quality filtering
    process = runCommand("bcftools mpileup --threads {0} -Q 0 --count-orphans --no-BAQ --annotate 'Format/AD,FORMAT/DP,FORMAT/SP,INFO/AD' --max-depth 100000 -O b -f {1} {2} > {3}".format( \
        threads, ref, inDir + sample + '.bam', sample + '-no-base-qual-pileup.bcf'))
    # Creats a VCF from the unfiltered pileup
    process = runCommand('bcftools call --threads {0} -A -m -v -Ov -o {1} {2}'.format(threads, sample + '-no-base-qual.vcf', sample + '-no-base-qual-pileup.bcf'))
    
    # Opens the unfiltered and filtered (previously created) VCF files
    unFiltered = vcf.Reader(open(sample+'-no-base-qual.vcf', 'r'))
   
    # Defines what the minimum percent alternate abundance to show will be.
    # Defaults to 30, unless the user specified a lower minimum percentage
    percShow = 30.00 / 100.00
    if minPerc < percShow:
        percShow = minPerc / 100.00

    # Loops over every record in the unfiltered vcf file to create
    # a stats file.
    for record in unFiltered:
        # Grabs the chromosome, position, reference base,
        # alternate base, unfiltered reference count,
        #  and unfiltered alternate count.
        Chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = record.ALT
        unFiltRefCount = record.INFO['AD'][0]
        unFiltAltCount = record.INFO['AD'][1]
        unFiltAltPerc =  unFiltAltCount / (sum(record.INFO['AD']))
        # If the unfiltered alternate base percentage is below the
        # minimum percentage to be shown in the file, this record is skipped.
        if unFiltAltPerc >= percShow:
            # Opens the vcf file created using a pileup with basequality filtering
            # and loops over the records in the file.
            filtered = vcf.Reader(open(sample+'.vcf', 'r'))
            for record2 in filtered:
                # Checks whether the chromosome and position of the filtered
                # record match those from the unfiltered record.
                if Chrom == record2.CHROM and pos == record2.POS:
                    filtRefCount = record2.INFO['AD'][0]
                    filtAltCount = record2.INFO['AD'][1]
                    filtAltPerc = filtAltCount / (filtAltCount + filtRefCount)
                    
                    # Writes to the stats file the record information
                    statFile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(Chrom, pos, ref, alt, unFiltRefCount, unFiltAltCount, unFiltAltPerc, filtRefCount, filtAltCount, filtAltPerc))

    statFile.close()     
   

# callVariants - uses lofreq to call variants based on the alignments
#
# Parameters:
# threads - the number of CPU threads available
# sample - the name of the sample
# ref - the path to the reference genome 
# minCov - the coverage cutoff for an alternative base to be considered a variant
# minPerc - the percent abundance cutoff for an alternative base to be considered a variant
# minQual - the base-quality cutoff to be considered in the variant calculations
#
# Returns:
# The number of Variants applied to the sample
def callVariants(threads, inDir, sample, ref, minCov, minQual, minPerc):

    # Indexes the sorted bam file
    process = runCommand('samtools index -@ {0} {1}'.format(threads, inDir + sample + '.bam'))
   
    # Creates a pileup based on the clipped and sorted bam file. 
    # The pileup uses the minQuality to remove low quality bases, additionally the
    # --count-orphas and --no-BAW options are provided to ensure that no base score
    # recalibration occurs and no unstranded alignments are excluded. The --max-depth
    # parameter ensures that all reads are used at each site (bcftools defaultly
    # considers 250 reads per position max). Finally, the AD option in included in the
    # annotate to grab allele (base) counts after filtering (DP does not appear to 
    # take into account bases filtered by the quality filter).
    process = runCommand("bcftools mpileup --threads {0} -Q {1} --count-orphans --no-BAQ --annotate 'Format/AD,FORMAT/DP,FORMAT/SP,INFO/AD' \
                         --max-depth 100000 -O b -f {2} {3} > {4}" \
                         .format(threads, minQual, ref, inDir + sample + '.bam', sample + '-pileup.bcf'))

    # The pileup is converted to a vcf file. The -v option
    # ensures that only variant sites are called. The -A
    # option ensures that all alternative alleles are considered.
    # And the -m option initiates multivariant calling
    process = runCommand('bcftools call --threads {0} -A -m -v -Ov -o {1} {2}'.format(threads, sample + '.vcf', sample + '-pileup.bcf'))

    # The variants in the called vcf are filtered based on depth of coverage (AD[0] + AD[1])
    # as well as whether the alternative allele is present at a percentage
    # higher than the specified minimum percentage threshold.
    process = runCommand("bcftools filter --threads {0} -i '(INFO/AD[0]+INFO/AD[1])>={1} && (INFO/AD[1] / (INFO/AD[0] + INFO/AD[1])) > {2}' {3} > {4}" \
                         .format(threads, minCov, minPerc/100.00, sample + '.vcf', sample + '-filtered.vcf'))

    # Zips and indexes the filtered vcf file
    process = runCommand('bgzip {0}'.format(sample + '-filtered.vcf'))
    process = runCommand('tabix {0}'.format(sample + '-filtered.vcf.gz'))

    # Calls the create variant/filtering satistics method 
    createVarStats(threads, inDir, sample, ref, minCov, minQual, minPerc)    

    # Creates a stats file for the filtered vcf file which will be
    # used to calculate the number of variants in this sample.
    process = runCommand('bcftools stats {0} > {1}'.format(sample + '-filtered.vcf.gz', sample + '-vcf-stats.txt'))
    
    # Opens the vcf stats file, and looks for the line  
    # containing "numbe rof records". The number of variants 
    # is extracted from this line via regex
    numVariants = ''
    statsFile = open(sample + '-vcf-stats.txt')
    line = statsFile.readline()
    while line:
        pattern = 'number of records:\t(\d+)'
        m = re.search(pattern, line)
        if m:
            numVariants = m.group(1)
        line = statsFile.readline()
    return numVariants

# findMaskCoordinates - Creates a BED file that will contain positions that
# are to be masked in the final consensus file. The method starts by creating
# a BED file of every coordinate in the reference genome and then creates
# a BED file for every coordinate that we wish to keep in the final consensus
# based on the depth of coverage within the pileup file. Finally, the bedtools
# intersect command is used to identify coordinates that are found in the reference
# but not in the sequencing data coordinates. 
#
# Parameters:
# threads - the number of CPU threads available
# sample - the name of the sample being analyzed
# ref - the path to the reference fasta file
# refName - the name of the reference
# minCov - the minimum coverage threshold to be included in consensus
#
# Returns:
# the name of a bed file containing coordinates to be masked.
def findMaskCoordinates(threads, sample, ref, refName, minCov):
    
    # Creates a reference bed file to write to and grabs the reference sequences 
    # from the reference fasta
    refPositions = open(refName + '-positions.bed', 'w+')
    refSeqs = list(SeqIO.parse(ref, "fasta"))

    # Loops over every sequenc in the reference fasta (incase there are multiple chromosomes
    for refSeq in refSeqs:
        # For each position in the sequence, it writes a bed entry with the
        # id of the sequence (either the reference name or the chromosome)
        # followed by the position.
        for i in range(0, len(refSeq)):
            refPositions.write("{0}\t{1}\t{2}".format(refSeq.id, i, i + 1))
            if i != len(refSeq) - 1:
                refPositions.write("\n")

    refPositions.close()
    
    # Creates a VCF format file containing every point based on the pileup created earlier.
    process = runCommand('bcftools call --threads {0} -A -m -Ov -o {1} {2}'.format(threads, sample + '-all.vcf', sample + '-pileup.bcf'))
    
    
    # Parses the VCF and creates a BED file for every position with a coverage (based on the sum of the AD)
    # that fits the minimum coverage threshold. This bed will be interested with the reference positions,
    # and any position that is not found in the intersection is either not found in the alignment or
    # falls below the coverage threshold. Thus, these sites should be masked.
    reader = vcf.Reader(open(sample + '-all.vcf', 'r'))
    bedOut = open(sample + "-noMask.bed", 'w+')

    # Loops over every record in the vcf file and grabs only those that 
    # fit the minimum coverage threshold
    for record in reader:
        if (sum(record.INFO['AD']) >= minCov):
            # Bed is in format CHR, START_POS, END_POS
            bedOut.write("{0}\t{1}\t{2}\n".format(record.CHROM, record.POS - 1, record.POS)) 
    bedOut.close()

    # Uses bedtools intersect to find the positions present in the reference that are not found in the 
    # file containing positions that meet minimum coverage threshold.
    process = runCommand('bedtools intersect -a {0} -b {1} -v > {2}'.format(refName + '-positions.bed', sample + '-noMask.bed', sample + '-mask.bed'))

    return (sample + "-mask.bed")

# makeConsensus - Uses the variants to generate a consensus genome 
# and then masks areas of low coverage in the genome. If a base is thrown 
# out due to poor quality, it is still considered in the calculated using 
# some tools. But, the AD field in a VCF contains exact counts of bases 
# (after filtering). Thus, this can be use don every position in the file 
# to get the coverage after filtering.
#
# Parameters:
# threads - the number of CPU threads available
# sample - the name of the sample
# ref - the path to the reference genome
# minCov - the minimum depth of coverage threshold for a position to be masked
# sMask - the position to mask bases before to account for regions not covered by amplicon sequencing
# eMask - the position to mask bases after to account for region not covered by amplicon sequencing
def makeConsensus(threads, sample, ref, refName, minCov):

    # Applies the variants to the reference
    process = runCommand('bcftools consensus -o {0} -f {1} {2}'.format(sample + '-unmasked.fasta', ref, sample + '-filtered.vcf.gz'))

    maskBed = findMaskCoordinates(threads, sample, ref, refName, minCov)

    # Masks the fasta using the bed file of low coverage spots created previously.
    process = runCommand('bedtools maskfasta -fi {0} -bed {1} -fo {2}'.format(sample + '-unmasked.fasta', maskBed, sample + '-masked-consensus.fasta'))

    # Unwraps the fasta to be all on one line.
    process = runCommand("awk '{{if(NR==1) print \">{0}\"; else printf(\"%s\", $0); next}}' {1} > {2}".format(sample, sample + '-masked-consensus.fasta', sample + '-final.fasta'))

# calcCovPerc - calculates the percent coverage of the final consensus
# fasta file based on the given reference.
# 
# Parameters:
# sample - the name of the sample
# refPath - the path to the reference file.
# 
# Returns:
# The percent coverage in the format '###.##%'
def calcCovPerc(sample, refPath):
    
    # The length of the refernce is calculated
    # by adding the length of each line (excluding
    # the header) together.
    ref = open(refPath)
    refLength = 0
    ref.readline()
    line = ref.readline()
    while line:
        refLength = refLength + len(line)
        line = ref.readline()

    # Opens the consensus fasta and skips the header
    fasta = open(sample + '-final.fasta')
    fasta.readline()
    
    # The consensus should be unwrapped (thus on one line),
    # so the sequence line is read, and the number of N's and '-'s is calculated 
    # using python's string methods.
    sequence = fasta.readline()
    ncount = sequence.count('N')
    dashCount = sequence.count('-')

    # The coverage is calculated by subtracting the number of Ns and '-'s from
    # the length of the consensus sequence and dividing this
    # by the length of the reference
    coverage = "{:.2f}".format(((len(sequence) - (ncount + dashCount)) / refLength) * 100) + "%"

    return coverage

def generateConsensus():

    parser = argparse.ArgumentParser(usage="ampLify generateConsensus -i DIR -r FILE [options]", \
                                     description = "Performs variant calling and generates a consensus \
                                                    genome based on alignments.")
    parser.add_argument('-i', '--input', required=True, type=str, help='[Required] - Path to directory with input BAM files', \
             action='store', dest='dir')
    parser.add_argument('-r', '--reference', type=str, required=True, help='[Required] - Path to a reference genome fasta file', action='store', dest='ref')
    parser.add_argument('-o', '--output', type=str, required=False, help='Name of the directory to send results [Default = \'./results\']', \
        action='store', dest='out')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads available', action='store', dest='threads')
    parser.add_argument('--min-depth', type=int, required=False, \
        help='Minimum Depth of Coverage for a variant/region to be included in final consensus [Default = 10]', \
        action='store', dest='minCov')
    parser.add_argument('--min-qual', type=int, required=False, \
        help='Minimum Base Quality required for base to be included in variant calling [Default = 20]', \
        action='store', dest='minQual')
    parser.add_argument('--percentCutoff', type=float, required=False, \
        help='Percent abundance threshold of an alternative base to be considered a variant [Default = 50.00]', \
        action='store', dest='minPerc')

    currentDir = os.getcwd() + '/'

    # Grab the arguments
    args = parser.parse_args()
    
    inDir = ''
    threads = 0
    refPath = ''
    refName = ''
    outName = ''

    inDir, threads, refPath, refName, outName = parseCommonOptions(args)

    print(" Running Generate Consensus Module ".center(consoleSize,'#'))
    # Sets the minimum coverage threshold at 10 (default) or a value
    # specified by the user
    minCov = 10
    if (args.minCov):
        if (args.minCov < 0):
            sys.exit("{0} is not a valid coverage threshold".format(args.minCov))
        else: 
            minCov = args.minCov

    # Sets the minimum base quality threshold at 20 (default) or a
    # value specified by the user
    minQual = 20
    if (args.minQual):
        if (args.minQual < 0):
            sys.exit("{0} is not a valid quality score".format(args.minQual))
        else: 
            minQual = args.minQual

    # Sets the minumum alt-base percent threshold at 50% (default) or a
    # value specified by the user.
    minPerc = 50.00
    if (args.minPerc):
        if (args.minPerc < 0):
            sys.exit("{0} is not a valid percentage".format(args.minPerc))
        else: 
            minPerc = args.minPerc

    resultsDir = currentDir + outName + '/'
    workingDir = resultsDir + 'Working' + '/'
    
    # Creates the results and working directory
    try:
        os.mkdir(resultsDir)
    except FileExistsError as err:
            #Directory already exists
            print("Directory {0} already exists, proceeding anyway...".format(resultsDir))
    try:    
        os.mkdir(workingDir)
    except FileExistsError as err:
            #Directory already exists
            print("Directory {0} already exists, proceeding anyway...".format(workingDir))

    samples = getSamples(inDir, ".+?(?=\.bam)")

    stats = []
    
    pbar = tqdm(samples)

    os.chdir(workingDir)
    process = runCommand('cp {0} .'.format(refPath))
    refPath = workingDir + os.path.basename(refPath)

    for sample in pbar:
        pbar.set_description("Sample: {0}".format(sample))
        samplewd = os.getcwd() + '/' + sample + "-Working/"
        try:
            os.mkdir(samplewd)
        except FileExistsError as err:
            #Directory already exists
            print("Directory {0} already exists, proceeding anyway...".format(sample + "Working/"))
        os.chdir(samplewd)
        
        try:
            # Creates the list to store the statistics calculated for this sample
            # which will be written to the stats file after the sample has completed
            sampleStats = []
            sampleStats.append(str(sample))

            # Runs the variant calling step, appends the variant count to the statistics list,
            # and moves the sample variant statistics file to the results directory
            pbar.write('Calling Variants - Minimum Coverage Cutoff: {1} Minimum Quality Cutoff: {2} Minimum Alt Percent Cutoff: {3}'.format(sample, minCov, minQual, minPerc))
            varCount = callVariants(threads, inDir, sample, refPath, minCov, minQual, minPerc)
            sampleStats.append(str(varCount))
            runCommand('mv {0} {1}'.format(sample + '-varStats.tsv', resultsDir))
        
            # Runs the consensus generation step
            pbar.write('Generating, Filtering, and Masking Consensus for sample {0}. Minimum Coverage Cutoff: {1}'.format(sample, minCov))
            makeConsensus(threads, sample, refPath, refName, minCov)
        
            # Calculates the coverage percent of the conesensus and appends the 
            # percentage to the statistics file
            covPerc = calcCovPerc(sample, refPath)
            sampleStats.append(str(covPerc))

            # Moves the consensus fasta to the results directory
            runCommand('mv {0} {1}'.format(sample + '-final.fasta', resultsDir))

            # Appends the sample statistics list to the complete list 
            stats.append(sampleStats)
        
            pbar.write('\n')
        except subprocess.CalledProcessError:
            print(" [!] Error encountered while analyzing sample {0} [!] ".format(sample).center(consoleSize,'-'))
            print("\n\n")
       
        os.chdir(workingDir)

    print(" Generate Consensus Module Completed - Results located at {0} ".format(resultsDir).center(consoleSize, "#"))

    # Creates that tab-delimited statistics file and writes each sample's statistics. 
    statsFile = open(resultsDir + "/sampleStats.csv", "w+")   
    statsFile.write("Sample,num_variants,coverage\n")
    for x in stats:
        statsFile.write(','.join(x) + '\n')
    statsFile.close()

if __name__ == "__main__":
    generateConsensus()
