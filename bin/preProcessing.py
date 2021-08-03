import argparse
import os
import re
import subprocess
import shutil
import sys
from utilities import runCommand
from utilities import getSamples
from utilities import parseCommonOptions
from tqdm import tqdm

consoleSize = shutil.get_terminal_size().columns

# buildIndex - builds a bowtie2 index and a fasta index for 
# provided reference
#
# Parameters:
# threads - the number of CPU threads available
# path - the path to the reference file
# refName - the name of the reference file
#
# Returns:
# the directory where the bowtie2 index has been created
def buildIndex(threads, path, refName):
    
    print(" Reference {0} detected ".format(refName).center(consoleSize,'-'))
    print(" Building Reference Genome Index ".center(consoleSize, '-'))
    idxDir = os.getcwd() + "/" + refName + "/"
    
    # Attempts to make the directory to store the index.
    # If directory exists, the script continues.
    try:
        os.mkdir(idxDir)
    except FileExistsError as err:
        # Directory already exists
        print("Directory {0} already exists, proceeding...".format(idxDir))
    
    # Builds bowtie2 index for the provided reference.
    process = runCommand('bowtie2-build --threads {0} {1} {2}'.format(threads, path, idxDir + "idx"))
    
    # Builds a fasta index for hte provided references
    process = runCommand('samtools faidx {0}'.format(path))
    

    print(' Reference Building Complete '.center(consoleSize,'#'))
    print('\n')
    return idxDir

# trimReads - Uses trimmomatic to performing quality and adapter trimming of the reads
# 
# Parameters:
# threads - the number of CPU threads available
# sample - the name of the sample
# inDir - the path to the input fasta files
# adapterloc - the path to the file containing adapters to be
# trimmed 
# fcrop - number of bases to crop off of head of reads
# ecrop - number of bases to crop off of end of reads
#
# Returns:
# A list of values in the format:
#  [ number of forward reads before trimming,
#    number of reverse reads before trimming,
#    number of forward reads after trimming,
#    number of reverse reads after trimming ]
def trimReads(threads, sample, inDir, adapterloc, fcrop, ecrop):

    # Creates variables storing the paths to the forward
    # and reverse reads
    reads1 = inDir + "/" + sample + '_R1_*.fastq.gz'
    reads2 = inDir + "/" + sample + '_R2_*.fastq.gz'
    stats = []

    # Loops over the forward and reverse fastq files
    # and calculates the number of reads in either
    for f in [reads1, reads2]:
        process=subprocess.run('zcat {0} | wc -l'.format(f), stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
        stats.append(int(process.stdout.decode('ascii'))/4)
    
    os.mkdir(sample+"QC")
    
    # Runs fastqc to generate a QC file for the untrimmed reads
    process = runCommand('fastqc -t {0} {1} {2} -o {3}QC'.format(threads, reads1, reads2, sample))
    output = process.stdout

    # If the user specified a number of bases to be cropped off of the start of each read, the
    # HEADCROP option will be included in trimming
    headCrop = ''
    if (fcrop > 0):
        headCrop = "HEADCROP:{0}".format(fcrop)

    # If the user specified a number of bases to be cropped off of the end of each read, the 
    # CROP option will be included in trimming
    endCrop = ''
    if (ecrop > 0):
        endCrop = "CROP:{0}".format(ecrop)
        

    # Performs adapter and quality trimming on the reads
    process = runCommand('trimmomatic PE {0} {1} {2} {3} {4} {5} {6} {7} ILLUMINACLIP:{8}:2:30:10:1:true \
        LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:75'.format(reads1, reads2, sample + '_1.trimmed.fq', \
             sample + '_1.unpaired.fq', sample +"_2.trimmed.fq", sample + "_2.unpaired.fq", headCrop, endCrop, adapterloc))

    # Runs fastqc to generate a QC file for the trimmed reads
    process = runCommand('fastqc -t {0} {1} {2} -o {3}QC'.format(threads, sample + '_1.trimmed.fq', sample + '_2.trimmed.fq', sample))
    output=process.stdout

    # Loops over the trimmed forward and reverse fastq
    # files and calculates the number of reads in either
    for f in ["./" + sample + '_1.trimmed.fq', "./" + sample + '_2.trimmed.fq']:
        process = runCommand('cat {0} | wc -l'.format(f))
        stats.append(int(process.stdout.decode('ascii'))/4)
    return stats

# alignReads - Uses bowtie2 to align the reads to the reference genome
# 
# Parameters:
# threads - the number of CPU threads available
# sample - the name of the sample
# idxDir - the path to the bt2 index
# refName - the name of the reference
#
# Returns:
# The percentage of reads that aligned to the genome
def alignReads(threads, sample, idxDir, refName):

    # Uses the trimmed forward and reverse reads created in the either step
    reads1 = sample+ '_1.trimmed.fq'
    reads2 = sample + '_2.trimmed.fq'

    # Aligns the trimmed reads to the reference genome using bowtie2.
    # Alignment is local to ensure that no gaps are created in the reads
    process=runCommand('bowtie2 -p {0} -x {1} -1 {2} -2 {3} --local -S {4}'.format(threads, idxDir + "idx", reads1, reads2, sample + '-align.sam'))
 
    # Bowtie2 writes to the stderr. Grabs the alignment percentage 
    # from the bowtie2 output message 
    output = process.stderr.decode('ascii')
    outLines = output.split('\n')
    pattern = '(\d{1,3}\.?\d{0,2}%) overall alignment rate'
    alignRate = ''
    for line in outLines:
        m = re.match(pattern, line)
        if m:
            alignRate = m.group(1) 
    return alignRate

# clipPrimers - clips any primers present in the alignments, which will influence
# the overall consesnsus
#
# Parameters:
# threads - the number of CPU threads available
# sample - the name of the sample
# masterFile - the file containing a list of primers
# to be trimmed
def clipPrimers(threads, sample, masterFile):

    # Converts the alignment sam file into a bam file, sorts it
    # by read name, and outputs it to a sam file again.
    # This redundancy is due to the fact that primerclip requires
    # a sorted sam file, but only bam files can be sorted using samtools
    # sort.
    process = runCommand("samtools view -@ {0} -bS {1} | samtools sort -@ {0} -n -O sam > {2}".format( \
        threads, sample + '-align.sam', sample + '.sort.sam'))
    output=process.stdout

    # Uses swift's primerclip to clip any primers from the sorted
    # sam file.
    process = runCommand('primerclip {0} {1} {2}'.format(masterFile, sample + '.sort.sam', sample + '.clip.sam'))

    # Converts the clipped sam file into bam and then sorts
    process = runCommand('samtools view -@ {0} -bS {1} | samtools sort -@ {0} > {2}'.format(threads, sample + '.clip.sam', sample + '.bam'))

def PreProcessing():
    
    parser = argparse.ArgumentParser(usage="ampSeq.sh preProcess -i DIR -r FILE -m FILE [options]", \
                                     description = "Performs Quality Control and Adapter Trimming, aligns reads \
                                         to provided reference,and clips amplicon sequencing primers.")
    parser.add_argument('-i', '--input', required=True, type=str, \
        help='[Required] - Path to directory with input fastq files', \
             action='store', dest='dir')
    parser.add_argument('-r', '--reference', type=str, required=True, help='[Required] - Path to a reference genome fasta file', action='store', dest='ref')
    parser.add_argument('-m', '--master-file', required=True, type=str, help='[Required] Path to swift primer masterfile [Required for PreProcessing', \
        action='store', dest='master')
    parser.add_argument('-o', '--output', type=str, required=False, help='Name of the directory to send results [Default = results]', \
        action='store', dest='out')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads available', action='store', dest='threads')
    #parser.add_argument('--pre-process', help='Runs the PreProcessing Module', action='store_true', dest='preProc')

    parser.add_argument('--head-crop', type=int, required=False, \
        help='Number of bases to trim off the front of reads [Default = 0] (Suggested to run the pipeline or FASTQ prior to setting this option)', \
        action='store', dest='hCrop')
    parser.add_argument('--end-crop', type=int, required=False, \
        help='Number of bases to trim off the end of reads [Default = 0] (Suggested to run the pipeline or FASTQ prior to setting this option)', \
        action = 'store', dest = 'eCrop')

    currentDir = os.getcwd() + '/'

    # Grab the arguments
    args = parser.parse_args()
    
    inDir = ''
    threads = 0
    refPath = ''
    refName = ''
    outName = ''

    inDir, threads, refPath, refName, outName = parseCommonOptions(args)

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

    print(" Running PreProcessing Module ".center(consoleSize,'#'))
    masterFile = ''
    if (args.master):
        if os.path.isabs(args.master):
            masterFile = args.master
        else:
            masterFile = os.getcwd() + "/" + args.master
    else:
        sys.exit('You must specify a master file for preprocessing')

    # Sets the staring mask position to 0 (defaultly does not mask anything) 
    # or to the position supplied by the user.
    hCrop = 0
    if (args.hCrop):
        if (args.hCrop < 0):
            sys.exit("{0} is not a valid number of bases. Continuing with default of 0.".format(args.hCrop))
        else: 
            hCrop = args.hCrop
    # as -1 corresponds to the end of a string in python) or to the 
    # position supplied by the user.
    eCrop = 0
    if (args.eCrop):
        if (args.eCrop < 0):
            sys.exit("{0} is not a valid number of bases. Continuing with default of 0.".format(args.eCrop))
        else: 
            eCrop = args.eCrop

    samples = getSamples(inDir, ".+?(?=_R[12]_.*\.fastq.*)")
    stats = []

    os.chdir(workingDir)
    idxDir = buildIndex(threads, refPath, refName)
    process = runCommand('cp {0} .'.format(refPath))
    refPath = workingDir + os.path.basename(refPath)

    pbar = tqdm(samples)
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
        
            pbar.write("Running QC checks and performing quality/adaptor trimming for sample {0}".format(sample))
            # Runs the trimming step of the pipeline and appends the read counts to the statistics list
            readCounts = trimReads(threads, sample, inDir, os.path.dirname(__file__) + "/Swift_Adapters.fa", hCrop, eCrop)
            for x in readCounts:
                sampleStats.append(str(x))

            # Runs the alignment step and appends the alignment percentage to the sratistics list
            pbar.write("Aligning reads for {0} to reference {1}".format(sample, refName))
            alignPerc = alignReads(threads, sample, idxDir, refName)
            sampleStats.append(str(alignPerc))

            # Runs the primer clipping step  
            pbar.write('Clipping Primers for sample {0} based on primers from {1}'.format(sample, masterFile))
            clipPrimers(threads, sample, masterFile)

            # Moves clipped and sorted bam to the results directory
            runCommand('mv {0} {1}'.format(sample + '.bam', resultsDir))

            # Appends the sample statistics list to the complete list 
            stats.append(sampleStats)
        
            pbar.write('\n')
        except subprocess.CalledProcessError:
            print(" [!] Error encountered while analyzing sample {0} [!] ".format(sample).center(consoleSize,'-'))
            print("\n\n")
       
        os.chdir(workingDir)

    print(" PreProcessing Complete - Results located at {0} ".format(resultsDir).center(consoleSize, "#"))

    # Creates that tab-delimited statistics file and writes each sample's statistics. 
    statsFile = open(resultsDir + "/sampleStats.csv", "w+")   
    statsFile.write("Sample,pre_trim_forward_reads,pre_trim_reverse_reads,trimmed_forward_reads,trimmed_reverse_reads,percent_alignment\n")
    for x in stats:
        statsFile.write(','.join(x) + '\n')
    statsFile.close()

if __name__ == "__main__":
    PreProcessing()