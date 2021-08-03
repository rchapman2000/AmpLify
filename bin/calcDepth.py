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

# calcDepth - calculates either the maximum depth of coverage
# or the average depth of coverage if the maximum is too larger 
# using the samtools depth command based on a reference bed file
#
# Parameters: 
# threads - the number of CPU threads available
# sample - the name of the sample
# refPath - the path to the reference file
# refName - the name of the reference file
#
# Returns:
# the depth of coverage calculated.
def maxDepth(sample):
	
    process = runCommand("awk 'BEGIN{{n=0}}{{if($3>n){{n=$3}}}}; END{{print n}}' {0}".format(sample+'-coverage.txt'))

    return process.stdout.decode('ascii').replace('\n', '')

def meanDepth(sample):

    process = runCommand("awk 'BEGIN{{n=0; lines=0}}{{n+=$3; lines+=1}}; END{{print n/lines}}' {0}".format(sample+'-coverage.txt'))

    return process.stdout.decode('ascii').replace('\n', '')

def calcDepth():

    parser = argparse.ArgumentParser(usage = "ampLify calcDepth -i DIR -r FILE [options]", \
                                     description = "Calculates the maximum and average depth of coverage across samples")
    parser.add_argument('-i', '--input', required=True, type=str, \
        help='[Required] - Path to directory with input BAM files ', \
             action='store', dest='dir')
    parser.add_argument('-r', '--reference', type=str, required=True, help='[Required] - Path to a reference genome fasta file', action='store', dest='ref')
    parser.add_argument('-o', '--output', type=str, required=False, help='Name of the directory to send results [Default = \'./results\']', \
        action='store', dest='out')
    #parser.add_argument('-t', '--threads', type=int, help='Number of threads available', action='store', dest='threads')

    currentDir = os.getcwd() + '/'

    # Grab the arguments
    args = parser.parse_args()
    
    inDir = ''
    threads = 0 # Included because of conserved method, but no program in this module offers multithreading
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

    samples = getSamples(inDir, ".+?(?=\.bam)")
    stats = []
    
    pbar = tqdm(samples)

    os.chdir(workingDir)
    process = runCommand('cp {0} .'.format(refPath))
    refPath = workingDir + os.path.basename(refPath)
    process = runCommand('samtools faidx {0}'.format(refPath))

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
            process = runCommand("awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS $2}}' {0}.fai > {1}.bed".format(refPath, "./" + refName))
	
            process = runCommand('samtools depth -a -b {0} {1} -o {2}'.format(refName + '.bed', inDir + '/' + sample +'.bam', sample + '-coverage.txt'))

            maxDep = maxDepth(sample)
            sampleStats.append(maxDep)
            avgDep = meanDepth(sample)
            sampleStats.append(avgDep)

            stats.append(sampleStats)

        except subprocess.CalledProcessError:
            print(" [!] Error encountered while analyzing sample {0} [!] ".format(sample).center(consoleSize,'-'))
            print("\n\n")

    print(" CalcDepth Complete - Results located at {0} ".format(resultsDir).center(consoleSize, "#"))

    statsFile = open(resultsDir + "/sampleDepths.csv", "w+")   
    statsFile.write("Sample,max_depth,mean_depth\n")
    for x in stats:
        statsFile.write(','.join(x) + '\n')
    statsFile.close()


if __name__ == "__main__":
    calcDepth()
