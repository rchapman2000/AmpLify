import os
import subprocess
import re

# runCommand - wrapper function for subprocess to simplify calls
#
# Parameters:
# command - string representation of the command to run
#
# Returns:
# the subprocess result objectu
def runCommand(command):
   
   process = ''
   try:
       process = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
   except subprocess.CalledProcessError as e:
       print("\n\n")
       print(command)
       print(e.stderr.decode('utf-8'))
       raise
   
   return process

def getSamples(dir, regex):
    # Grabs the list of files in the input directory.
    fileList = os.listdir(dir)
    samples = []
    # Uses a regular expression to grab the sample name from the files.
    for f in fileList:
        result = re.match(regex, f)
        if result:
            samples.append(result.group(0))
    # Duplicates will be present (from grabbing sample from _1 and _2),
    # so deduplicate the sample list)
    samples = list(dict.fromkeys(samples))

    if (len(samples) == 0):
        sys.exit("Directory entered contains no input files")
    
    return samples

def parseCommonOptions(args):
    
    # Grabs the input directory specified by the user
    # and adds the path to the current working directory before.
    inDir = ''
    if os.path.isabs(args.dir):
        inDir = args.dir
    else:
        inDir = os.getcwd() + '/' + args.dir

    if inDir[-1] != '/':
        inDir = inDir + '/'
    
    print("Directory = ", inDir)
    
    # If the threads argument is not provided, it is 
    # defaultly set at 1. If it is provided, then the provided
    # number overrides the default value
    threads = 1
    if (args.threads):
        if (args.threads < 0):
            sys.exit("{0} is not a valid number of threads".format(args.threads))
        else: 
            threads = args.threads

    # Grabs the path to the reference file
    refPath = ''
    if os.path.isabs(args.ref):
        refPath = args.ref
    else:
        refPath = os.getcwd() + "/" + args.ref

    # Removes the file extension and path from the path
    # to get the name of the reference.
    refName = ".".join(os.path.basename(refPath).split('.')[:-1])

    # The output folder will be named 'results' by default,
    # but if the user provided a string to the results argument,
    # then the default will be overridden.
    outName = 'results'
    if (args.out):
        outName = args.out

    return inDir, threads, refPath, refName, outName