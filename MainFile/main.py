#!/usr/bin/env python3
''' 
MainFile

Usage:
    main.py --reference=grefs (<sample>|(--paired=<sampleA,sampleB>))...

TODO:
1. NCBI querying
2. Refactor and clean

Possible Improvements:
1. All I/O filenames declared here by user
2. Object-orientation

'''

import subprocess as sp
from docopt import docopt
import os.path
import pexpect

'''
0. User -> FASTA ('')
1. FASTA -> BWA => SAM ('twoReferences.sam')
2. SAM -> parser.py => CSV ('parsed_SAM_v4.csv')
3. CSV -> grammy.cpp => GRA ('results.gra')
4. GRA -> visualization => User
''' 

# Step 0: Get the user's input
options = docopt(__doc__, version='mti-vis 1.0')


# Split paired-end strings into two filenames
single = options['<sample>']
paired = [tuple(s.split(',')) for s in options['--paired']]
if any(not len(e) == 2 for e in paired):
    print('Error: paired samples must have exactly two files')
    raise SystemExit(0)

# Ensure everything ends with '.fasta', 'fna', .fastq'
if any((not s.endswith('.fasta')       # ends with .fasta
        and (not s.endswith('.fastq')) # ends with .fastq
        and (not s.endswith('.fna')))  # ends with .fna
        for s in single + [x for p in paired for x in p]):
    print('Error: Please only valid provide .fasta or .fastq files')
    raise SystemExit(0)
'''
# Check if .fasta/fastq files exist.
if any(not os.path.isfile(s)
        for s in single + [x for p in paired for x in p]):
    print('Error: The provided .fasta/fastq files do not exist in /bwa-0.7.15 and /MainFile')
    raise SystemExit(0)
'''
#If the user provided a reference, index it. Otherwise use the pre-indexed complete bacteria genome set.
grefs = options['--reference']
if(grefs):
    sp.call(["module load bwa && bwa index " + grefs], shell=True)
else:
    print("Error: Must provide a valid reference file.")
    raise SystemExit(0)
   
# Step 1 (single-end): Execute BWA
for read in single:
    nameBase = read[:-6]
    filenameBase = '../GRAMMy/' + nameBase
    filename = filenameBase + ".sam"  
    with open(filename, 'w') as f:
         sp.call(
	     ["/usr/bin/modulecmd bash load bwa && bwa mem " + grefs + " " + read], shell=True, stdout=f)
    # Step 2: Execute parser
    sp.call(["python3", "parser.py", nameBase + ".sam"], cwd="../GRAMMy")
	
    # Step 3: Execute GRAMMy
    sp.call(["./grammy", nameBase + '.csv'], cwd="../GRAMMy")


# Step 1 (paired-ends): Execute BWA
for reads in paired:
    nameBase = reads[0][:-6] + '_' + reads[1][:-6]
    filenameBase = '../GRAMMy/' + nameBase	
    filename = filenameBase + ".sam"
    with open(filename, 'w') as f:
         sp.call(
             ["module load bwa && bwa mem -v3 " + grefs + " " +  reads[0] + " " +  reads[1]], shell=True , stdout = f)

    # Step 2: Execute parser
    print(filename)
    sp.call(["python3", "parser.py", nameBase + ".sam"], cwd="../GRAMMy")
    
    # Step 3: Execute GRAMMy
    sp.call(["./grammy", nameBase + '.csv'], cwd="../GRAMMy")
    

#os.chdir("../GRAMMy/")
#gra_file_name = nameBase + '.gra'



