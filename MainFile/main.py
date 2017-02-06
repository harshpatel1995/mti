''' 
Felix Sosa
02/05/2017

MainFile

TODO:
1. NCBI querying
2. Refactor and clean

Possible Improvements:
1. All I/O filenames declared here by user
2. Object-orientation

'''

import subprocess as sp

'''
0. User -> FASTA ('')
1. FASTA -> BWA => SAM ('twoReferences.sam')
2. SAM -> parser.py => CSV ('parsed_SAM_v4.csv')
3. CSV -> grammy.cpp => GRA ('results.gra')
4. GRA -> visualization => User
'''

# Step 0: Ask User for FASTA input and visualization arguments
visList = ["tree"]

# Step 1: Execute BWA

# Step 2: Excecute parser
sp.call(["python3", "parser.py"])

# Step 3.a: Compile GRAMMy
sp.call(["g++", "-std=c++11", "grammy.cpp", "-o", "grammy"])

# Step3.b: Execute GRAMMy
sp.call("./grammy")

# Step 4: Execute visualization
for arg in visList:
	sp.call(["python3", "visualize.py", arg, "results.gra"])