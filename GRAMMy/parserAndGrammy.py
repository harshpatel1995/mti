# Fall 2016 - Taxonomic Inference
# Austin Vo, Felix Sosa, Harsh Patel & Trevor Ballard
# University of Central Florida
# Objective : Parse out read length and mismatches from a SAM File and run GRAMMy.

# Parse the MD Tag in the SAM file to extract and return the number of mismatches.
# Let's say the MD Tag is MD:Z:54A15 for a read length of 70. 
# The read matches the genome for first 54 bp followed by a mismatch at location of 'A' and 15 more matches. 
# Total matches = 54 + 15 = 69 and mismatches = 70 - 69 = 1

def parseMDTag(md_tag):

    # Extract the field value for the MD Tag i.e 54A15 for the example above
    md_field_value = md_tag.rsplit(':', 1)[1]
    # Count the number of upper case characters in the field  value i.e 1 for the example above 
    uppers = [l for l in md_field_value if l.isupper()]
    return len(uppers)  

# GRAMMy    
def grammy(R):

    sigma = 0.05
    # Number of sample reads
    N = len(R)
    #Number of references
    J = 1 
    
 '''
    # Setting up fake data.
    # This is a matrix of dictionaries.
    # Shows each the read length and number of mismatches between read and reference genome.
    
    R = [
        {'id': 1, 'read_length': 100, 'mismatches': {
            'g1': 3,
            'g2': 1
        }},
        {'id': 2, 'read_length': 200, 'mismatches': {
            'g1': 0,
            'g2': 5
        }},
        {'id': 3, 'read_length': 150, 'mismatches': {
            'g1': 8,
            'g2': 4
        }}
    ]
'''

    # Mixing coefficient.
    # "Size" of probability of reference genome in relation to the sample size.
    # All values of PI should sum up to 1.
    PI = {'g1': 1.0}
    
'''
    # Responsibility/Adjacency/Probability matrix.
    # Shows the probability that a reference genome is responsible for a read.
    # "Randomly" initialize each probability to have equal parts.
    
Z = {
	1: {
		'g1': 0.5, 
		'g2': 0.5
	}, 
	2: {
		'g1': 0.5,
		'g2': 0.5,
	}, 
	3: {
		'g1': 0.5,
		'g2': 0.5
	}
}
'''

    Z = {}
    genome1 = {'g1': 1.0}
    
    #Populate Z for all reads -> 200 for this example
    for read in R:
        print(read['id'])
        Z[read['id']] = genome1
        
    print(len(Z))

    # Iterating the EM step 10 times.
    # Will add threshold or log-likelihood checking later on.
    for i in range(10):
        
        # E-step: Get probability matrix.
        # This implements Algorithm (3) in the GRAMMy paper for the E-step.
        # However, it uses Algorithm (1) in the Sigma paper to calculate read probabilities instead of Algorithm (5) in GRAMMy paper.
        for read in R:
            for genome_id, num_mismatches in read['mismatches'].items():
                this_read_genome_prob = PI[genome_id] * ((sigma**num_mismatches) * ((1-sigma)**(read['read_length'] - num_mismatches)))
                all_read_genome_probs = 0
                for g_id, n_mis in read['mismatches'].items():
                    all_read_genome_probs += PI[g_id] * ((sigma**n_mis) * ((1-sigma)**(read['read_length'] - n_mis)))
                Z[read['id']][genome_id] = this_read_genome_prob / all_read_genome_probs
                
        # M-step: Update reference probability sizes.
        # This calculates Algorithm (4) in the GRAMMy paper.
        # Gets mixing coefficient for next iteration of EM.
        for genome_id, size in PI.items():
            prob_sum = 0
            for read_id, read_probs in Z.items():
                prob_sum += read_probs[genome_id]
            PI[genome_id] = prob_sum / N
        
        print('Genome prob sizes:')
        print(PI)
        print()
        print('Read probabilities')
        print(Z)
        print()
        print()
        
        # Add step for thresholding or log-likelihood convergence here.
    return

# Main performs parsing and calls GRAMMy. 
def main():

    # Open and read input SAM file
    file = open ("output.sam","r")
    lines = file.readlines()
    file.close()  
    # List of dictionary each storing information about a specific read
    read_matrix = []
    read_count = 1;
    for line in lines:
        text = line.split()
        # Exclude the header lines by using the fact that the header lines won't have more than 10 tags
        if len(text) > 10:
            # Store the read name
            read_id = text[0]
            # Store the reference that the read was mapped to
            reference_name = text[2]
            # Count the read length by parsing the read (not directly provided in SAM)
            read_length = len(text[9])
            # Count mismatches using the MD Tag
            mismatches = parseMDTag(text[12])
            
            # Add the information about the read in form of a dictionary to the list
            mismatch_dict = {'g1': mismatches}
            read_info = {'id': read_count,'read_length': read_length, 'mismatches': mismatch_dict}
            read_count = read_count + 1
            read_matrix.append(read_info)           
    
    #Call GRAMMy after populating all the reads
    grammy(read_matrix)
    
main()



