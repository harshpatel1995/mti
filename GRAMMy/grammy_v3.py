# GRAMMy

sigma = 0.05
N = 3
J = 2

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

# Mixing coefficient.
# "Size" of probability of reference genome in relation to the sample size.
# All values of PI should sum up to 1.
PI = {'g1': 0.5, 'g2': 0.5}

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