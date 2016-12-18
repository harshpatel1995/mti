# GRAMMy

sigma = 0.05
N = 3
J = 2
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
PI = {'g1': 0.5, 'g2': 0.5}
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

for i in range(10):
	# Get probability matrix
	for read in R:
		for genome_id, num_mismatches in read['mismatches'].items():
			this_read_genome_prob = PI[genome_id] * ((sigma**num_mismatches) * ((1-sigma)**(read['read_length'] - num_mismatches)))
			all_read_genome_probs = 0
			for g_id, n_mis in read['mismatches'].items():
				all_read_genome_probs += PI[g_id] * ((sigma**n_mis) * ((1-sigma)**(read['read_length'] - n_mis)))
			Z[read['id']][genome_id] = this_read_genome_prob / all_read_genome_probs
			
	# Update reference probability sizes
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