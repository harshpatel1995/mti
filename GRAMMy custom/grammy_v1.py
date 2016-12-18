# GRAMMy

sigma = 0.05
num_reads = 3
num_genomes = 2
reads = [
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
genome_prob_sizes = {'g1': 0.5, 'g2': 0.5}
read_probabilities = {
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
	for read in reads:
		for genome_id, num_mismatches in read['mismatches'].items():
			this_read_genome_prob = genome_prob_sizes[genome_id] * ((sigma**num_mismatches) * ((1-sigma)**(read['read_length'] - num_mismatches)))
			all_read_genome_probs = 0
			for g_id, n_mis in read['mismatches'].items():
				all_read_genome_probs += genome_prob_sizes[g_id] * ((sigma**n_mis) * ((1-sigma)**(read['read_length'] - n_mis)))
			read_probabilities[read['id']][genome_id] = this_read_genome_prob / all_read_genome_probs
			
	# Update reference probability sizes
	for genome_id, size in genome_prob_sizes.items():
		prob_sum = 0
		for read_id, read_probs in read_probabilities.items():
			prob_sum += read_probs[genome_id]
		genome_prob_sizes[genome_id] = prob_sum / num_reads
	
	print('Genome prob sizes:')
	print(genome_prob_sizes)
	print()
	print('Read probabilities')
	print(read_probabilities)
	print()
	print()