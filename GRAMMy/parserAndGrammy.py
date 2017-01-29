# 1-12-2017

# Metagenomic Taxonomic Inference (MTI) , University of Central Florida 2016-2017
# Austin Vo, Felix Sosa, Harsh Patel & Trevor Ballard
# Sponsor: Dr. Shibu Yooseph

# Parses the MD tag for a read mapped to a reference to find the number of mismatches
# A read may have multiple mappings to a reference (result of size); if so return the best match (least mismatches)
def parseMDTag(list):
    
    # If a read does not map to a reference, return -1 to indicate so.
    if(list == []):
        return -1
       
    mismatches = []
    
    # Else, return the least number of mismatches from the list
    for md_tag in list:
    
        # Extract the field value for the MD Tag i.e 54A15 for the example above
        md_field_value = md_tag.rsplit(':', 1)[1]
        # Count the number of upper case characters in the field  value i.e 1 for the example above 
        uppers = [l for l in md_field_value if l.isupper()]
        mismatches.append(len(uppers))

    return min(mismatches)  

    
def grammy(references, read_list):
    
    # Variables
    sigma = 0.05
    no_reads = len (read_list)
    no_references = len (references)
    
    fvo = open("check_pAG_v8.txt", 'w')
    
    # Initial probability that a reference is responsible for a read
    initial_probability = 1.0 / no_references
    
    # Mixing coefficient.
    # "Size" of probability of reference genome in relation to the sample size.
    # All values of PI should sum up to 1.
    PI = {}
    
    # Responsibility/Adjacency/Probability matrix.
    # Shows the probability that a reference genome is responsible for a read.
    # "Randomly" initialize each probability to have equal parts.
    Z = {}
    
    # Populate PI
    for reference in references:
        PI[reference] = initial_probability
    print(PI)
        
    # Populate Z
    for read in read_list:
        Z[read['id']] = {}
    for key in Z:
        for reference in references:
            Z[key][reference] = initial_probability
           
    # Iterating the EM step 10 times.
    for i in range(2):
        # E-step: Get probability matrix.
        # This implements Algorithm (3) in the GRAMMy paper for the E-step.
        # However, it uses Algorithm (1) in the Sigma paper to calculate read probabilities instead of Algorithm (5) in GRAMMy paper.
        
#        fvo.write(str(i))
#        fvo.write('\nE-STEP\n')
        
        for read in read_list:
            for genome_id, num_mismatches in read['mismatches'].items():
            
                this_read_genome_prob = PI[genome_id] * ((sigma**num_mismatches) * ((1-sigma)**(read['read_length'] - num_mismatches)))
                
                fvo.write('\t')
#                fvo.write(str(this_read_genome_prob))
#                fvo.write('\n')
#                
                all_read_genome_probs = 0
                
                for g_id, n_mis in read['mismatches'].items():
                    all_read_genome_probs += PI[g_id] * ((sigma**n_mis) * ((1-sigma)**(read['read_length'] - n_mis)))
                    
                    if(all_read_genome_probs != 0):
                        Z[read['id']][genome_id] = this_read_genome_prob / all_read_genome_probs

#                    fvo.write('\t\t')
#                    fvo.write(str(all_read_genome_probs))
#                    fvo.write('\t')
#                    fvo.write(str(Z[read['id']][genome_id]))
#                    fvo.write('\t')
#                    fvo.write(str(PI[g_id]))
#                    fvo.write('\t')
#                    fvo.write(str(read['read_length']))
#                    fvo.write('\t')
#                    fvo.write(str(n_mis))
#                    fvo.write('\n')
        
        fvo.write('\n---------------Z\n')
        fvo.write(str(Z))
        
        fvo.write('\n---------------PI\n')
        fvo.write(str(PI))
        
        fvo.write('\nM-STEP\n')
        # M-step: Update reference probability sizes.
        # This calculates Algorithm (4) in the GRAMMy paper.
        # Gets mixing coefficient for next iteration of EM.            
        for genome_id, size in PI.items():
            prob_sum = 0
            
            for read_id, read_probs in Z.items():
                prob_sum += read_probs[genome_id]
                fvo.write('\t')
                fvo.write(str(read_probs[genome_id]))
                fvo.write('\t')
                fvo.write(str(prob_sum))
                fvo.write('\n')
                
            PI[genome_id] = prob_sum / no_reads
            fvo.write('\t\t')
            fvo.write(str(PI[genome_id]))
            fvo.write('\n')

    # Write the results to grammyOutput.txt
    f = open("grammyOutput_v8.txt", 'w')
    
    f.write('Parser and GRAMMy v8\n')
    
    f.write('Genome prob sizes:\n')
    f.write(str(PI))
    
    f.write('\n\n\nRead probabilities: \n')
    
    f.write('Read,')
    f.write('\tProb_G1, ')
    f.write('\tMM_G1, ')
    f.write('\tProb_G2, ')
    f.write('\tMM_G2, ')
    f.write('\n')
    
    for r in read_list:
        
        #Read.
        f.write(str(r['id']))
        f.write(', ')
        
        for ref in references:
            #Prob_G#
            f.write('\t')
            f.write(str(Z[r['id']][ref]))
            f.write(', ')
            
            #MM_G#
            f.write('\t')
            f.write(str(r['mismatches'][ref]))
            f.write(', ')
        
        #f.write('\tMismatches = ')
        #f.write(str(r['mismatches']))
        f.write('\n')
    
    f.close()

# Uses the information about a read to generate and return a dictionary 
def helpParser(md_tags, read_info):
    
    for key, value in md_tags.items():
                    
        md_value = parseMDTag(value)
        # If no a read doesn't map to a reference, no of matches = 0 & no of mismatches = read length
        if(md_value == -1):
            read_info['mismatches'][key] = read_info['read_length']
        # Otherwise, use the least number of mismatches
        else:
            read_info['mismatches'][key] = md_value
        md_tags[key] = []
        
    temp  = {}
    for key in read_info['mismatches']:
        temp[key] = read_info['mismatches'][key]
        
    # Generate information about the reads
    read = {'id': read_info['id'], 'read_length' : read_info['read_length'], 'mismatches' : temp}
    
    return read
    
def parser():

    # Stores the information about the reference and the reads
    references = []   
    read_matrix = [] 
    
    # Open SAM file from BWA for parsing
    filename = "twoReferences.sam"        
    with open(filename , 'r+') as f:
    
        for line in f:
            text = line.split()
            
            # Parse just the headers with the @SQ Tag : gives us the count and the names of the reference genomes
            if (text[0] == "@SQ"):
                references.append(text[1][3:])
            
        # Stores the mismatch information for a read mapped to all reference genomes
        mismatch_dict = {}
        md_tags = {}
        for reference in references:
            mismatch_dict[reference] = 0
            md_tags[reference] = []
               
        # What a read should look like
        read_info = {'id' : "", 'read_length' : 0, 'mismatches' : mismatch_dict}
        
        # Rewind the file pointer to the beginning for another pass through the file
        f.seek(0,0)     
        
        for line in f:
                
            text = line.split()
            
            # Ignore the header files (we already looked at those)
            if(text[0] != "@SQ" and text[0] !=  "@PG"):
                            
                # If this is a new read 
                if len(text[9]) > 1:
                    # If this is not the very first read
                    if read_info['id'] != "":
                        #Add what we have already stored
                        read_matrix.append(helpParser(md_tags, read_info)) 
                    
                    # Store the ID (Read itself), Read Length and MD Tag 
                    read_info['id'] = text[9]
                    read_info['read_length'] = len(text[9])
                    md_tags[text[2]].append(text[12])
                
                # If this is not a new read
                elif len(text[9]) == 1:
                    # Add the MD Tag Information
                    md_tags[text[2]].append(text[12])  
        
        # We hit EOF, add the information about the last read too
        read_matrix.append(helpParser(md_tags, read_info)) 
 
    f.close()
    
    f = open("py_grammy_test_v8.txt", 'w')
    
    f.write('References\n')
    f.write(str(references))
    
    f.write('Read Matrix\n')
    f.write(str(read_matrix))
    
    grammy(references, read_matrix)

def main():
    parser()
    
main()


