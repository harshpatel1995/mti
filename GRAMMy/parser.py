"""
Parser

Usage:
    parser.py <samfile>
"""

from docopt import docopt

# 1-29-2017

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

# Checks the integrity of the read.
# If a bad read is found change the tags to correctly parse.
# text: a line of text parsed from SAM file as a list.
def checkReadIntegrity(text):
    
    # If missing genome id at text[2] because of a bad read, get it from text[0] instead.
    if text[2] == '*':
        text[2] = text[0].split('_')[0]
    
    # If missing MD tag at text[12] because of a bad read, create MD tag using entire read of base-pairs from text[9].
    #This identifies bad read as total mismatch.
    if 'MD:Z:' not in text[12]:
        text[12] = 'MD:Z:' + text[9]
        
    return text
    
def parser():

    # Stores the information about the reference and the reads
    references = []   
    read_matrix = [] 
    
    # Open SAM file from BWA for parsing
    filename = options['<samfile>']
    with open(filename , 'r+') as f:
    
        for line in f:
            text = line.split()
            
            # Parse just the headers with the @SQ Tag : gives us the count and the names of the reference genomes
            try:
                if (text != [] and text[0] == "@SQ"):
                    try:
                        references.append(text[1][3:])
                    except:
                        print('Bad read for', text[0])
            except:
                print('Error for the following text in the SAM file:', text)
            
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
            if(text != [] and text[0] != "@SQ" and text[0] !=  "@PG"):
                text = checkReadIntegrity(text)
                
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
    
    # Output parsed reads into a tabular format
    f = open(options['<samfile>'][:-4] + ".csv", 'w')
    
    # Header row of table.
    # Sort and list out reference genomes for column names.
    f.write('read_length')
    f.write(',')
    for ref in sorted(list(references)):
        f.write(ref)
        f.write(',')
        
    f.write('\n')

    # List out reads, read length, and mismatches of each reference.
    for read in read_matrix:
        f.write(str(read['read_length']))
        
        # List out number of mismatches sorted by reference genomes.
        for gid, mm in sorted(list(read['mismatches'].items())):
            f.write(',')
            f.write(str(mm))
            
        f.write('\n')

def main():
    print('Parsing ' + options['<samfile>'] + ' into ' + options['<samfile>'][:-4] + ".csv")
    parser()
    print('Parsing completed.')
    print('\t'+ options['<samfile>'][:-4] + ".csv" + ' is ready for GRAMMy.')
    
options = docopt(__doc__, version='mti-vis 1.0')
main()

