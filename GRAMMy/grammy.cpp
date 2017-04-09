/*
Austin Vo
2-4-2017
GRAMMy implementation in C++ v8
	Finally fixed getting data.
	But it should be really slow on large datasets because
	I am reading the file 3 times.

g++ -std=c++11 grammy.cpp -o grammy
	The -std flag is for using c++11 so that stoi() could work.
sudo apt-get install libboost-all-dev
	Use this to install Boost for C++ on Ubuntu.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <tuple>
#include <math.h>
#include <stdlib.h>
#include <algorithm> //std::unique, std::distance

#define CONVERGENCE 0.00000000001

using namespace std;

//Custom struct to hold metadata of genome references.
struct genome_reference{
	string name; //Proper name of genome.
	int taxid; //Taxonomic ID.
	string gbid; //ID for GenBank only.
	long length; //Length of genome in base-pairs.
	
	//Hierarchical taxonomy of bacteria.
	//Starts from current bacteria and ascends to higher level organisms.
	//String with semicolons (;) that separate each level.
	string taxonomy;
};

//Custom struct to hold data on reads from sam -> parsed csv file.
struct mapped_reads{
	string gbid; //ID for GenBank only.
	
	//If gbid == "read_length", val holds read length.
	//Else, val holds number of mismatches for matching gref.
	vector<int> val;
};

//Custom struct to hold genomic relative abundance (gra)
//for a particular taxid.
struct genome_rel_abund{
	int taxid; //Taxonomic ID.
	string name; //Name of organism.
	string taxonomy; //Taxonomy of organism.
	float rel_abund; //Relative abundance, in %.
	float err; //% error when checked against actual gra.
};

/*
Custom function that splits a string by a delimiter into
a vector of strings.
*/
void split(const string &s, char delim, vector<string> &elems) {
	stringstream ss;
	ss.str(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
}

//Helper function that calls custom split string function.
vector<string> split(const string &s, char delim) {
	vector<string> elems;
	split(s, delim, elems);
	return elems;
}

//Helper function to help find unique taxids.
bool compare_unique(int i, int j){
	return (i == j);
}

/*
Get the names of Genome References.
They should be the header row in the parsed csv file.
*/
vector <string> getGRefs(const char* inputfile){
	ifstream infile(inputfile);
	string line;
	vector <string> genome_refs;
	
	getline(infile, line);
	genome_refs = split(line, ','); // Split each line by commas.
	return genome_refs;
}

/*
Gets metadata of genome references from complete bacteria info.

Reads complete_bacteria_info.csv and aggregates data based on [Tax ID].
	For every unique [Tax ID], sum the [Molecule Length].
Ensures that I get the correct genome length for every unique bacteria.
	Accounts for genome references separated by chromosomes and plasmids.
Maintains book-keeping and accurate genetic data.

Returns a 2D int array in format of:
	[Tax ID][Complete Bacterial Genome Length]

Assumed column structure of complete_bacteria_info.csv:
	GenBank Account #, Molecule Type, Molecule Length, Tax ID, Organism Name
Assume that each row is completely unique.
*/
vector<genome_reference> getGRefMetaData(vector<string> grefs)
{
	ifstream infile("bacteria_summary.csv");
	string line; //Holds full line of csv file.
	vector<string> fields;
	vector<genome_reference> metadata; //Will return this.
	
	//For each gref, read through bacteria info for metadata.
	//Read file this way to keep sorted order of grefs.
	//Sorted order of grefs is important for grammy.
	for(int i = 0; i < grefs.size(); i++)
	{
		getline(infile,line); //Skips header row.
		while(getline(infile,line)) //Reads non-header rows.
		{
			//Parses line into individual fields.
			fields = split(line, ',');

			//Add genome reference metadata if there is a matching
			//[GenBank Account #] and gbid of current grefs.
			if(fields[0].compare(grefs[i]) == 0)
			{
				metadata.push_back(genome_reference());
				
//				metadata[metadata.size()-1].name = fields[4];
//				metadata[metadata.size()-1].taxid = stoi(fields[3]);
//				metadata[metadata.size()-1].gbid = fields[0];
//				metadata[metadata.size()-1].length = stol(fields[2]);
				
				metadata[i].name = fields[4];
				metadata[i].taxid = stoi(fields[3]);
				metadata[i].gbid = fields[0];
				metadata[i].length = stol(fields[2]);
				metadata[i].taxonomy = fields[5];
			}
		}

		//Reset fstream to actually read the data.
		infile.clear();
		infile.seekg(0, ios::beg);		
	}
	
	return metadata;
}

/*
Gets the reads of from the parsed SAM file.
Prepare data for GRAMMy algorithm.
Returns a vector in the format of:
	read_length, gref 1 mismatches, gref 1 mismatches, etc.
Need to know number of rows (reads) and columns (grefs) first.
This function assumes that the reference genome columns
are sorted from the python script.
*/
vector<mapped_reads> getMappedReads(const char* inputfile, vector<string>& grefs)
{
	ifstream infile(inputfile);
	string line;	
	vector<string> fields;	// Temp var for splitting strings.

	/*
	mr (for mapped reads) will hold the number of mismatches
	for each gref.
	mr[0] will hold the read lengths of each read.
	mr[1 to grefs.size()] will hold the number of mismatches
	for its corresponding gref.
	*/
	vector<mapped_reads> mr;
	
	//Add column dedicated to read lengths only.
	mr.push_back(mapped_reads());
	mr[mr.size()-1].gbid = "read_length";
	
	//Add columns for the rest of the grefs.
	//Ensures that gbid order of gref matches order of mr.
	for (int i = 0; i < grefs.size(); i++)
	{
		mr.push_back(mapped_reads());
		mr[i].gbid = grefs[i];
	}
	
	getline(infile, line); // Read and ignore header row.
	
	while(getline(infile, line)){	// Read parsed SAM csv file line by line.
		fields = split(line, ',');
		
		//Add each cell from infile to its matching cell in mr.
		//This works because of matching order of gbid.
		for(int c = 0; c < mr.size(); c++){
			try{
				mr[c].val.push_back(stoi(fields[c]));
			}
			catch (exception& e)
			{
				cout << "Standard exception inside getMappedReads():\n";
				cout << e.what() << endl;
				cout << "mr.size()\t" <<  mr.size() << endl;
				cout << "mr.capacity()\t" <<  mr.capacity() << endl;
				cout << "mr.max_size()\t" <<  mr.max_size() << endl;
				cout << "mr[" << c << "].val.size()\t" <<  mr[c].val.size() << endl;
				cout << "mr[" << c << "].val.capacity()\t" <<  mr[c].val.capacity() << endl;
				cout << "mr[" << c << "].val.max_size()\t" <<  mr[c].val.max_size() << endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	
	return mr; //All mr[].val vectors should have the same size.
}

/*
Executes GRAMMy framework from Xia et. al. paper.
Assume that grefs<> and gref_meta[].length sorted and match the same order.
*/
void grammy(const char* outputfile, vector<mapped_reads>& reads
	, vector <string>& grefs, vector<genome_reference>& gref_meta){
	double gref_prob, read_prob;

	//row = i, col = j in the GRAMMy paper.
	int row, col;
	
	//This is from Sigma paper.
	//Sigma is the uniform probility of assumed for any mismatch between a read and a genome.
	//May stem of genome variability of sequencing errors.
	//Defaults at 5%.
	double sigma = 0.05;
	
	//Instantiate data containers to match size of grefs<>.
	//Ensures that gbid order of gref matches order of Z.
	vector <double> PI; //Probability of read responsibility.
	vector <double> prevPI; //Tracks previous PI.
	vector <mapped_reads> Z; //Matrix indicating whether a read is from a corresponding genome.
	for (int i = 0; i < grefs.size(); i++)
	{
		// Initializing PI: size of probabilities.
		// Set each gref to be EQUALLY abundant in the sample.
		PI.push_back(1.0 / grefs.size());
		prevPI.push_back(0.0);
		
		//Create the columns of responsibility matrix.
		Z.push_back(mapped_reads());
		Z[i].gbid = grefs[i];
		
		//Create the rows of responsibility matrix.
		//Start at 1 to skip read_length column in reads<>.
		for (int j = 0; j < reads[i].val.size(); j++){
			Z[i].val.push_back(0);
		}
	}
	
	//-------------------------------------------------
	//EM Algorithm stage.
	int iter = 0; //Track number of iterations of EM.
	double rmse = 1.0; //Root Means Squared Error: for convergence.
	
	//Keep doing EM algorithm until RMSE is really low.
	while(rmse > CONVERGENCE)
	{
		//Always do EM at least once. Helps coverge.
		iter++; //Track number of iterations of EM.
	
		//E-step: Get probability matrix.
		//Implements Algorithm (3) in the GRAMMy paper for the E-step.
		//However, it uses Algorithm (1) in the Sigma paper to calculate
		//read probabilities instead of Algorithm (5) in GRAMMy paper.
		//Remember to skip read_length for reads.
		//Will manually call read_length as needed.
		//reads[col].val[row] = mismatch for a gref.
		//reads[0].val[row] = read length.
		for (col = 0; col < grefs.size(); col++)
		{
			for (row = 0; row < reads[col].val.size(); row++)
			{
				//This is the probability obtaining the read reads[c+1] with reads[c+1].val[row] mismatches in the alignment.

				//This is the probability of a gref being responsible for a read.
				//This is numerator of E-Step, which is Algorithm (5) of Sigma paper.
				gref_prob = PI[col] * ( 
					(pow(sigma, reads[col+1].val[row])) 
					* (pow(
						(1-sigma), 
						(reads[0].val[row] - reads[col+1].val[row])
					))
				);
				
				//This is the probability of the entire read.
				//Comes from denominator of E-step.
				read_prob = 0.0;
				for (int c = 0; c < grefs.size(); c++)
				{
					read_prob += PI[c] * ( 
						(pow(sigma, reads[c+1].val[row]))
						* (pow(
							(1-sigma), 
							(reads[0].val[row] - reads[c+1].val[row])
						))
					);
				}
				
				//Final responsibility Z.
				if(read_prob != 0)
				{
					try {
						Z[col].val[row] = gref_prob / read_prob;
					}
					catch (exception& e)
					{
						cout << "Standard exception inside grammy():\n";
						cout << e.what() << endl;
						cout << "Could not insert value into responsibility matrix Z.\n";
						exit(EXIT_FAILURE);
					}
				}
			}
		}
				
		//M-step: Update reference probability sizes.
		//This calculates Algorithm (4) in the GRAMMy paper.
		//Gets mixing coefficient for next iteration of EM.
		for (col = 0; col < grefs.size(); col++)
		{
			double prob_sum = 0.0;
			for (row = 0; row < reads[0].val.size(); row++){
				prob_sum += Z[col].val[row];
			}
			
			//Track current PI.
			prevPI[col] = PI[col];
			
			//Update current PI.
			PI[col] = prob_sum / reads[0].val.size();
		}
		
		//-------------------------------------------------
		//Convergence checking stage.
		
		//Check RMSE using current and previous PI values.  
		//Ensure that EM runs at least twice before checking RMSE.
		if (iter > 1)
		{
			double err = 0.0;
			for(int j = 0; j < grefs.size(); j++){
				err += pow(PI[j] - prevPI[j], 2);
			}
			rmse = sqrt(err);
		}
		
		cout << "\tIteration: " << iter << "\tConvergence: " << rmse << endl;
	}
	
	//-------------------------------------------------
	//Abundance Estimation stage.
	//Algorithm (1) from Xia paper.
	//Make sure that abundance[] sums up to 1.
	double abundance[grefs.size()];
	
	for (int j = 0; j < grefs.size(); j++)
	{
		//Calculate summation part of denominator.
		double summ = 0.0;
		for (int k = 0; k < grefs.size(); k++){
			summ += PI[k] / gref_meta[k].length;
		}
		
		abundance[j] = PI[j] / (gref_meta[j].length * summ);
	}
	
	//-------------------------------------------------
	//Aggregate unqieu organisms by summing relative abundance.
	
	//Copy and sort Tax IDs of the genome references.
	vector<int> taxids;
	for(int i = 0; i < gref_meta.size(); i++){
		taxids.push_back(gref_meta[i].taxid);
	}
	
	sort(taxids.begin(), taxids.end());
	
	//Find unique Tax IDs consecutively on sorted Tax IDs.
	vector<int>::iterator it;
	it = unique(taxids.begin(), taxids.end());
	taxids.resize(distance(taxids.begin(), it));
	unique(taxids.begin(), taxids.end(), compare_unique);
	
	//Aggregate sums by unique taxids.
	vector<genome_rel_abund> gra;
	for(int i = 0; i < taxids.size(); i++)
	{
		gra.push_back(genome_rel_abund());
		gra[i].taxid = taxids[i];
		
		//Loop through gref_meta for aggregation.
		for(int j = 0; j < grefs.size(); j++)
		{
			//For matching taxids:
				//sum relative abundance,
				//get the organism full name, and taxonomy.
			//Organisms under the same taxid also have
			//the same fill organism name and taxonomy.
			if(taxids[i] == gref_meta[j].taxid)
			{
				//Order of abundance must match gref_meta.
				gra[i].rel_abund += abundance[j];
				gra[i].name = gref_meta[j].name;
				gra[i].taxonomy = gref_meta[j].taxonomy;
			}
		}
	}
	
	//-------------------------------------------------
	
	//Output results of GRAMMy to GRA file.
	ofstream of (outputfile);
	
	//Output taxon id (gref name) on 1 line.
	for (col = 0; col < gra.size(); col++)
	{
		of << gra[col].name;
		
		if (col < gra.size()-1)
			of << "\t";
	}
	of << "\n";
	
	//Output taxonomy on 2 line.
	for (col = 0; col < gra.size(); col++)
	{
		of << gra[col].taxonomy;
		
		if (col < gra.size()-1)
			of << "\t";
	}
	of << "\n";
	
	//Output relative abundance of each gref on 3 line.
	//Remember, the data is sorted by gref name.
	for (col = 0; col < gra.size(); col++)
	{
		of << gra[col].rel_abund;
		
		if (col < gra.size()-1)
			of << "\t";
	}
	of << "\n";
}

int main (int argc, char* argv[])
{
	cout << "Start GRAMMy C++\n";
	const char* parsedFile = argv[1];
	
	cout << "Getting genome references from " << parsedFile << endl;
	vector<string> genome_refs;
	try {
		genome_refs = getGRefs(parsedFile);
	}
	catch (exception& e)
	{
		cout << "Standard exception when calling getGRefs():\n";
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	cout << "Getting mapped reads from " << parsedFile << endl;
	vector<mapped_reads> parsed_reads;
	try {
		parsed_reads = getMappedReads(parsedFile, genome_refs);
	}
	catch (exception& e)
	{
		cout << "Standard exception when calling getMappedReads():\n";
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	cout << "Getting genome reference metadata from bacteria_summary.csv\n";
	vector<genome_reference> gref_metadata;
	try {
		gref_metadata = getGRefMetaData(genome_refs);
	}
	catch (exception& e)
	{
		cout << "Standard exception when calling getGRefMetaData():\n";
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}

	//Convert filename arg from .csv to .sam file.
	string s(argv[1]);
	char gra[s.length()+1];
	memcpy(gra, &parsedFile[0], s.length()-4);
	gra[s.length()-4] = '.';
	gra[s.length()-3] = 'g';
	gra[s.length()-2] = 'r';
	gra[s.length()-1] = 'a';
	gra[s.length()] = '\0';
	
	//Execute GRAMMy calculations.
	cout << "Executing GRAMMy\n";
	try{
		grammy(gra, parsed_reads, genome_refs, gref_metadata);
	}
	catch (exception& e)
	{
		cout << "Standard exception when calling grammy():\n";
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}

	cout << "GRAMMy completed successfully\n";
	return 0;
}
