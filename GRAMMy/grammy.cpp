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

using namespace std;

//Custom struct to hold metadata of genome references.
struct genome_reference{
	string name; //Proper name of genome.
	int taxid; //Taxonomic ID.
	string gbid; //ID for GenBank only.
	long length; //Length of genome in base-pairs.
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
	ifstream infile("complete_bacteria_info.csv");	
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
				metadata[metadata.size()-1].name = fields[4];
				metadata[metadata.size()-1].taxid = stoi(fields[3]);
				metadata[metadata.size()-1].gbid = fields[0];
				metadata[metadata.size()-1].length = stol(fields[2]);
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
vector<mapped_reads> getMappedReads(const char* inputfile, vector<string> grefs)
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
		mr[mr.size()-1].gbid = grefs[i];
	}
	
	getline(infile, line); // Read and ignore header row.
	
	// Read parsed SAM csv file line by line.
	while(getline(infile, line)){
		fields = split(line, ',');
		
		//Add each cell from infile to its matching cell in mr.
		//This works because of matching order of gbid.
		for(int c = 0; c < mr.size(); c++){
			mr[c].val.push_back(stoi(fields[c]));
		}
	}
	
	//All mr[].val vectors should have the same size.
	return mr;
}

/*
Executes GRAMMy framework from Xia et. al. paper.
Assume that grefs<> and gref_meta[].length sorted and match the same order.
*/
void grammy(const char* outputfile, vector<mapped_reads> reads, vector <string> grefs, vector<genome_reference> gref_meta){
	double sigma = 0.05;
	double init_prob = 1.0 / grefs.size();
	int row, col;
	double gref_prob,read_prob;
	
	//These [][] do not use read_length column.
	//All reads[].val.size() should be the same.
	double PI [grefs.size()];
	double Z [reads[0].val.size()][grefs.size()];
	double prevPI [grefs.size()];

	// Initializing PI: size of probabilities.
	// Set each gref to be EQUALLY abundant in the sample.
	for (col = 0; col < grefs.size(); col++)
	{
		PI[col] = init_prob;
		prevPI[col] = 0.0;
	}
	
	//-------------------------------------------------
	//EM Algorithm stage.
	int iter = 0; //Track number of iterations of EM.
	double rmse = 1.0; //Root Means Squared Error - use for convergence.
	
	//Keep doing EM algorithm until RMSE is really low.
	while(rmse > 0.00000000001)
	{
		//Always do EM at least once. Helps coverge.
		iter++; //Track number of iterations of EM.
	
		//E-step: Get probability matrix.
		//Implements Algorithm (3) in the GRAMMy paper for the E-step.
		//However, it uses Algorithm (1) in the Sigma paper to calculate
		//read probabilities instead of Algorithm (5) in GRAMMy paper.
		for (col = 0; col < grefs.size(); col++)
		{
	
			//Remember to skip read_length for data[][].
			//Will manually call read_length as needed.
			//reads[col].val[row] = mismatch for a gref.
			//reads[0].val[row] = read length.
			for (row = 0; row < reads[col].val.size(); row++)
			{

				//This is the probability of a gref being responsible for a read.
				//Comes from numerator of E-step.
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
						* (pow((1-sigma), 
							(reads[0].val[row] - reads[c+1].val[row]))
						) 
					);
				
					if(read_prob != 0){
						Z[row][col] = gref_prob / read_prob;
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
			for (row = 0; row < reads[0].val.size(); row++)
			{
				prob_sum += Z[row][col];
			}
			
			prevPI[col] = PI[col]; //Track current PI.
			PI[col] = prob_sum / reads[0].val.size(); //Update current PI.
		}
		
		//-------------------------------------------------
		//Convergence checking stage.
		
		//Check RMSE using current and previous PI values.  
		//Ensure that EM runs at least twice before checking RMSE.
		if (iter > 1)
		{
			rmse = 0.0;
			for(int j = 0; j < grefs.size(); j++){
				rmse += pow(PI[j] - prevPI[j], 2);
			}
			rmse = sqrt(rmse);
		}
		
		cout << iter << "\t" << rmse << "\n";
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
	//Aggregate Taxids by summing relative abundance.
	
	//Finding unique taxids.
	vector<int> taxids;
	for(int i = 0; i < gref_meta.size(); i++){
		taxids.push_back(gref_meta[i].taxid);
	}
	
	vector<int>::iterator it;
	it = unique(taxids.begin(), taxids.end());
	taxids.resize(distance(taxids.begin(), it));
	unique(taxids.begin(), taxids.end(), compare_unique);
	
	//Aggregate by unique taxids.
	vector<genome_rel_abund> gra;
	for(int i = 0; i < taxids.size(); i++)
	{
		gra.push_back(genome_rel_abund());
		gra[gra.size()-1].taxid = taxids[i];
		
		//Summing relative abundances for the same taxids.
		for(int j = 0; j < grefs.size(); j++)
		{
			if(taxids[i] == gref_meta[j].taxid)
			{
				//Order of abundance must match gref_meta.
				gra[gra.size()-1].rel_abund += abundance[j];
			}
		}
	}
	
	//-------------------------------------------------
	//Output results stage.
	//Output results of GRAMMy to GRA file.
	ofstream of (outputfile);
	
	//Output taxon id (gref name) on 1 line.
	for (col = 0; col < gra.size(); col++)
	{
		of << gra[col].taxid;
		
		if (col < gra.size()-1)
			of << "\t";
	}
	of << "\n";
	
	//Output relative abundance of each gref on 1 line.
	//Remember, the data is sorted by gref name.
	for (col = 0; col < gra.size(); col++)
	{
		of << gra[col].rel_abund;
		
		if (col < gra.size()-1)
			of << "\t";
	}
	of << "\n";
	
	//Output standard errors on 1 line.
	//These will be random until we figure out what the errors
	//are supposed to check against.
	srand (time(NULL));
	for (col = 0; col < gra.size(); col++)
	{
		of << 0.0 + ((double)rand() / RAND_MAX) * (1.0 - 0.0);
		
		if (col < gra.size()-1)
			of << "\t";
	}
	of << "\n";
}

int main (int argc, char* argv[])
{
	const char* parsedFile = argv[1];
	vector<string> genome_refs = getGRefs(parsedFile);
	vector<mapped_reads> parsed_reads = getMappedReads(parsedFile, genome_refs);
	vector<genome_reference> gref_metadata = getGRefMetaData(genome_refs);

	cout << "GRAMMy C++\n";

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
	grammy(gra, parsed_reads, genome_refs, gref_metadata);

	return 0;
}
