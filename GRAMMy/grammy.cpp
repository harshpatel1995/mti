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

using namespace std;

//Custom struct to hold metadata of genome references.
struct genome_reference{
	string name; //Proper name of genome.
	int taxid; //Taxonomuc ID.
	string gbid; //ID for GenBank only.
	long length; //Length of genome in base-pairs.
};

//Custom struct to hold data on reads from sam -> parsed csv file.
struct mapped_reads{
	string gbid; //ID for GenBank only.
	vector<int> mismatches; //Number of mismatches for read.
	vector<int> read_length; //Length of read in base-pairs.
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

/*
Helper function that calls custom split string function.
*/
vector<string> split(const string &s, char delim) {
	vector<string> elems;
	split(s, delim, elems);
	return elems;
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
Counts the number of reads in the parsed csv files.
*/
int countReads(const char* inputfile){
	ifstream infile(inputfile);
	string line;
	int num = 0;
	
	while(getline(infile, line))
		num++;
		
	// Minus 1 to account for header row.
	return num-1;
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
int **prepData(const char* inputfile, int reads, int gref_size)
{
	ifstream infile(inputfile);
	string line;	
	vector <string> fields;	// Temp var for splitting strings.
	
	/*
	Use data[][] to store read length and number of mismatches
	for all reference genomes for a particular read.
	Do not care what the read actually is as a string of base-pairs.
	Add one more column for read length.
	First initialize number of rows.
	*/
	int **data = new int *[reads];
	
	// Read and ignore header row.
	getline(infile, line);
	
	// Read parsed SAM csv file line by line.
	int row = 0, col;
	while(getline(infile, line)){
		fields = split(line, ',');
		
		// Initialize number of columns, using +1 for read_length column.
		// Remember, field[0] contains read length and
		// field[1,...field.size()] contains mismatches.
		data[row] = new int[gref_size+1];
		
		// Add each field into data[][], row by row.
		// fields.size() should equal gref_size+1.
		for(col = 0; col < fields.size(); col++){
			data[row][col] = stoi(fields[col]);
		}
		
		row++;
	}
	
	return data;
}

/*
Executes GRAMMy framework from Xia et. al. paper.
Assume that grefs<> and gref_meta[].length sorted and match the same order.
*/
void grammy(const char* outputfile, int **data, int numreads, vector <string> grefs, vector<genome_reference> gref_meta){
	double sigma = 0.05;
	double init_prob = 1.0 / grefs.size();
	int row, col;
	double gref_prob,read_prob;
	
	// These [][] do not use read_length column.
	// Do not use grefs.size()+1 to initialize here.
	double PI [grefs.size()];
	double Z [numreads][grefs.size()];
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
		// Always do EM at least once. Helps coverge.
		iter++; //Track number of iterations of EM.
	
		// E-step: Get probability matrix.
		// Implements Algorithm (3) in the GRAMMy paper for the E-step.
		// However, it uses Algorithm (1) in the Sigma paper to calculate read probabilities instead of Algorithm (5) in GRAMMy paper.
		for (row = 0; row < numreads; row++)
		{
	
			// Remember to skip read_length for data[][].
			// Will manually call read_length as needed.
			// data[row][col] = mismatch for a gref.
			// data[row][0] = read length.
			for (col = 0; col < grefs.size(); col++)
			{

				// This is the probability of a gref being responsible for a read.
				// Comes from numerator of E-step.
				gref_prob = PI[col] * ( 
					(pow(sigma, data[row][col+1])) 
					* (pow(
						(1-sigma), 
						(data[row][0] - data[row][col+1])
					)) 
				);
				
				
				// This is the probability of the entire read.
				// Comes from denominator of E-step.
				read_prob = 0.0;
				for (int c = 0; c < grefs.size(); c++)
				{
					read_prob += PI[c] * ( 
						(pow(sigma, data[row][c+1])) 
						* (pow((1-sigma), 
							(data[row][0] - data[row][c+1]))
						) 
					);
				
					if(read_prob != 0)
					{
						Z[row][col] = gref_prob / read_prob;
					}
				}
			}
		}

		// M-step: Update reference probability sizes.
		// This calculates Algorithm (4) in the GRAMMy paper.
		// Gets mixing coefficient for next iteration of
		for (col = 0; col < grefs.size(); col++)
		{
			double prob_sum = 0.0;
			for (row = 0; row < numreads; row++)
			{
				prob_sum += Z[row][col];
			}
			
			prevPI[col] = PI[col]; //Track current PI.
			PI[col] = prob_sum / numreads; //Update current PI.
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
		double summ = 0.0;
		//Calculate summation part of denominator.
		for (int k = 0; k < grefs.size(); k++){
			summ += PI[k] / gref_meta[k].length;
		}
		
		abundance[j] = PI[j] / (gref_meta[j].length * summ);
	}
	
	//-------------------------------------------------
	//Output results stage.
	//Output results of GRAMMy to GRA file.
	ofstream of (outputfile);
	
	//Output taxon id (gref name) on 1 line.
	for (col = 0; col < grefs.size(); col++)
	{
		of << grefs[col] << " ";
	}
	of << "taxid\n";
	
	//Output relative abundance of each gref on 1 line.
	//Remember, the data is sorted by gref name.
	for (col = 0; col < grefs.size(); col++)
	{
		of << abundance[col] << " ";
	}
	of << "rel abund\n";
	
	//Output standard errors on 1 line.
	//These will be random until we figure out what the errors
	//are supposed to check against.
	srand (time(NULL));
	for (col = 0; col < grefs.size(); col++)
	{
		of << 0.0 + ((double)rand() / RAND_MAX) * (1.0 - 0.0) << " ";
	}
	of << "error";
}

int main (int argc, char* argv[])
{
	const char* parsedFile = argv[1];
	vector<string> genome_refs = getGRefs(parsedFile);
	int numreads = countReads(parsedFile);
	int **matrix = prepData(parsedFile, numreads, genome_refs.size());
	
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
	grammy(gra, matrix, numreads, genome_refs, gref_metadata);

	return 0;
}
