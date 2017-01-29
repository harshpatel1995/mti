/*
Austin Vo
1-29-2017
GRAMMy implementation in C++ v5
	Finally fixed getting data.
	But it should be really slow on large datasets because
	I am reading the file 3 times.

g++ -std=c++11 grammy_cpp_v4.cpp -o grammy
	The -std flag is for using c++11 so that stoi() could work.
sudo apt-get install libboost-all-dev
	Use this to install Boost for C++ on Ubuntu.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <math.h>

using namespace std;

/*
Custom function that splits a string by a delimiter into
a vector of strings.
*/
void split(const string &s, char delim, vector<std::string> &elems) {
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
	std::ifstream infile(inputfile);
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
	std::ifstream infile(inputfile);
	string line;
	int num = 0;
	
	while(getline(infile, line))
		num++;
		
	// Minus 1 to account for header row.
	return num-1;
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

void grammy(const char* outputfile, int **data, int numreads, vector <string> grefs){
	double sigma = 0.05;
	double init_prob = 1.0 / grefs.size();
	int row, col;
	double gref_prob,read_prob;
	
	// These [][] do not use read_length column.
	// Do not use grefs.size()+1 to initialize here.
	double PI [grefs.size()];
	double Z[numreads][grefs.size()];

	// Initializing PI: size of probabilities.
	// Set each gref to be EQUALLY abundant in the sample.
	for (col = 0; col < grefs.size(); col++)
		PI[col] = init_prob;
	
	// EM Algorithm.
	// Iterate 10 times. Will add threshold later.
	for (int i = 0; i < 2; i++)
	{
	
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
			
			PI[col] = prob_sum / numreads;
		}
	}
	
	ofstream of (outputfile);
	
	of << "grammy_cpp_v7\n";
	of << "Probability that genome references are responsible for sample:\n";
	
	for(col = 0; col < grefs.size(); col++)
	{
		of << grefs[col] << ": " << PI[col] << "\n";
	}

	of << "\nProbability that genome references are responsible for reads:\n";
	for(col = 0; col < grefs.size(); col++)
	{
		of << grefs[col] << "\t";
	}

	of << "\n";

	for(row = 0; row < numreads; row++)
	{
		for(col = 0; col < grefs.size(); col++)
		{
			of << Z[row][col] << "\t";
		}
		of << "\n";
	}
}

int main ()
{
	const char* parsedFile = "parsed_SAM_v4.csv";
	vector <string> genome_refs = getGRefs(parsedFile);
	int numreads = countReads(parsedFile);
	int **matrix = prepData(parsedFile, numreads, genome_refs.size());

	cout << "GRAMMy C++ v7\n";
	
	grammy("output_grammy_cpp_v7.txt", matrix, numreads, genome_refs);
	
	return 0;
}
