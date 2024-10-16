#include <iostream>
#include <cstring>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <list>
#include "class.h"
#include "ghashfunction.h"
#include "cubehash.h"

using namespace std;

int main(int argc, char *argv[]) {
	int i, width = 128; // j, maxj, z, a, b, , TableSize, perc = 0, tperc = 0;
	
	// assign default values
	int k, probes, M, N; // int k=14, probes=2, M=10, N=1;
	float R; // float R=10000.0;

	// int* dim = new int[2];
	// int* dim2 = new int[2];
	// long int* gtest = new long int[2];
	
	// double disttemp=-1;
	// double* dist;
	bool inputfile = true, queryfile = true, outputfile = true; // , flag4;
	
	// string* neighbourname;
	string datapath, querypath, outputpath;
	vector<datasetarray*> datasetarraytest;
	vector<datasetarray*> queryarray;
	vector<hashclass*> datasethashtable;

	// read command line arguments
	if(argc > 1) {
		for(i=1; i<argc; i=i+2){
			if(strcmp(argv[i],"-i")==0) {       //-i for path of dataset
				datapath=argv[i+1];
				inputfile = false;
			}
			else if(strcmp(argv[i],"-q")==0) {    //-q for path to query file
				querypath=argv[i+1];
				queryfile = false;
			}
			else if(strcmp(argv[i],"-o")==0) {    //-o for path to output file
				outputpath=argv[i+1];
				outputfile = false;
			}
			else if(strcmp(argv[i],"-d")==0) {   //-d for dimensions. that means the number of coordinates in one line without counting the name id
				width=atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-k")==0) {   //-k (d') for dimension
 				k=atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-M")==0) {    //-M for max number of points to be checked
				M=atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-probes")==0) {    //-probes for max number of hypercube's edges
				probes=atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-N")==0) {    //-N for number of nearest neighbours
				N=atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-R")==0) {    //-R for number of search radius 
				R=atof(argv[i+1]);
			}
			else{
				cout << "Wrong input\n";
				return -1;
			}
		}
	}
	if(inputfile) {
		cout << "Give dataset file path\n";   //give path of the dataset file 
		cin >> datapath;                  // e.g. :   dataset.txt
	}
	if(queryfile) {
		cout << "Give query file path\n";   
		cin >> querypath;                  
	}
	if(outputfile) {
		cout << "Give output file path\n";   
		cin >> outputpath;                  
	}
	
}