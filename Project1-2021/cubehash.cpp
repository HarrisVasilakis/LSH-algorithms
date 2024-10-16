#include <iostream>
#include <random>
#include <time.h>
#include <cstring>

#include "datasetarray.h"
#include "ghashfunction.h"
#include "class.h"
#include "cubehash.h"
#include "hypercubefunctions.h"

using namespace std;

// ============ CUBE HASH ============ //
cubehash::cubehash(int buckets, int k, int d) : lshhash(int buckets, int k, int d) {
	this->k = k;
}

int cubehash::ffunction(datasetarray dataset, int width) {
	float* v = this->lshhash.getRandomV();
	int t = this->lshhash.getRandomShiftT();

	srand(time(NULL));

    int hvalue = hFunction(dataset.getcoordinates(), v, t, w, width);

	// check if value already exists in map and if it does, return it
	if (this->valuemapping.find(hvalue) != this->valuemapping.end() ) {
		return this->valuemapping[hvalue];
	}

	int fvalue = rand() % 2;
	return fvalue;
}

string cubehash::hashvalue(datasetarray dataset) {
	// this calculates the value of a vertex (bucket) of the cube
	// there is a loop to calculate the k (=d') digits of the value
	// and then the digits need to be concatenated into a string
	
	int* hvalue;
  	hvalue = (int*)malloc(this->k*sizeof(int));

	for (int i = 0; i < this->k; i++) { // Creates an array in the form [0, 0, 0]
		hvalue[i] = (int) this->ffunction(dataset, i);
	}

	string hvaluestring = "";
    for (int i = 0; i < this->k; i++) { // Changes the array to string ex. "000"
        hvaluestring[i] += to_string(hvalue[i]);
    }

	return hvaluestring;
}

// ============ CUBE ENTRY ============ //
cubeentry::cubeentry(datasetarray dataset, int index) {
	this->dataset = dataset;
	this->index = index;
}

// ============ CUBE HASH TABLE ============ //
cubehashtable::cubehashtable() {}

cubehashtable::insert(datasetarray dataset, cubehash &cubehash) {
	string hashvalue = cubehash.hashvalue(dataset);
  	cubeentry *entry = new cubeentry(dataset, dataset.getname());

  	this->hashtable.insert(make_pair(hashvalue, entry));
	// This will be inserted even if the hashvalue already exists as a key,
	// because we need to be able to have multiple entries under the same hashvalue.
	// This is why we use Unordered Multimap instead of Unordered Map.

}

cubehashtable::~cubehashtable() {
 // clear hashtable structure
}