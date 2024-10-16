#ifndef _DATASETARRAY_H_
#define _DATASETARRAY_H_

#include <iostream>
#include <cstring>
#include <vector>

using namespace std;

class datasetarray{
		vector<float> coordinates;
		string name;
		int dimension;
	public:
		datasetarray(string);
		~datasetarray();
		string getname();
		vector<float> getcoordinates();
		void adddimension(int);
		void printname();
		void addcoordinate(float);
		void printarray();
};

#endif