#include <iostream>
#include <cstring>
#include <vector>
#include "datasetarray.h"

using namespace std;

datasetarray::datasetarray(string name){
	//this->coordinates = (float*)malloc(sizeof(float)* dimension);
	this->name=name;
	
};

datasetarray::~datasetarray(){
	coordinates.clear();
};

void datasetarray::adddimension(int dimension){
	this->dimension=dimension;
};

void datasetarray::printname(){
	cout << name << "\n";
};

void datasetarray::addcoordinate(float coordinate){
	coordinates.push_back(coordinate);
	
};

void datasetarray::printarray(){
	int i;
	for (i=0;i<dimension;i++){
		cout << coordinates[i] << " ";
	}
	cout << "\n";
};

string datasetarray::getname(){
	return name;
};

vector<float> datasetarray::getcoordinates(){
	return coordinates;
};