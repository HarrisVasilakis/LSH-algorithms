#ifndef _GHASHFUNCTION_H_
#define _GHASHFUNCTION_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <list>
#include <cmath>
#include "class.h"
#include "datasetarray.h"

using namespace std;

#define w 20              //w = [2,6]

class datasetarray;
class hashclass;

long int* ghashfunction(datasetarray, int, int, int, hashclass);             // (dataset, tableSize, k /* number of h functions used */, width, hashclass)
long int* ghashfunction(float* dataset, int TableSize, int k, int width, hashclass hash);
long int hFunction(vector<float> p, float* v, int t,  int window, int d);    // ()
long int hFunction(float* p, float* v, int t,  int window, int d);
float innerproduct(vector<float> p, vector<float> v, int d);
float innerproduct(float* p, float* v, int d);
double calcEuclideanDist(vector<float> v1, vector<float> v2,int dimension);
double calcEuclideanDist(vector<float*> v1,int j, vector<float> v2,int dimension);
double calcEuclideanDist(float* v1, vector<float> v2,int dimension);
int* readfile(string datapath,vector<datasetarray*>& datasetarraytest);
vector<string> rangesearch(int L, int TableSize, int* dim,vector<hashclass*> datasethashtable,datasetarray queryarray,float R,int k);
vector<string> rangesearch(int L, int TableSize, int* dim,vector<hashclass*> datasethashtable,float* queryarray,float R,int k);
vector<string> rangesearch(int L, int TableSize, int* dim,vector<hashclass*> datasethashtable,float* queryarray,float Rmin, float Rmax,int k);   


#endif