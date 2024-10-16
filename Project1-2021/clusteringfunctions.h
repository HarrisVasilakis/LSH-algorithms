#ifndef _CLUSTERINGFUNCTIONS_H_
#define _CLUSTERINGFUNCTIONS_H_

#include <iostream>
#include <cstring>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <list>
#include <random>
#include <time.h>
#include "ghashfunction.h"

using namespace std;

int* readconfig(string path);
void kmeansplusplus(vector<float*> &centroids,vector<datasetarray*> datasetarraytest,int K,int* dim );


#endif