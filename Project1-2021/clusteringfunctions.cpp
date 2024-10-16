#include <iostream>
#include <cstring>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <list>
#include <random>
#include <time.h>
#include "clusteringfunctions.h"
// #include "datasetarray.h"
// #include "class.h"
// #include "ghashfunction.h"

using namespace std;

int* readconfig(string path){
	int* arg = new int[6];
	int i,j;
	string text;
	char str[3];
	fstream MyReadFile(path);
	j=0;
	while (getline (MyReadFile, text)) { 
		memset (str,0,3);
		for (i=0 ; i<(int)text.length(); i++){ 
			if (isspace(text[i])!=0){
				str[0]=text[i+1];
				if(isspace(text[i+2])==0){
					str[1]=text[i+2];
				}
				else if(isspace(text[i+3])==0){
					str[1]=text[i+2];
					str[2]=text[i+3];
				}
				arg[j]=atoi(str);
				j++;
				break;
			}
		}
	}
	return arg;
}

void kmeansplusplus(vector<float*> &centroids,vector<datasetarray*> datasetarraytest,int K,int* dim ){
	int i,j,z;
	double temp;
	double P[dim[0]]={-1};
	float* tempcentroid=new float[dim[1]];
	srand(time(NULL)); 
	int r=rand()%dim[0];                                        ///choose on random first centroid
	if(r<0){r=-r;}
	for(j=0;j<dim[1];j++){
		tempcentroid[j]=datasetarraytest[r]->getcoordinates()[j];
	}
	centroids.push_back(tempcentroid);
	for(z=0;z<K;z++){
		for(i=0;i<dim[0];i++){
			for(j=0;j<z;j++){
				temp =calcEuclideanDist(centroids,j,datasetarraytest[i]->getcoordinates(),dim[1]);               ///make the distances
				//cout << centroids[j][dim[1]-1] << " " << datasetarraytest[i]->getcoordinates()[dim[1]-1] << endl;
				if(P[i]==-1){
					P[i]=temp;
				}
				else if(P[i]>temp){
					P[i]=temp;
				}
			}
			P[i]=sqrt(P[i]);
			if(i>0){
				P[i]=P[i]+P[i-1];
			}

		}
		r=remainder((double)rand(),P[dim[0]-1]);               ///pick a random number
		if(r<0){r=-r;}
		for(i=0;i<dim[0];i++){
			if(P[i]>=r){
				float* tempcentroid=new float[dim[1]];             //initialize every time or it gets overwritten
				for(j=0;j<dim[1];j++){                             ///get next random centroid
					tempcentroid[j]=datasetarraytest[i]->getcoordinates()[j];
				}
				centroids.push_back(tempcentroid);
				break;
			}
		}
		for(j=0;j<dim[0];j++){
			P[j]=-1;
		}
	}
}