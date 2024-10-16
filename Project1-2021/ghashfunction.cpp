#include <iostream>
#include <cstring>
#include <vector>
#include <list>
#include <cmath>
#include <stdlib.h>     
#include <time.h>
#include <fstream>
#include <functional>
#include "ghashfunction.h"

using namespace std;

long int* ghashfunction(datasetarray dataset, int TableSize, int k, int width, hashclass hash) { 
	float* v = hash.getRandomV();
	int* r = hash.getRandomR();
	int t = hash.getRandomShiftT();
	// srand (time(NULL));          //random root
	long int M = 2;                 //big numbers so I need a log int or else I cant store
	for (int i = 1; i < 32; i++) {      //Make M = pow(2,32)-5
		M = M * 2; 
	}
	M = M - 5;
	long int g = 0;                  //start making g(p)
	for (int i = 0; i < k; i++) {	       // k is number for how many h functions
		
		g = g + r[i] * hFunction(dataset.getcoordinates(), v, t, w, width);		//k1h(p)+k2h(p)+......+kkH(p)
	}
	g= g % M; 
	long int* ID = new long int[2];
	ID[0] = g;
	g= g % TableSize; 
	if (g < 0) {         //for negative g
		g = -g;
	}
	// v.clear();
	ID[1] = g;
	//cout << g << endl;
	return ID;
}

long int* ghashfunction(float* dataset, int TableSize, int k, int width, hashclass hash) { 
	float* v = hash.getRandomV();
	int* r = hash.getRandomR();
	int t = hash.getRandomShiftT();
	// srand (time(NULL));          //random root
	long int M = 2;                 //big numbers so I need a log int or else I cant store
	for (int i = 1; i < 32; i++) {      //Make M = pow(2,32)-5
		M = M * 2; 
	}
	M = M - 5;
	long int g = 0;                  //start making g(p)
	for (int i = 0; i < k; i++) {	       // k is number for how many h functions
		
		g = g + r[i] * hFunction(dataset, v, t, w, width);		//k1h(p)+k2h(p)+......+kkH(p)
	}
	g= g % M; 
	long int* ID = new long int[2];
	ID[0] = g;
	g= g % TableSize; 
	if (g < 0) {         //for negative g
		g = -g;
	}
	// v.clear();
	ID[1] = g;
	//cout << g << endl;
	return ID;
}

float innerproduct(vector<float> p, float* v, int d) { // p vector, v random gen vector, d dimension of both vectors
	float inpr = 0;
	for (int i = 0; i < d; i++) {
		inpr = inpr + v[i] * p[i];
	}
	return inpr;
}

float innerproduct(float* p, float* v, int d) { // p vector, v random gen vector, d dimension of both vectors
	float inpr = 0;
	for (int i = 0; i < d; i++) {
		inpr = inpr + v[i] * p[i];
	}
	return inpr;
}


long int hFunction(vector<float> p, float* v, int t,  int window, int d) {
  float inpr = innerproduct(p, v, d);
  long int h = floor( (inpr + t) / window );
  return h;
}

long int hFunction(float* p, float* v, int t,  int window, int d) {
  float inpr = innerproduct(p, v, d);
  long int h = floor( (inpr + t) / window );
  return h;
}

double calcEuclideanDist(vector<float> v1, vector<float> v2,int dimension) {
  double sumOfSquares = 0;
  for (int i = 0; i < dimension; i++) {
    double diff = v1[i] - v2[i];
	sumOfSquares += diff * diff;
  }
  double dist = sqrt(sumOfSquares);
  return dist;
}

double calcEuclideanDist(vector<float*> v1,int j, vector<float> v2,int dimension) {
  double sumOfSquares = 0;
  for (int i = 0; i < dimension; i++) {
    double diff = v1[j][i] - v2[i];
	sumOfSquares += diff * diff;
  }
  double dist = sqrt(sumOfSquares);
  return dist;
}

double calcEuclideanDist(float* v1, vector<float> v2,int dimension) {
  double sumOfSquares = 0;
  for (int i = 0; i < dimension; i++) {
    double diff = v1[i] - v2[i];
	sumOfSquares += diff * diff;
  }
  double dist = sqrt(sumOfSquares);
  return dist;
}

int* readfile(string datapath,vector<datasetarray*>& datasetarraytest){       //function for reading files
	string text;											//I will use this for reading the lines of the file
	int i,j,length,s,width=0;
	char datasetstring[10];                                //I will use this for reading each word/number of the line
	fstream MyReadFile(datapath);                         //reads the file
	length=0;												  // k is size of lines in the file
	while (getline (MyReadFile, text)) {                  //one line at a time
  	  j=0;            
  	  s=0;   
  	  width=0;                    
  	  memset (datasetstring,0,10);
  	  for (i=0 ; i < (int)text.length(); i++){                   //then read the line by character to separate the tabs and get each set separate 
  	  	if (isspace(text[i])==0){
  	  		datasetstring[s]=text[i];               //then save it in an array , where the 0 position is the id name
  	  		s++;                                     //s counts the position of the char array for saving the words/numbers
  	  	}
  	  	else{
  	  		if(j==0){
  	  			datasetarray *d=new datasetarray(datasetstring);
  	  			datasetarraytest.push_back(d);
  	  		}
  	  		else{
  	  			datasetarraytest[length]->addcoordinate(atof(datasetstring));     //this is to turn a string into an int and add it to the allocated space of the class
  	  			width++;
  	  		}
  	  		s=0;
  	  		memset (datasetstring,0,10);                                            //delete space
  	  		j++;
  	  	}
  	  }
  	  datasetarraytest[length]->adddimension(width-1);
  	  length++;
	}
	datasetarraytest[length-1]->adddimension(width);  //it doesnt read the final \n so the last set comes off one value
	int* dim= new int[2];
	dim[0]=length;
	dim[1]=width;
	MyReadFile.close();
	return dim;
}

vector<string> rangesearch(int L, int TableSize, int* dim,vector<hashclass*> datasethashtable,datasetarray queryarray,float R,int k){
	int z,bv;
	double disttemp;
	bool flag5,flag4;
	long int* gtest=new long int[2];
	long int trueg;  
	vector<string> namebank,namelist;
	for(z=0;z<L;z++){       ///every hash table
		flag5 =true;
		bv=1;
		gtest = ghashfunction(queryarray,TableSize,k,dim[1],*datasethashtable[z]);
		trueg=gtest[1];            ///save g
		while(flag5){
			flag5=false;
			if(datasethashtable[z]->ishashlistempty(trueg)){              //this is to ensure it goes in a bucket
				do{
					trueg=trueg+bv;
					if(trueg>=TableSize){
						trueg=gtest[1]-1;
						bv=-1;
					}
					if(trueg<0){
						trueg=0;
						break;
					}
				}while(datasethashtable[z]->ishashlistempty(trueg));
			} 
			for (auto temp:datasethashtable[z]->gethashlist(gtest[1])){        //check the bucket
				disttemp = calcEuclideanDist(queryarray.getcoordinates(), temp.dataset->getcoordinates(),dim[1]);   //distance
				flag4=true;
				for(auto tempnamebank: namebank){                                             //here check if the item is already in the list of already checked names
					if(temp.dataset->getname()==tempnamebank){
						flag4=false;
						break;
					}
				}
				if(flag4){
					namebank.push_back(temp.dataset->getname());                                 //save name for future
					if(disttemp<=R){															//if inside range print
						flag5=true;
						namelist.push_back(temp.dataset->getname()); 
					}
				}
			}   
			trueg=trueg+bv;
			if(trueg<0){
				flag5=false;
			}
		}
	}
	return namelist;
}

vector<string> rangesearch(int L, int TableSize, int* dim,vector<hashclass*> datasethashtable,float* queryarray,float R,int k){
	int z,bv;
	double disttemp;
	bool flag5,flag4;
	long int* gtest=new long int[2];
	long int trueg;  
	vector<string> namebank,namelist;
	for(z=0;z<L;z++){       ///every hash table
		flag5 =true;
		bv=1;
		gtest = ghashfunction(queryarray,TableSize,k,dim[1],*datasethashtable[z]);
		trueg=gtest[1];            ///save g
		while(flag5){
			flag5=false;
			if(datasethashtable[z]->ishashlistempty(trueg)){              //this is to ensure it goes in a bucket
				do{
					trueg=trueg+bv;
					if(trueg>=TableSize){
						trueg=gtest[1]-1;
						bv=-1;
					}
					if(trueg<0){
						trueg=0;
						break;
					}
				}while(datasethashtable[z]->ishashlistempty(trueg));
			} 
			for (auto temp:datasethashtable[z]->gethashlist(gtest[1])){        //check the bucket
				disttemp = calcEuclideanDist(queryarray, temp.dataset->getcoordinates(),dim[1]);   //distance
				flag4=true;
				for(auto tempnamebank: namebank){                                             //here check if the item is already in the list of already checked names
					if(temp.dataset->getname()==tempnamebank){
						flag4=false;
						break;
					}
				}
				if(flag4){
					namebank.push_back(temp.dataset->getname());                                 //save name for future
					if(disttemp<=R){															//if inside range print
						flag5=true;
						namelist.push_back(temp.dataset->getname()); 
					}
				}
			}   
			trueg=trueg+bv;
			if(trueg<0){
				flag5=false;
			}
		}
	}
	return namelist;
}

vector<string> rangesearch(int L, int TableSize, int* dim,vector<hashclass*> datasethashtable,float* queryarray,float Rmin, float Rmax,int k){
	int z,bv;
	double disttemp;
	bool flag5,flag4;
	long int* gtest=new long int[2];
	long int trueg;  
	vector<string> namebank,namelist;
	for(z=0;z<L;z++){       ///every hash table
		flag5 =true;
		bv=1;
		gtest = ghashfunction(queryarray,TableSize,k,dim[1],*datasethashtable[z]);
		trueg=gtest[1];            ///save g
		while(flag5){
			flag5=false;
			if(datasethashtable[z]->ishashlistempty(trueg)){              //this is to ensure it goes in a bucket
				do{
					trueg=trueg+bv;
					if(trueg>=TableSize){
						trueg=gtest[1]-1;
						bv=-1;
					}
					if(trueg<0){
						trueg=0;
						break;
					}
				}while(datasethashtable[z]->ishashlistempty(trueg));
			} 
			for (auto temp:datasethashtable[z]->gethashlist(gtest[1])){        //check the bucket
				disttemp = calcEuclideanDist(queryarray, temp.dataset->getcoordinates(),dim[1]);   //distance
				flag4=true;
				for(auto tempnamebank: namebank){                                             //here check if the item is already in the list of already checked names
					if(temp.dataset->getname()==tempnamebank){
						flag4=false;
						break;
					}
				}
				if(flag4){
					namebank.push_back(temp.dataset->getname());                                 //save name for future
					if(disttemp<Rmax && disttemp>=Rmin){															//if inside range print
						flag5=true;
						namelist.push_back(temp.dataset->getname()); 
					}
				}
			}   
			trueg=trueg+bv;
			if(trueg<0){
				flag5=false;
			}
		}
	}
	return namelist;
}