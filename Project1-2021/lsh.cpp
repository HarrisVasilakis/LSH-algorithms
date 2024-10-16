#include <iostream>
#include <cstring>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <list>
#include <chrono>
#include "class.h"
#include "ghashfunction.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]) {
	int i,j,maxj,z,a,b,k=4,L=5,N=1,TableSize,perc=0,tperc=0;        //default values so you dont have to put everything in the command line
	int* dim=new int[2];
	int* dim2=new int[2];
	long int* gtest=new long int[2];
	ofstream outfile;
	double disttemp=-1;
	double* dist;
	bool flag1=true,flag2=true,flag3=true,flag4;
	float R=10000.0;
	string* neighbourname;
	string datapath,querypath,outputpath;
	vector<string> namebank;
	vector<datasetarray*> datasetarraytest;
	vector<datasetarray*> queryarray;
	vector<hashclass*> datasethashtable;
	if(argc>1){
		for(i=1; i<argc; i=i+2){
			if(strcmp(argv[i],"-i")==0){       //-i for path of dataset
				datapath=argv[i+1];
				flag1=false;
			}
			else if(strcmp(argv[i],"-q")==0){    //-q for path to query file
				querypath=argv[i+1];
				flag2=false;
			}
			else if(strcmp(argv[i],"-o")==0){    //-o for path to output file
				outputpath=argv[i+1];
				flag3=false;
			}
			else if(strcmp(argv[i],"-k")==0){   //-k for number of LSH functions
 				k=atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-L")==0){    //-L for number of hash tables
				L=atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-N")==0){    //-N for number of nearest neighbours
				N=atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-R")==0){    //-R for number of search radius 
				R=atof(argv[i+1]);
			}
			else{
				cout << "Wrong input\n";
				return -1;
			}
		}
	}
	if(flag1){
		cout << "Give dataset file path\n";   //give path of the dataset file 
		cin >> datapath;                  // e.g. :   dataset.txt
	}
	if(flag2){
		cout << "Give query file path\n";   
		cin >> querypath;                  
	}
	if(flag3){
		cout << "Give output file path\n";   
		cin >> outputpath;                  
	}
	vector<int> min(N,-1);
	vector<double> truedist(N,-1);
	outfile.open(outputpath);
	dim=readfile(datapath,datasetarraytest);  //here I read a file and save it in a vector of arrays of the class
	TableSize=dim[0]/32;                       ///////////////////////////////////////////////////////////////////////////Initialize vector   dim/4 is default
	auto start = high_resolution_clock::now();
	for (i=0; i<L; i++){
		hashclass* temphashclass=new hashclass(TableSize, k, dim[1]);
		datasethashtable.push_back(temphashclass);
	} 
	for(z=0; z<L; z++){
		for (i=0; i<dim[0]; i++){                 
			datasethashtable[z]->addtohashtable(*datasetarraytest[i], TableSize, k, dim[1], *datasethashtable[z]);
		}
		//datasethashtable[z]->printhashtable();
	}
	dim2 = readfile(querypath,queryarray);
	dist =new double[N*L];
	neighbourname = new string[N*L];      //make L arrays of int and string
	for (i=0; i<dim2[0]; i++){          //take each querry array
		namebank.clear();
		for(z=0;z<L;z++){               //check all of hashtables to find the best results
			for(j=0;j<N;j++){
				dist[N*z+j]=-1;              
			}
			gtest = ghashfunction(*queryarray[i],TableSize,k,dim[1],*datasethashtable[z]);
			if(datasethashtable[z]->ishashlistempty(gtest[1])){              //this is to ensure it goes in a bucket
				do{
					gtest[1]++;
					if(gtest[1]>=TableSize){
						gtest[1]=0;
					}
				}while(datasethashtable[z]->ishashlistempty(gtest[1]));
			}    
			for (auto temp:datasethashtable[z]->gethashlist(gtest[1])){        //check the bucket
				if (gtest[0]==temp.IDp){                      //check if they are close 
					disttemp = calcEuclideanDist(queryarray[i]->getcoordinates(), temp.dataset->getcoordinates(),dim[1]);   //distance
					flag4=true;
					maxj=-1;
					for(j=0;j<N*z;j++){                                             //here check if the item is already in the list from other hashtable
						if(temp.dataset->getname()==neighbourname[N*z+j]){
							flag4=false;
							break;
						}
					}
					if(flag4){
						flag4=false;
						for(j=0;j<N;j++){                                   //this for is for checking the N closest neighbours
							if(dist[N*z+j]==-1){                                     //now compare and save the closest neighbour so i can print them on the output file later
								dist[N*z+j]=disttemp;                             //if we haven't found N closest neighbours then count this one and stop the loop
								neighbourname[N*z+j]=temp.dataset->getname();
								flag4=false;
								break;
							}
							else if(disttemp<dist[N*z+j]){                              //if there is a closest neighbour that is further than the one we have
								flag4=true;
								if(maxj==-1){                                            //save him
									maxj=j;
								}
								else if(dist[N*z+j]>dist[N*z+maxj]){                     //and see if one of the closest neighbours we have has longer distance than the one we saved
									maxj=j;													//so if it does save him instead so we only save the closest of the closest neighbours
								}
							}
						}
						if(flag4){															//in the end if all the N closest neighbours positions are filled and we found an even closer neighbour
							dist[N*z+maxj]=disttemp;                                             //change the position
							neighbourname[N*z+maxj]=temp.dataset->getname();
						}
					}
				}
			}
			if(dist[N*z+N-1]==-1){
				for (auto temp:datasethashtable[z]->gethashlist(gtest[1])){        //if didn't find same IDp then do it again but check everything
					disttemp = calcEuclideanDist(queryarray[i]->getcoordinates(), temp.dataset->getcoordinates(),dim[1]);   //distance
					flag4=true;
					maxj=-1;
					for(j=0;j<N*z;j++){                                             //here check if the item is already in the list from other hashtable
						if(temp.dataset->getname()==neighbourname[N*z+j]){
							flag4=false;
							break;
						}
					}
					if(flag4){
						flag4=false;
						for(j=0;j<N;j++){                                   //this for is for checking the N closest neighbours
							if(dist[j]==-1){
								if(dist[N*z+j]==-1){                                     //now compare and save the closest neighbour so i can print them on the output file later
									dist[N*z+j]=disttemp;                             //if we haven't found N closest neighbours then count this one and stop the loop
									neighbourname[N*z+j]=temp.dataset->getname();
									flag4=false;
									break;
								}
								else if(disttemp<dist[N*z+j]){                              //if there is a closest neighbour that is further than the one we have
									flag4=true;
									if(maxj==-1){                                            //save him
										maxj=j;
									}
									else if(dist[N*z+j]>dist[N*z+maxj]){                     //and see if one of the closest neighbours we have has longer distance than the one we saved
										maxj=j;													//so if it does save him instead so we only save the closest of the closest neighbours
									}
								}
							}
						}
						if(flag4){															//in the end if all the N closest neighbours positions are filled and we found an even closer neighbour
							dist[N*z+maxj]=disttemp;                                             //change the position
							neighbourname[N*z+maxj]=temp.dataset->getname();
						}
					}
				}
				//outfile << "Problem with closest neighbour for " << queryarray[i]->getname() << "   " <<  disttemp <<  endl;
			}
			disttemp=-1;
		}
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		/*cout << i << endl;
		for(int hag=0;hag<N*L;hag++){
			cout << dist[hag] << ", ";
		}
		cout << endl;*/ ////testing
		//auto starttrue = high_resolution_clock::now();
		outfile << "Query: " << queryarray[i]->getname() << endl;
		for(a=0;a<N;a++){
			min[a]=-1;
			for(j=0;j<L;j++){         //find the minimum distance from all hashtables
				flag4=false;
				for(b=0;b<a;b++){                 //check if i have already included the distance
					if(dist[min[b]]==dist[j]){             ///what if some have the same distance???????????????
						flag4=true;
						break;
					}
				}
				if(flag4){
					continue;
				}
				if(dist[j]==-1){                     //fail proof
					continue;
				}
				if(min[a]==-1){                       //initialize
					min[a]=j;
				}
				else if(dist[min[a]]>dist[j] ){   //if i find shorter distance 
					min[a]=j;
				}
			}
			truedist[a]=-1;
			for(j=0;j<dim[0];j++){                    ///exhaustive search
				flag4=false;
				disttemp=calcEuclideanDist(queryarray[i]->getcoordinates(), datasetarraytest[j]->getcoordinates(),dim[1]);
				for(b=0;b<a;b++){
					if(truedist[b]==disttemp){
						flag4=true;
						break;
					}
				}
				if(flag4){
					continue;
				}
				if(truedist[a]==-1){
					truedist[a]=disttemp;
				}
				else if(truedist[a]>disttemp ){
					truedist[a]=disttemp;
				}
			}
			//cout << 100-i << endl;
			if(dist[min[a]]==truedist[a]){
				perc++;
			}
			tperc++;
			//cout << min[a] << " " << a << "out\n";
			//cout << neighbourname[min[a]] << endl;
			//cout << dist[min[a]] << endl;
			outfile << "Nearest neighbor-" << a+1 << ": " << neighbourname[min[a]] << endl;
			outfile << "distanceLSH: " << dist[min[a]] << endl;
			outfile << "distanceTrue: " << truedist[a] << endl;      //write in output file
			//auto stoptrue = high_resolution_clock::now();
			//auto durationtrue = duration_cast<microseconds>(stoptrue - starttrue);
			//cout << "out out\n";
		}
		outfile << "tLSH: " << duration.count()/1000 << endl;
		outfile << "tTrue: -" << endl;
		outfile << R << "-near neighbors: " << endl;           ///R range search
		namebank=rangesearch(L,TableSize, dim,datasethashtable,*queryarray[i],R,k);
		for(auto tempnamebank:namebank){
			outfile << tempnamebank << endl;
		}
		outfile << endl;
	}
	cout << "Success rate: " << perc*100/tperc <<endl;
	outfile.close();
}
