#include <iostream>
#include <cstring>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <list>
#include <random>
#include <time.h>
#include <chrono>
#include "clusteringfunctions.h"
//#include "class.h"
//#include "ghashfunction.h"
//#include "datasetarray.h" 

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]) {
	int i,z,j,TableSize,c,min;
	int* arg=new int[6];                                    ///{K,L,k,M,k,probes}
	bool flag1=true,flag2=true,flag3=true,flag4=true,flag5=true,flag6=true,flag7=false;
	string datapath,configpath,outputpath,method,dataname;
	ofstream outfile;
	double distance,mindistance;
	float R;
	vector<datasetarray*> datasetarraytest;
	vector<hashclass*> datasethashtable;
	vector<float> tempcoordinates;
	srand(time(NULL)); 
	int* dim=new int[2];
	if(argc>1){
		for(i=1; i<argc; i=i+2){
			if(strcmp(argv[i],"-i")==0){       //-i for path of dataset
				datapath=argv[i+1];
				flag1=false;
			}
			else if(strcmp(argv[i],"-c")==0){    //-c for path to configuration file
				configpath=argv[i+1];
				flag2=false;
			}
			else if(strcmp(argv[i],"-o")==0){    //-o for path to output file
				outputpath=argv[i+1];
				flag3=false;
			}
			else if(strcmp(argv[i],"-m")==0){    //-m for method
				method=argv[i+1];
				flag4=false;
			}
			else if(strcmp(argv[i],"-complete")==0){    //-complete for printing
				flag7=true;
			}
		}
	}
	if(flag1){
		cout << "Give dataset file path\n";   //give path of the dataset file 
		cin >> datapath;                  // e.g. :   dataset.txt
	}
	if(flag2){
		cout << "Give configuration file path\n";   
		cin >> configpath;                  
	}
	if(flag3){
		cout << "Give output file path\n";   
		cin >> outputpath;                  
	}
	if(flag4){
		cout << "Give method\n";   
		cin >> method;                  
	}
	outfile.open(outputpath);
	dim=readfile(datapath,datasetarraytest);         ///makw dataset array
	vector<float*> centroids;
	float* newcentroid = new float[dim[1]];
	float* oldcentroid = new float[dim[1]];
	arg=readconfig(configpath);
	auto start = high_resolution_clock::now();
	kmeansplusplus(centroids,datasetarraytest,arg[0],dim );
	vector<string> clusterednames[arg[0]],totalclusterednames;
	vector<string> clusteredvaluesnames[arg[0]];
	if(method=="lsh"){
		TableSize=dim[0]/32;
		for (i=0; i<arg[1]; i++){
		hashclass* temphashclass=new hashclass(TableSize, arg[2], dim[1]);              ///make hashclass for lsh
		datasethashtable.push_back(temphashclass);
		} 
		for(z=0; z<arg[1]; z++){
			for (i=0; i<dim[0]; i++){                 
				datasethashtable[z]->addtohashtable(*datasetarraytest[i], TableSize, arg[2], dim[1], *datasethashtable[z]);
			}
		}
		while(flag5){
			R=0;
			while(flag6){
				c=0;
				for(auto tempos:centroids){
					clusterednames[c].clear();
					clusterednames[c]=rangesearch(arg[1],TableSize,dim,datasethashtable,tempos,R,R+10,arg[2]);       //search for range [i,i+x]
					for(auto name:clusterednames[c]){
						for(auto totalnames:totalclusterednames){
							if(totalnames==name){                                                                                    //if i find a name that has been counted before end lsh implementation
								flag6=false;
								break;                                                                                         
							}
						}
						totalclusterednames.push_back(name);
					}
					c++;

				}
				if(flag6){
					for(c=0;c<arg[0];c++){
						for(auto name:clusterednames[c]){
							clusteredvaluesnames[c].push_back(name);                  //here i put the range names in the right name array of the right cluster
						}
					}
				}
				if(flag6){
					R=R+10;
				}
			}
			flag6=true;
			do{
				totalclusterednames.clear();
				for(auto tempos2:centroids){
					clusterednames[0].clear();
					clusterednames[0]=rangesearch(arg[1],TableSize,dim,datasethashtable,tempos2,R,R+100,arg[2]); //like lsh get the range 
					for(auto anothername:clusterednames[0]){                 
						totalclusterednames.push_back(anothername);
					}
				}
				for(auto name:totalclusterednames){                                                               //but for the names you get , do loyd
					tempcoordinates=datasetarraytest[stoi(name)-1]->getcoordinates();
					for(c=0;c<arg[0];c++){                              //and for every centroid
						if(c==0){
							mindistance=calcEuclideanDist(centroids,c,tempcoordinates,dim[1]);           //count distance and put it in the cluster that has the minimum distance from the centroid
							min=c;
						}
						else{
							distance=calcEuclideanDist(centroids,c,tempcoordinates,dim[1]);
							if(mindistance>distance){
								mindistance=distance;
								min=c;
							}
						}
					}
					for(auto clusname:clusteredvaluesnames[min]){
						if(clusname==name){
							flag6=false;
						}
					}
					if(flag6){
						clusteredvaluesnames[min].push_back(name);
					}
					flag6=true;
				}
				R=R+100;
			}while(totalclusterednames.empty()==false);                ///do this until range returns nothing
			flag6=true;
			for(c=0;c<arg[0];c++){                              //for every centroid                      
				for(i=0;i<dim[1];i++){                         ///make a new centroid
					newcentroid[i]=0;
					oldcentroid[i]=centroids[c][i];
				}
				
				for(auto datasetname:clusteredvaluesnames[c]){
					tempcoordinates.clear();
					tempcoordinates=datasetarraytest[stoi(datasetname)-1]->getcoordinates();
					for(j=0;j<dim[1];j++){
						newcentroid[j]=newcentroid[j]+tempcoordinates[j];                       //that is the mean of all the datasets in the cluster
						
					}
					
				}
				for(j=0;j<dim[1];j++){
					newcentroid[j]=newcentroid[j]/(float)clusteredvaluesnames[c].size();
				}
				for(j=0;j<dim[1];j++){  
					centroids[c][j]=newcentroid[j];
				}
				if(flag6){
					for(j=0;j<dim[1];j++){                                //see if the new centroid changes or else end function
						if(centroids[c][j]==oldcentroid[j]){
							flag5=false;                  //it must find that everything has remained the same
						}
						else{
							flag5=true;                   //flag5 is to end the algorithm if the centroids don't change
							flag6=false;                  ///flag6 is if it finds one that has changed then ensure that the whole loop will continue
							break;
						}
					}
				}
			}
		}
	}
	else if(method=="loyd"){
		while(flag5){                //flag5
			for(i=0;i<arg[0];i++){
				clusteredvaluesnames[i].clear();               ///clear the list so I can put new names in the cluster as the centroids change
			}
			for(auto data:datasetarraytest){                      //for every dataset
				for(i=0;i<arg[0];i++){                              //and for every centroid
					if(i==0){
						mindistance=calcEuclideanDist(centroids,i,data->getcoordinates(),dim[1]);           //count distance and put it in the cluster that has the minimum distance from the centroid
						min=i;
						dataname=data->getname();
					}
					else{
						distance=calcEuclideanDist(centroids,i,data->getcoordinates(),dim[1]);
						if(mindistance>distance){
							mindistance=distance;
							min=i;
							dataname=data->getname();
						}
					}
				}
				clusteredvaluesnames[min].push_back(dataname);
			}
			flag6=true;
			for(c=0;c<arg[0];c++){                              //for every centroid                      
				for(i=0;i<dim[1];i++){                         ///make a new centroid
					newcentroid[i]=0;
					oldcentroid[i]=centroids[c][i];
				}
				
				for(auto datasetname:clusteredvaluesnames[c]){
					tempcoordinates.clear();
					tempcoordinates=datasetarraytest[stoi(datasetname)-1]->getcoordinates();
					for(j=0;j<dim[1];j++){
						newcentroid[j]=newcentroid[j]+tempcoordinates[j];                       //that is the mean of all the datasets in the cluster
						
					}
					
				}
				for(j=0;j<dim[1];j++){
					newcentroid[j]=newcentroid[j]/(float)clusteredvaluesnames[c].size();
				}
				for(j=0;j<dim[1];j++){  
					centroids[c][j]=newcentroid[j];
				}
				if(flag6){
					for(j=0;j<dim[1];j++){                                //see if the new centroid changes or else end function
						if(centroids[c][j]==oldcentroid[j]){
							flag5=false;                  //it must find that everything has remained the same
						}
						else{
							flag5=true;                   //flag5 is to end the algorithm if the centroids don't change
							flag6=false;                  ///flag6 is if it finds one that has changed then ensure that the whole loop will continue
							break;
						}
					}
				}
			}
			// for(c=0;c<arg[0];c++){ 
			// 	for(j=0;j<dim[1];j++){  
			// 		cout << centroids[c][j] << " ";
			// 	}
				
			// 	cout <<endl;
			// }
		}
	}
	else{
		cout << "There is no such method \n";
		return 0;
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	for(c=0;c<arg[0];c++){ 
		outfile << "CLUSTER-" << c << " {size:" << clusteredvaluesnames[c].size() << ", centroid: ";
		for(j=0;j<dim[1];j++){  
			outfile << centroids[c][j] << " ";
		}
		outfile << endl;
	}
	outfile << "clustering_time: " << duration.count()/1000 << endl;
	double a;
	double b;
	double s;
	double stotal=0;
	outfile << "Silhouette: [";
	for(c=0;c<arg[0];c++){ 
		s=0;
		for(auto namecl:clusteredvaluesnames[c]){
			a=0;
			b=0;
			for(auto namecl2:clusteredvaluesnames[c]){
				a=a+calcEuclideanDist(datasetarraytest[stoi(namecl2)-1]->getcoordinates(),datasetarraytest[stoi(namecl)-1]->getcoordinates(),dim[1]);
			}
			for(i=0;i<arg[0];i++){
				if(i!=c){
					if(i==0 || (c==0 && i==1)){
						mindistance=calcEuclideanDist(centroids,i,datasetarraytest[stoi(namecl)-1]->getcoordinates(),dim[1]); //mindistance here is actually second minimum distance
						min=i;
					}
					else{
						distance=calcEuclideanDist(centroids,i,datasetarraytest[stoi(namecl)-1]->getcoordinates(),dim[1]);
						if(mindistance>distance){
							mindistance=distance;
							min=i;
						}
					}
				}
			}
			for(auto namecl3:clusteredvaluesnames[min]){   
				b=b+calcEuclideanDist(datasetarraytest[stoi(namecl3)-1]->getcoordinates(),datasetarraytest[stoi(namecl)-1]->getcoordinates(),dim[1]);
			}  
			a=a/(clusteredvaluesnames[c].size()-1);
			b=b/clusteredvaluesnames[min].size();
			if(a>b){
				s=s+((b-a)/a);
			}       
			else{
				s=s+((b-a)/b);
			}
		}
		s=s/clusteredvaluesnames[c].size();
		outfile << s << ",";
		stotal=stotal+s;
	}
	stotal=stotal/arg[0];
	outfile << stotal << "]\n";
	if(flag7){                ///-complete
		for(c=0;c<arg[0];c++){ 
			outfile << "CLUSTER-" << c << " {";
			for(i=0;i<(int)clusteredvaluesnames[c].size();i++){
				outfile << clusteredvaluesnames[c][i] << ", ";
			}
			outfile<< "} \n";
		}
	}
	outfile.close();
}
