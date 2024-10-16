#include <iostream>
#include <cstring>
#include <vector>
#include <list>
#include <random>
#include <time.h>
#include "class.h"
#include "ghashfunction.h"

using namespace std;

float* nordist(int d) {
  float* v;
  v=(float*)malloc(d*sizeof(float));
  random_device device_random_;
  default_random_engine generator(device_random_()); 	//random root
  normal_distribution<float> distribution(0, 1);        //normal distribution function μ=0 σ=1
  // vector<float> v;
  for (int i = 0 ; i < d ; ++i) {
    v[i]= distribution(generator);   //generate random vector
    /*if (number<0){
      number=-number;
    }*/
    // v.push_back(number);
  }
  return v;
}

hashlist::hashlist(datasetarray& addingvalue){
	this->dataset=&addingvalue;
};

hashlist::~hashlist(){};

void hashlist::addIDp(long int IDp){         // this is (k1h(p)+k2h(p)+......+kkH(p))modM for when I put it in the hashtable
	this->IDp=IDp;
};

hashclass::hashclass(int buckets, int k, int d) {
	this->buckets=buckets;
	hashtable = new list<hashlist>[this->buckets];
	v=(float*)malloc(d*sizeof(float)); // make a v with a dimension of d
	v = nordist(d);
	srand(time(NULL));          //random root
	r=(int*)malloc(k*sizeof(int));
	for(int i = 0; i < k; i++) {
		r[i] = rand() % 2000000000; // approximate max int value
	}
	t = rand() % w;
};

hashclass::~hashclass() {
	hashtable->clear();
	//free(r);free(v);
};

float* hashclass::getRandomV() {
	return v;
}

int* hashclass::getRandomR() {
	return r;
}

int hashclass::getRandomShiftT() {
	return t;
}

void hashclass::addtohashtable(datasetarray& addingvalue, int TableSize, int k, int d, hashclass hash) {
	long int* gvalues = new long int[2];
	gvalues = ghashfunction(addingvalue, TableSize, k, d, hash);
	hashlist* temp = new hashlist(addingvalue);
	temp->addIDp(gvalues[0]);
	hashtable[gvalues[1]].push_back(*temp);
};

void hashclass::printhashtable(){   ////make it so i can see how much in every brucket
	int i,star[buckets]={0};
	for(i=0;i<buckets; i++){
		cout << "Bucket:" << i << endl;
		for(auto temp:hashtable[i]){
			temp.dataset->printname();
			//cout << temp.IDp << endl;
			temp.dataset->printarray();
			star[i]++;
		}
	}
	for(i=0;i<buckets;i++){      //number of datasets in each bucket
		cout << i << ": ";
		cout << star[i];
		cout << endl;
	}
};
void hashclass::printhashtablebucket(int position){
	//g(position)    ///here the hash function
	if (hashtable[position].empty()){
		cout << "there is nothing in this bucket\n";
	}
	else{
		for(auto temp:hashtable[position]){
			temp.dataset->printarray();
		}
	}
};

list<hashlist> hashclass::gethashlist(int position){
	return hashtable[position];
};

bool hashclass::ishashlistempty(int position){
	bool flag=true;
	for (auto temp:this->gethashlist(position)){
		flag=false;
		break;
	}
	return flag;
}