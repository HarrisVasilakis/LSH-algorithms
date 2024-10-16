#ifndef _HYPERCUBE_FUNCTIONS_H_
#define _HYPERCUBE_FUNCTIONS_H_
#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

vector<string> findneighborvertices(string, int); // arguments are a vertex ex. "101010" and an int indicating the hamming distance
void fillvertices(string, int , int, vector<string> &);
char reversebit(char);


#endif