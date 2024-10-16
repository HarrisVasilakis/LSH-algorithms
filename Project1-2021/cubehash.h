#ifndef _CUBE_HASH_H_
#define _CUBE_HASH_H_

#include <unordered_map>
#include <iostream>
#include <cstring>
#include "class.h"

using namespace std;

struct cubeentry {                  // This struct defines the structure of the entries of the hypercube hash table
    datasetarray dataset;           // The point that is added to the hash table
    int index;                      // The index that the point has in the input file ex. item_idN

    cubeentry(datasetarray, int);
};

class cubehash {
  public:
    int k;
    unordered_map<int, int> valuemapping; // map values produced by hfunction to random num (0/1) produced by ffunction
    hashclass lshhash;

    cubehash(int, int, int);

    int ffunction(datasetarray dataset, int);
    string hashvalue(datasetarray dataset);
};

class cubehashtable {
  public:
    unordered_multimap<string, cubeentry*> hashtable;

    /* Hashtable has this form:

        [   
            "000": {
                index: 0, 
                datasetarray: {...}
            },
            "001": [{ 
                index: 1, 
                datasetarray: {...}
            },
            { 
                index: 5, 
                datasetarray: {...}
            }],
            ...
        ]
    */
    cubehashtable();
    ~cubehashtable();

    void insert(datasetarray  , cubehash& );
};

#endif