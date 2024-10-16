#include <iostream>
#include <vector>
#include "hypercubefunctions.h"

using namespace std;

// ============= FIND NEIGHBOR VERTICES ============= //
vector<string> findneighborvertices(string vertex, int hammingdist) { 
    // This function returns a vector of the vertices that have max hamming distance = hammingdist from the vertex we provided
    vector<string> vertices;

    fillvertices(vertex, hammingdist, vertex.size() - 1,  vertices);

    return vertices;
}

// ============= FILL VERTICES VECTOR RECURSIVELY============= //
void fillvertices(string vertex, int diff, int index,  vector<string> &vertices) {
  if (diff == 0) {                  // zero hamming distance = same vertex
    vertices.push_back(vertex);
    return;
  }
  if (index < 0) return;           // no more bits, all checked

  vertex[index] = reversebit(vertex[index]);                    
  fillvertices(vertex, diff - 1, index - 1,  vertices);

  vertex[index] = reversebit(vertex[index]);
  fillvertices(vertex, diff, index - 1,  vertices);
}

char reversebit(char bit) { 
    if ( bit == '1' ) {
        return '0';
    } else if ( bit == '0' ) {
        return '1';
    }
     // 1 -> 0, 0 -> 1
} 