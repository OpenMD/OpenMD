#include "Utility.hpp"
#include <vector>
#include <iostream>
#include <math.h>
using namespace std;

double dotProduct(vector<double>& v1, vector<double>& v2){
  double sum;

  sum = 0;

  if(v1.size() != v2.size()){
    cerr << "Utility Error: dimension of two vectors are not matched" << endl;
    exit(-1);
  }
  
   for(int  i = 0; i < v1.size(); i++)
     sum += v1[i]*v2[i];
   return sum;
}

double norm2(vector<double>& x){
  return sqrt(dotProduct(x, x));
}	
