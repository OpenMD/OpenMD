#include <iostream>
#include <stdlib.h>

#include "Exclude.hpp"

Exclude* Exclude::_instance = 0;

Exclude* Exclude::Instance() {
  if (_instance == 0) {
    _instance = new Exclude;
  }
  return _instance;
}

Exclude::Exclude(){   
  exPairs = NULL;
  newFortranArrayNeeded = 1;
}

Exclude::~Exclude() {
  if (exPairs != NULL) {
    delete[] exPairs;
  }
  delete _instance;
}
  
int* Exclude::getFortranArray(){
  
  set<pair<int, int> >::iterator  i;
  int j;

  if (newFortranArrayNeeded != 0) {
    delete[] exPairs;
    exPairs = new int[2*getSize()];  
    j = 0;
    for(i = excludeSet.begin(); i != excludeSet.end(); ++i) {   
      exPairs[j] = (*i).first;
      j++;
      exPairs[j] = (*i).second;
      j++;
    }
    newFortranArrayNeeded = 0;
  }
 
  return exPairs;  
}


void Exclude::addPair(int i, int j) {
  
  if (!hasPair(i, j)) {
 
    if (i != j) {
      
      if (i < j) 
        excludeSet.insert(make_pair(i, j));
      else 
        excludeSet.insert(make_pair(j, i));
    }

    newFortranArrayNeeded = 1;
  }

}


void Exclude::printMe( void ){

  set<pair<int, int> >::iterator  i;
  int index;
  
  index = 0;
  for(i = excludeSet.begin(); i != excludeSet.end(); ++i) {   

    std::cerr << "exclude[" << index << "] i, j: " << (*i).first << " - " << (*i).second << "\n";
    index++;
  
  }  
}

int Exclude::hasPair(int i, int j) {

  set<pair<int, int> >::iterator  position;

  if (i != j) {   
    if (i < j) 
      position = excludeSet.find(make_pair(i, j));
    else 
      position = excludeSet.find(make_pair(j, i));

    if (position != excludeSet.end()) 
      return 1;
    else
      return 0;
  } else 
    return 0;
}

int Exclude::getSize() {
  return excludeSet.size();
}
