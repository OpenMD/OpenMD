#include <iostream>
#include <stdlib.h>

#include "brains/SkipList.hpp"

SkipList* SkipList::_instance = 0;

SkipList* SkipList::Instance() {
  if (_instance == 0) {
    _instance = new SkipList;
  }
  return _instance;
}

SkipList::SkipList(){   
}

SkipList::~SkipList() {
  delete _instance;
}
  
void SkipList::addAtom(int i) {
  
  if (!hasAtom(i)) 
    skipSet.insert(i);

}


void SkipList::printMe( void ){

  set<int>::iterator  i;
  int index;
  
  index = 0;
  for(i = skipSet.begin(); i != skipSet.end(); ++i) {   

    std::cerr << "SkipList[" << index << "] i: " << *i << "\n";
    index++;  
  }  
}

int SkipList::hasAtom(int i) {

  set<int>::iterator  position;

  position = skipSet.find(i);
  
  if (position != skipSet.end()) 
    return 1;
  else
    return 0;

}

int SkipList::getSize() {
  return skipSet.size();
}
