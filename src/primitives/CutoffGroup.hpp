#ifndef _CUTOFFGROUP_H_
#define _CUTOFFGROUP_H_
#include "primitives/Atom.hpp"

class CutoffGroup{
public:
  
  CutoffGroup() { 
    haveTotalMass = false; 
    totalMass = 0.0;
  }
  
  void addAtom(Atom* atom) {cutoffAtomList.push_back(atom);}
  
  Atom* beginAtom(vector<Atom*>::iterator& i){
    i = cutoffAtomList.begin();
    return i != cutoffAtomList.end()? *i : NULL;
  }
  
  Atom* nextAtom(vector<Atom*>::iterator& i){
    i++;
    return i != cutoffAtomList.end()? *i : NULL;
  }
  
  double getMass(){
    vector<Atom*>::iterator i;
    Atom* atom;
    double mass;
    
    if (!haveTotalMass) {

      totalMass = 0;
      
      for(atom = beginAtom(i); atom != NULL; atom = nextAtom(i)){
        mass = atom->getMass();
        totalMass += mass;
      }

      haveTotalMass = true;
    }
    
    return totalMass;
  }
  
  void getCOM(double com[3]){

    vector<Atom*>::iterator i;
    Atom* atom;
    double pos[3];
    double mass;
    
    com[0] = 0;
    com[1] = 0;
    com[2] = 0;
    totalMass = getMass();    
    
    if (cutoffAtomList.size() == 1) {
      
      beginAtom(i)->getPos(com);
      
    } else {
      
      for(atom = beginAtom(i); atom != NULL; atom = nextAtom(i)){
        mass = atom->getMass();
        atom->getPos(pos);
        com[0] += pos[0] * mass;
        com[1] += pos[1] * mass;
        com[2] += pos[2] * mass;
      }
      
      com[0] /= totalMass;
      com[1] /= totalMass;
      com[2] /= totalMass;
    }
    
  }
  
  int getNumAtom() {return cutoffAtomList.size();}

  int getGlobalIndex() {return globalIndex;}
  void setGlobalIndex(int id) {this->globalIndex = id;}
private:
  vector<Atom*> cutoffAtomList;
  bool haveTotalMass;
  double totalMass;
  int globalIndex;

};

#endif
