#include "UseTheForce/ForceField.hpp"

AtomType* ForceField::getMatchingAtomType(const string &at) {

  map<string, AtomType*>::iterator iter;
  
  iter = atomTypeMap.find(at);
  if (iter != atomTypeMap.end()) {
    return iter->second; 
  } else {
    return NULL;
  }
}

BondType* ForceField::getMatchingBondType(const string &at1, 
                                          const string &at2) {

  map<pair<string,string>, BondType*>::iterator iter;
  vector<BondType*> foundTypes;

  iter = bondTypeMap.find(pair<at1, at2>);
  if (iter != bondTypeMap.end()) {
    // exact match, so just return it
    return iter->second;
  } 

  iter = bondTypeMap.find(pair<at2, at1>);
  if (iter != bondTypeMap.end()) {
    // exact match in reverse order, so just return it
    return iter->second;
  } 

  iter = bondTypeMap.find(pair<at1, wildCardAtomTypeName>);
  if (iter != bondTypeMap.end()) {
    foundTypes.push_back(iter->second);
  }

  iter = bondTypeMap.find(pair<at2, wildCardAtomTypeName>);
  if (iter != bondTypeMap.end()) {
    foundTypes.push_back(iter->second);
  }

  iter = bondTypeMap.find(pair<wildCardAtomTypeName, at1>);
  if (iter != bondTypeMap.end()) {
    foundTypes.push_back(iter->second);
  }

  iter = bondTypeMap.find(pair<wildCardAtomTypeName, at2>);
  if (iter != bondTypeMap.end()) {
    foundTypes.push_back(iter->second);
  }
  
  if (foundTypes.empty()) {
    return NULL;
  } else {
    

 



  


BendType* ForceField::getMatchingBendType(const string &at1, const string &at2,
                                          const string &at3);
TorsionType* ForceField::getMatchingTorsionType(const string &at1, const string &at2,
                                                const string &at3, const string &at4);

double ForceField::getRcutForAtomType(AtomType* at);


 vector<vector<string> > generateWildcardSequence(const vector<string> atomTypes) { 
   
   vector<vector<string> > results;

   


   vector<vector< string> > getAllWildcardPermutations(const vector<string> myAts) {
     
     int nStrings;
     vector<string> oneResult;
     vector<vector<string> > allResults;

     nStrings = myAts.size();

     if (nStrings == 1) {
       oneResult.push_back(wildcardCharacter);
       allResults.push_back(oneResult);
       return allResults;
     } else {
       
       for (i=0; i < nStrings; i++) {
         oneResult = myAts;
         replace(oneResult.begin(), oneResult.end(), 
