#ifndef __EXCLUDE_H__
#define __EXCLUDE_H__

#include <set>
#include <utility>

using namespace std;

class Exclude{
  
 public:

  ~Exclude();
   
  void addPair(int i, int j);
  int  hasPair(int i, int j);
  void printMe( void );
  int  getSize( void );  
  int* getFortranArray( void );
  static Exclude* Instance();
 
 protected:

  set<pair<int, int> > excludeSet;
  int* exPairs;
  bool newFortranArrayNeeded;
  Exclude();

 private:
  static Exclude* _instance;  

};     

#endif // __EXCLUDE_H__
