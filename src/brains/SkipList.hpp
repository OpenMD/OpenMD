#ifndef __SKIPLIST_H__
#define __SKIPLIST_H__

#include <set>
#include <utility>

using namespace std;

class SkipList{
  
 public:

  ~SkipList();
   
  void addAtom(int i);
  int  hasAtom(int i);
  void printMe( void );
  int  getSize( void );  
  static SkipList* Instance();
 
 protected:

  set<int> skipSet;
  SkipList();

 private:
  static SkipList* _instance;  

};     

#endif // __SKIPLIST_H__
