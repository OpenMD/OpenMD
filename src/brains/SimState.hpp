#ifndef __SIMSTATE_H__
#define __SIMSTATE_H__


class SimState{

public:
  SimState();
  ~SimState();

  void createArrays (int the_nElements);
  void destroyArrays(void);

  bool isAllocated(void) { return arraysAllocated; }
  int  getNelements(void){ return nElements; }
    
  void getAtomPointers( int index,
			double** pos, 
			double** vel, 
			double** frc, 
			double** trq, 
			double** Amat,
			double** mu,  
			double** ul);
  

  double* getFrcArray ( void ) { return frc; }
  double* getPosArray ( void ) { return pos; }
  double* getTrqArray ( void ) { return trq; }
  double* getAmatArray( void ) { return Amat; } 
  double* getUlArray  ( void ) { return ul; }
  
private:
  int nElements;        // the number of elements in the arrays
  bool arraysAllocated; // lets us know the arrays have been allocated.

  double* pos;  // the position array
  double* vel;  // the velocity array
  double* frc;  // the forc array
  double* trq;  // the torque vector  ( space fixed )
  double* Amat; // the rotation matrix
  double* mu;   // the dipole moment array
  double* ul;   // the lab frame unit directional vector
};




#endif
