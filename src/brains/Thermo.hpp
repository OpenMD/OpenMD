#ifndef __THERMO_H__
#define __THERMO_H__

#include "Atom.hpp"
#include "SRI.hpp"
#include "SimInfo.hpp"
#include "randomSPRNG.hpp"

class Thermo{

public:
  
  Thermo( SimInfo* the_info );
  ~Thermo();

  // note: all the following energies are in kcal/mol

  double getKinetic(); // the total kinetic energy 
  double getPotential(); // the total potential energy
  double getTotalE(); // gets the total energy

  double getTemperature(); // gives the instant temp. in K

  double getPressure(); // gives the instant pressure in atm;
  double getPressureX(); // gives the instant pressure in atm;
  double getPressureY(); // gives the instant pressure in atm;
  double getPressureZ(); // gives the instant pressure in atm;

  void   getPressureTensor(double press[3][3]); // gives the pressure 
                                                // tensor in 
                                                // amu*fs^-2*Ang^-1
  double getVolume();   // gives the volume in Ang^3 

  int getNDF();    // get the number of degrees of freedom in the system
  int getNDFraw(); // get the number of raw degrees of freedom in the system
                   // i.e. don't subtract constraints or system COM.
  
  void velocitize(); // set the temperature to the target temp in SimInfo
                     // NOTE: srand48 should be seeded before calling.
  void getCOMVel(double vdrift[3]);
  void getCOM(double COM[3]);
  void removeCOMdrift();

private:
  SimInfo* info;
  gaussianSPRNG *gaussStream;
};
#endif
