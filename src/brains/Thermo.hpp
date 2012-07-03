/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#ifndef BRAINS_THERMO_HPP
#define BRAINS_THERMO_HPP

#include "primitives/Atom.hpp"
#include "brains/SimInfo.hpp"

namespace OpenMD {

  class Thermo{

  public:

    Thermo( SimInfo* info ) : info_(info) {}

    // note: all the following energies are in kcal/mol

    RealType getTranslationalKinetic(); // the translational kinetic energy 
    RealType getRotationalKinetic(); // the rotational kinetic energy 
    RealType getKinetic(); // the total kinetic energy 
    RealType getPotential(); // the total potential energy
    RealType getTotalEnergy(); // gets the total energy

    RealType getTemperature(); // Gives the instant temp. in K
    RealType getElectronicTemperature(); // gives the instant electronic temperature in K

    RealType getPressure(); // gives the instant pressure in atm;

    // gives the pressure tensor in amu*fs^-2*Ang^-1
    Mat3x3d  getPressureTensor(); 
    RealType getVolume();   // gives the volume in Ang^3 

    // accumulate and return the simulation box dipole moment in C*m
    Vector3d getSystemDipole(); 
    Vector3d getHeatFlux();

    /** Returns the center of the mass of the whole system.*/
    Vector3d getCom();

    /** Returns the velocity of center of mass of the whole system.*/
    Vector3d getComVel();

    /** Returns the center of the mass and Center of Mass velocity of
        the whole system.*/ 
    void getComAll(Vector3d& com,Vector3d& comVel);

    /** Returns intertia tensor for the entire system and system
        Angular Momentum.*/
    void getInertiaTensor(Mat3x3d &intertiaTensor,Vector3d &angularMomentum);
    
    /** Returns system angular momentum */
    Vector3d getAngularMomentum();

    /** Returns volume of system as estimated by an ellipsoid defined
        by the radii of gyration */
    RealType getGyrationalVolume();

    /** Overloaded version of gyrational volume that also returns
        det(I) so dV/dr can be calculated */
    void getGyrationalVolume(RealType &vol, RealType &detI);

    RealType getHullVolume();

    RealType getTaggedAtomPairDistance();
    
  private:    
    SimInfo* info_;
  };
  
} //end namespace OpenMD
#endif
