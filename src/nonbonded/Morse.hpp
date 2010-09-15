/*
 * Copyright (c) 2009 The University of Notre Dame. All Rights Reserved.
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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
 
#ifndef NONBONDED_MORSE_HPP
#define NONBONDED_MORSE_HPP

#include "types/AtomType.hpp"
#include "UseTheForce/ForceField.hpp"
#include "math/Vector3.hpp"

using namespace std;
namespace OpenMD {

  enum MorseInteractionType {
    shiftedMorse,
    repulsiveMorse,
  };

  struct MorseInteractionData {
    RealType De;
    RealType Re;
    RealType beta;
    MorseInteractionType interactionType;
  };

  class Morse {
    
  public:    
    static Morse* Instance();
    static void setForceField(ForceField *ff) {forceField_ = ff;};
    static void initialize();
    static void addExplicitInteraction(AtomType* atype1, AtomType* atype2, RealType De, RealType Re, RealType beta, MorseInteractionType mit);
    static void calcForce(AtomType* at1, AtomType* at2, const Vector3d d, const RealType rij, const RealType r2, const RealType rcut, const RealType sw, const RealType vdwMult, RealType &vpair, RealType &pot, Vector3d &f1);
    
    // Fortran support routines;
    static void do_morse_pair(int *atid1, int *atid2, RealType *d, RealType *rij, RealType *r2, RealType *rcut, RealType *sw, RealType *vdwMult, RealType *vpair, RealType *pot, RealType *f1);
    static void setMorseDefaultCutoff(RealType *thisRcut, int *sP, int *sF);

  private:
    virtual ~Morse() { }
    // singleton pattern, prevent reconstruction
    Morse() { }
    Morse(Morse const &) {};
    Morse& operator=(Morse const&) {};
    static Morse* _instance;      

    static bool initialized_;
    // MorseMap will be used for providing access from Fortran.
    // All of the C++ native classes should use AtomType*, but the
    // fortran calls use int, so we need to look these up first before
    // we can do anything else;
    static map<int, AtomType*> MorseMap;
    static map<pair<AtomType*, AtomType*>, MorseInteractionData> MixingMap;
    static bool shiftedPot_;
    static bool shiftedFrc_;
    static ForceField* forceField_;    
    static map<string, MorseInteractionType> stringToEnumMap_;

  };
}

                               
#endif
