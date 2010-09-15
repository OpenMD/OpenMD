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
 
#ifndef NONBONDED_SHAPES_HPP
#define NONBONDED_SHAPES_HPP

#include "types/ShapeAtomType.hpp"
#include "UseTheForce/ForceField.hpp"
#include "math/SquareMatrix3.hpp"

using namespace std;
namespace OpenMD {

  class SHAPES {
    
  public:    
    static SHAPES* Instance();
    static void setForceField(ForceField *ff) {forceField_ = ff;};
    static void initialize();
    static void addType(AtomType* atomType);

    static void calcForce(AtomType* at1, AtomType* at2, const Vector3d d, const RealType rij, const RealType r2, const RealType sw, RealType &vpair, RealType &pot, const RotMat3x3d A1, const RotMat3x3d A2, Vector3d &f1, Vector3d &t1, Vector3d &t2);

    // Fortran support routines;
    static RealType getShapeCut(int atid);
    static void do_shape_pair(int *atid1, int *atid2, RealType *d, RealType *rij, RealType *r2, RealType *sw, RealType *vpair, RealType *pot, RealType *A1, RealType *A2, RealType *f1, RealType *t1, RealType *t2);
    
  private:
    virtual ~SHAPES() { }
    // singleton pattern, prevent reconstruction
    SHAPES() { }
    SHAPES(SHAPES const &) {};
    SHAPES& operator=(SHAPES const&) {};
    static SHAPES* _instance;
  
    static bool initialized_;
    static map<int, AtomType*> ShapesMap;
    static map<pair<AtomType*, AtomType*>, SHAPESInteractionData> MixingMap;
    static ForceField* forceField_;
    static int lMax_;
    static int mMax_;
  };
}

                               
#endif
