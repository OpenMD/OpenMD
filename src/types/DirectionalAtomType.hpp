/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#ifndef TYPES_DIRECTIONALATOMTYPE_HPP
#define TYPES_DIRECTIONALATOMTYPE_HPP

#include "types/AtomType.hpp"
#include "math/SquareMatrix3.hpp"

namespace oopse {

/**
 * @class DirectionalAtomType 
 *
 * DirectionalAtomType is what OOPSE looks to for unchanging data
 * about a directional atoms. 
 */
class DirectionalAtomType : public AtomType {

    public:

        DirectionalAtomType() : AtomType() { atp.is_Directional = 1; }

        Mat3x3d getI() {return I;}

        void    setI(Mat3x3d theI) {I = theI;}

        RotMat3x3d getElectroBodyFrame() {
            return electroBodyFrame_;
        }

        void setElectroBodyFrame(const RotMat3x3d& electroBodyFrame) {
            electroBodyFrame_ =electroBodyFrame;
        }

        void setDipole() { atp.is_Dipole = 1; }
        void setSplitDipole() { atp.is_SplitDipole = 1; atp.is_Dipole=1;}
        void setQuadrupole() { atp.is_Quadrupole = 1; }
        void setGayBerne() { atp.is_GayBerne = 1; }
        void setSticky() { atp.is_Sticky = 1; }
        void setShape() { atp.is_Shape = 1;}

        virtual void complete();

    private:

        Mat3x3d I;
        RotMat3x3d electroBodyFrame_;
};


struct StickyParam {
double w0;
double v0;
double v0p;
double rl;
double ru;
double rlp;
double rup;
};

typedef SimpleTypeData<StickyParam> StickyParamGenericData;

typedef SimpleTypeData<Vector3d> Vector3dGenericData;
  
}
#endif
