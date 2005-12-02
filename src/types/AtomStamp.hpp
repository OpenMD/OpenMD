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
 
#ifndef TYPES_ATOMSTAMP_HPP
#define TYPES_ATOMSTAMP_HPP

#include "types/DataHolder.hpp"
namespace oopse {

class AtomStamp  : public DataHolder {
    DeclareParameter(Type, std::string);
    public:
        AtomStamp(int index);
    public:

      bool setPosition(const std::vector<double>& pos);
      bool setOrientation(const std::vector<double>& ort);
      bool havePosition() { return havePos_; }
      bool haveOrientation() { return haveOrt_; }      
      double getPosX() { return position_[0]; }
      double getPosY() { return position_[1]; }
      double getPosZ() { return position_[2]; }
      double getEulerPhi()   { return orientation_[0]; }
      double getEulerTheta() { return orientation_[1]; }
      double getEulerPsi()   { return orientation_[2]; }
      int getIndex() { return index_;}
      virtual void validate();

      AtomStamp* getNextBondedAtom();
      
    private:
        Vector3d position_;
        Vector3d orientation_;
        bool havePos_;
        bool haveOrt_;
        int index_;
        std::vector<int> bonds_;
};

}
#endif
