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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#ifndef TYPES_ANALYZERSTAMP_HPP
#define TYPES_ANALYZERSTAMP_HPP

#include "types/DataHolder.hpp"

namespace OpenMD {
  class AnalyzerStamp : public DataHolder {
    DeclareParameter(UseAnalyzer, bool);
    DeclareParameter(Method, std::string);
    DeclareParameter(QueryTime, RealType);
    DeclareParameter(NBins, RealType);
    DeclareParameter(NBinsX, RealType);
    DeclareParameter(NBinsY, RealType);
    DeclareParameter(NBinsZ, RealType);
    DeclareParameter(NRBins, RealType);
    DeclareParameter(NAngleBins, RealType);
    DeclareParameter(RCut, RealType);
    DeclareParameter(OOCut, RealType);
    DeclareParameter(ThetaCut, RealType);
    DeclareParameter(OHCut, RealType);
    DeclareParameter(Dz, RealType);
    DeclareParameter(Length, RealType);
    DeclareParameter(ZLength, RealType);
    DeclareParameter(ZOffSet, RealType);
    DeclareParameter(Sele1, std::string);
    DeclareParameter(Sele2, std::string);
    DeclareParameter(Sele3, std::string);
    DeclareParameter(ComSele, std::string);
    DeclareParameter(SeleOffSet, RealType);
    DeclareParameter(SeleOffSet2, RealType);
    DeclareParameter(MolName, std::string);
    DeclareParameter(Begin, RealType);
    DeclareParameter(End, RealType);
    DeclareParameter(Radius, RealType);
    DeclareParameter(VoxelSize, RealType);
    DeclareParameter(GaussWidth, RealType);
    DeclareParameter(PrivilegedAxis, std::string);
    DeclareParameter(PrivilegedAxis2, std::string);
    
  public:
    AnalyzerStamp();
    virtual ~AnalyzerStamp();
    virtual void validate();
  };
}
#endif
