/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#ifndef APPLICATIONS_DYNAMICPROPS_HBONDJUMP_HPP
#define APPLICATIONS_DYNAMICPROPS_HBONDJUMP_HPP

#include "applications/dynamicProps/TimeCorrFunc.hpp"
namespace OpenMD {

  class HBondJump : public TimeCorrFunc<RealType> {
  public:
    HBondJump(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2, double OOCut,
              double thetaCut, double OHCut);
    
  protected:
    virtual void correlation();
    virtual void computeFrame(int istep);
    
    virtual void computeProperty1(int frame) {return;}
    virtual void computeProperty2(int frame) {return;}
    virtual int computeProperty1(int frame, Molecule* mol) {return -1;}
    virtual int computeProperty1(int frame, StuntDouble* sd) {return -1;}
    virtual int computeProperty1(int frame, Bond* bond) {return -1;}
    virtual int computeProperty2(int frame, Molecule* mol) {return -1;}
    virtual int computeProperty2(int frame, StuntDouble* sd) {return -1;}
    virtual int computeProperty2(int frame, Bond* bond) {return -1;}

    virtual RealType calcCorrVal(int frame1, int frame2, int id1, int id2) {return 0.0;}
    virtual RealType calcCorrVal(int frame1, int frame2) { return 0.0; }
       
    virtual void postCorrelate();

    virtual int  registerHydrogen(int frame, int hIndex);
    virtual void findHBonds( int frame );
    bool isHBond(Vector3d donorPos, Vector3d hPos, Vector3d acceptorPos);
    void registerHydrogenBond(int frame, int index, int hIndex, int aIndex);
    void processNonOverlapping(int frame, SelectionManager& sman1, 
                               SelectionManager& sman2);
    void processOverlapping(int frame, SelectionManager& sman);
    
    std::vector<std::vector<int> >  GIDtoH_;
    std::vector<std::vector<int> >  hydrogen_;
    std::vector<std::vector<int> >  acceptor_;
    std::vector<std::vector<int> >  lastAcceptor_;
    std::vector<std::vector<bool> > selected_;
    std::vector<std::vector<int> >  acceptorStartFrame_;
    
    RealType OOCut_;
    RealType thetaCut_;
    RealType OHCut_;

    SelectionManager sele1_minus_common_;
    SelectionManager sele2_minus_common_;
    SelectionManager common_;    
  };

  class HBondJumpZ : public HBondJump {
  public:
    HBondJumpZ(SimInfo* info, const std::string& filename,
               const std::string& sele1, const std::string& sele2, double OOCut,
               double thetaCut, double OHCut, int nZbins, int axis=2);
    virtual int  registerHydrogen(int frame, int hIndex);
    virtual void findHBonds( int frame );
    virtual void correlation();
    virtual void postCorrelate();
    virtual void writeCorrelate();
    
  private:
    std::vector<std::vector<RealType> > histogram_;
    std::vector<std::vector<int> > counts_;
    std::vector<std::vector<int> > zbin_;
    unsigned int nZBins_;
    int axis_;
    std::string axisLabel_;

  };
}
#endif

