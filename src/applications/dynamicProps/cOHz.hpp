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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
#ifndef APPLICATIONS_DYNAMICPROPS_COHZ_HPP
#define APPLICATIONS_DYNAMICPROPS_COHZ_HPP

#include "math/Polynomial.hpp"
#include "math/Vector3.hpp"
#include "applications/dynamicProps/ParticleTimeCorrFunc.hpp"

using namespace std;
namespace OpenMD {

  class COHZ : public ParticleTimeCorrFunc {
  public:
    COHZ(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, int order, int nZbins, long long int memSize);   

  private:
    virtual void correlateFrames(int frame1, int frame2);
    virtual RealType calcCorrVal(int frame1, int frame2, StuntDouble* sd1,  StuntDouble* sd2) { return 0.0; }
    virtual Vector3d calcCorrVals(int frame1, int frame2, StuntDouble* sd1, StuntDouble* sd2);
    virtual void writeCorrelate();

    virtual void validateSelection(const SelectionManager& seleMan);


    void setOutputName1(const std::string& filename) {
      outputFilename1_ = filename;
    }
    const std::string& getOutputFileName1() const {
      return outputFilename1_;
    }

    void setOutputName2(const std::string& filename) {
      outputFilename2_ = filename;
    }
    const std::string& getOutputFileName2() const {
      return outputFilename2_;
    }


    int order_;
    DoublePolynomial legendre_;

  protected:
    virtual void preCorrelate();
    virtual void postCorrelate();
    vector<vector<Vector3d> > histogram_;
    vector<vector<int> > counts_;
    int nZBins_;
    std::string outputFilename1_;
    std::string outputFilename2_;

  };

}
#endif
