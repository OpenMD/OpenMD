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
 
#ifndef TYPES_EAMADAPTER_HPP
#define TYPES_EAMADAPTER_HPP

#include <memory>

#include "utils/GenericData.hpp"
#include "types/AtomType.hpp"
#include "math/CubicSpline.hpp"

namespace OpenMD {

  string const EAMtypeID = "EAM";
  string const FuncflTypeID = "FUNCFL";
  string const ZhouTypeID = "ZHOU";
  string const ZhouRoseTypeID = "ZHOUROSE";
  

  enum EAMType{
    eamFuncfl,
    eamZhou2001,
    eamZhou2004,
    eamZhou2005,
    eamZhou2005Oxygen,
    eamZhouRose,
    eamOxygenFuncfl,
    eamUnknown
  };
  
  struct EAMParameters {
    EAMType eamType;
    std::string latticeType;
    RealType latticeConstant;  
  };
  
  struct FuncflParameters {
    // This first set is for parameters read from DYNAMO 86 funcfl files:
    int atomicNumber;
    RealType atomicMass;
    int nrho;
    RealType drho;
    int nr;
    RealType dr;
    RealType rcut;
    std::vector<RealType> Z;   // Z(r) 
    std::vector<RealType> rho; // rho(r)
    std::vector<RealType> F;   // F[rho]
  };

  struct ZhouParameters {
    // This set is for parameters for the parameterization of EAM
    // described in: Acta mater 49, 4005 (2001)
    RealType re;
    RealType fe;
    RealType rhoe;
    RealType alpha;
    RealType beta;
    RealType A;
    RealType B;
    RealType kappa;
    RealType lambda;
    std::vector<RealType> F;
    std::vector<RealType> Fn;
    RealType eta;
    RealType Fe;
    // These additional parameters are for the parameterization of EAM
    // described in: in: X. W. Zhou, R. A. Johnson, and
    // H. N. G. Wadley, Phys. Rev. B, 69, 144113 (2004).
    RealType rhos;
    RealType rhol;
    RealType rhoh;
    // These additional parameters are for the parammeterization of
    // EAM described in:  X. W. Zhou and H. N. G. Wadley,
    // J. Phys. Condens. Matter 17 (2005) 3619-3635.
    RealType F3minus;
    RealType F3plus;
    RealType gamma;
    RealType nu;
    std::vector<RealType> OrhoLimits;
    std::vector<RealType> OrhoE;
    std::vector<std::vector<RealType> > OF;
    // For a Rose-type functional
    RealType F0;
  };
  
  
  typedef SimpleTypeData<EAMParameters> EAMData;
  typedef SimpleTypeData<FuncflParameters> FuncflData;
  typedef SimpleTypeData<ZhouParameters> ZhouData;

  class EAMAdapter {
  public:
    EAMAdapter(AtomType* AT) { at_ = AT; };

    void makeFuncfl(RealType latticeConstant,
                    std::string latticeType,
                    int nrho,
                    RealType drho,
                    int nr,
                    RealType dr,
                    RealType rcut,
                    std::vector<RealType> Z,
                    std::vector<RealType> rho,
                    std::vector<RealType> F);

    void makeZhou2001(std::string latticeType,
                      RealType re,
                      RealType fe,
                      RealType rhoe,
                      RealType alpha,
                      RealType beta,
                      RealType A,
                      RealType B,
                      RealType kappa,
                      RealType lambda,
                      std::vector<RealType> Fn,
                      std::vector<RealType> F,
                      RealType eta,
                      RealType Fe);
    
    void makeZhou2004(std::string latticeType,
                      RealType re,
                      RealType fe,
                      RealType rhoe,
                      RealType rhos,
                      RealType alpha,
                      RealType beta,
                      RealType A,
                      RealType B,
                      RealType kappa,
                      RealType lambda,
                      std::vector<RealType> Fn,
                      std::vector<RealType> F,
                      RealType eta,
                      RealType Fe,
                      RealType rhol,
                      RealType rhoh);
    
    void makeZhou2005(std::string latticeType,
                      RealType re,
                      RealType fe,
                      RealType rhoe,
                      RealType rhos,
                      RealType alpha,
                      RealType beta,
                      RealType A,
                      RealType B,
                      RealType kappa,
                      RealType lambda,
                      std::vector<RealType> Fn,
                      std::vector<RealType> F,
                      RealType F3minus,
                      RealType F3plus,
                      RealType eta,
                      RealType Fe);
    
    void makeZhou2005Oxygen(RealType re,
                            RealType fe,
                            RealType alpha,
                            RealType beta,
                            RealType A,
                            RealType B,
                            RealType kappa,
                            RealType lambda,
                            RealType gamma,
                            RealType nu,
                            std::vector<RealType> OrhoLimits,
                            std::vector<RealType> OrhoE,
                            std::vector<std::vector<RealType> > OF);
    
    void makeZhouRose(RealType re,
                      RealType fe,
                      RealType rhoe,
                      RealType alpha,
                      RealType beta,
                      RealType A,
                      RealType B,
                      RealType kappa,
                      RealType lambda,
                      RealType F0);
    void makeOxygenFuncfl(RealType re,
                      RealType fe,
                      RealType alpha,
                      RealType beta,
                      RealType A,
                      RealType B,
                      RealType kappa,
                      RealType lambda,
                      RealType drho,
                      RealType nrho,
                      std::vector<RealType> F);

    
    bool isEAM();
    bool hasSplines();
    EAMType getEAMType();
    std::string getLatticeType();
    RealType getLatticeConstant();
    int getNr();
    RealType getDr();
    int getNrho();
    RealType getDrho();
    RealType getRcut();
    
    RealType getRe();
    RealType get_fe();
    RealType getRhoe();
    RealType getAlpha();
    RealType getBeta();
    RealType getA();
    RealType getB();
    RealType getKappa();
    RealType getLambda();
    RealType getGamma();
    RealType getNu();
    std::vector<RealType> getFn();
    std::vector<RealType> getF();
    RealType getF3plus();
    RealType getF3minus();
    RealType getEta();
    RealType getFe();
    RealType getRhos();
    RealType getRhol();
    RealType getRhoh();
    std::vector<RealType> getOrhoLimits();
    std::vector<RealType> getOrhoE();
    std::vector<std::vector<RealType> > getOF();
    RealType getF0();
    CubicSplinePtr getZSpline();
    CubicSplinePtr getRhoSpline();
    CubicSplinePtr getFSpline();

  private:
    AtomType* at_;
    EAMParameters      getEAMParam();
    FuncflParameters   getFuncflParam();
    ZhouParameters     getZhouParam();
  };
  
}
#endif
