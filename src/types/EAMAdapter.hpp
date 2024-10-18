/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#ifndef TYPES_EAMADAPTER_HPP
#define TYPES_EAMADAPTER_HPP

#include <memory>

#include "math/CubicSpline.hpp"
#include "types/AtomType.hpp"
#include "utils/GenericData.hpp"

namespace OpenMD {

  string const EAMtypeID      = "EAM";
  string const FuncflTypeID   = "FUNCFL";
  string const ZhouTypeID     = "ZHOU";
  string const ZhouRoseTypeID = "ZHOUROSE";

  enum EAMType {
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
    std::vector<RealType> Z;    // Z(r)
    std::vector<RealType> rho;  // rho(r)
    std::vector<RealType> F;    // F[rho]
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
    std::vector<std::vector<RealType>> OF;
    // For a Rose-type functional
    RealType F0;
  };

  using EAMData    = SimpleTypeData<EAMParameters>;
  using FuncflData = SimpleTypeData<FuncflParameters>;
  using ZhouData   = SimpleTypeData<ZhouParameters>;

  class EAMAdapter {
  public:
    EAMAdapter(AtomType* AT) { at_ = AT; };

    void makeFuncfl(RealType latticeConstant, std::string latticeType, int nrho,
                    RealType drho, int nr, RealType dr, RealType rcut,
                    std::vector<RealType> Z, std::vector<RealType> rho,
                    std::vector<RealType> F);

    void makeZhou2001(std::string latticeType, RealType re, RealType fe,
                      RealType rhoe, RealType alpha, RealType beta, RealType A,
                      RealType B, RealType kappa, RealType lambda,
                      std::vector<RealType> Fn, std::vector<RealType> F,
                      RealType eta, RealType Fe);

    void makeZhou2004(std::string latticeType, RealType re, RealType fe,
                      RealType rhoe, RealType rhos, RealType alpha,
                      RealType beta, RealType A, RealType B, RealType kappa,
                      RealType lambda, std::vector<RealType> Fn,
                      std::vector<RealType> F, RealType eta, RealType Fe,
                      RealType rhol, RealType rhoh);

    void makeZhou2005(std::string latticeType, RealType re, RealType fe,
                      RealType rhoe, RealType rhos, RealType alpha,
                      RealType beta, RealType A, RealType B, RealType kappa,
                      RealType lambda, std::vector<RealType> Fn,
                      std::vector<RealType> F, RealType F3minus,
                      RealType F3plus, RealType eta, RealType Fe);

    void makeZhou2005Oxygen(RealType re, RealType fe, RealType alpha,
                            RealType beta, RealType A, RealType B,
                            RealType kappa, RealType lambda, RealType gamma,
                            RealType nu, std::vector<RealType> OrhoLimits,
                            std::vector<RealType> OrhoE,
                            std::vector<std::vector<RealType>> OF);

    void makeZhouRose(RealType re, RealType fe, RealType rhoe, RealType alpha,
                      RealType beta, RealType A, RealType B, RealType kappa,
                      RealType lambda, RealType F0);
    void makeOxygenFuncfl(RealType re, RealType fe, RealType alpha,
                          RealType beta, RealType A, RealType B, RealType kappa,
                          RealType lambda, RealType drho, RealType nrho,
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
    std::vector<std::vector<RealType>> getOF();
    RealType getF0();
    CubicSplinePtr getZSpline();
    CubicSplinePtr getRhoSpline();
    CubicSplinePtr getFSpline();

  private:
    AtomType* at_;
    EAMParameters getEAMParam();
    FuncflParameters getFuncflParam();
    ZhouParameters getZhouParam();
  };
}  // namespace OpenMD

#endif
