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
 * * 2. Redistributions in binary form must reproduce the above copyright
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
#ifndef ANALYSIS_INTERACTIONTYPE_HPP
#define ANALYSIS_INTERACTIONTYPE_HPP

  /**
   * @class InteractionType
   * @brief Interaction Type
   *
   * This class is the interface for other classes to interact with Analyzers.
   * This class defines the types (ie, atom, stuntDouble, pairs, triples)
   * of analysis.
   *
   */

namespace OpenMD {

  class InteractionType {

  public:
    InteractionType();
    virtual ~InteractionType();

    void setEvaluator()=0;
    void setSeleMan()=0;
    void setSeleScript()=0;
    

  protected:

  };
  

  class SingletType : public InteractionType {
    
  public:
    SingletType();
    virtual ~SingletType();

    void setEvaluator(SimInfo* info) { evaluator_(info); }
    void setSeleMan(SimInfo* info) { seleMan_(info); }
    void setSeleScript(const std::string& sele1) {selectionScript_(sele1); }

  protected:
    std::string selectionScript_;
    SelectionManager seleMan_;
    
  };
  

  class DoubletType : public InteractionType {

  public:
    DoubletType();
    virtual ~DoubletType();

    void setEvaluator(SimInfo* info) { evaluator1_(info); evaluator2_(info); }
    void setSeleMan(SimInfo* info) { seleMan1_(info); seleMan2_(info); }
    void setSeleScript(const std::string& sele1, const std::string& sele2) {
      selectionScript1_(sele1); selectionScript2_(sele2); }
    
  protected:
    std::string selectionScript1_;
    SelectionManager seleMan1_;    
    SelectionEvaluator evaluator1_;
    std::string selectionScript2_;
    SelectionManager seleMan2_;    
    SelectionEvaluator evaluator2_;
    
  };
  
  // PairType will handle analysis dependent on pairs of things
  // e.g., GofR
  class PairType : public InteractionType {
    
  public:
    PairType();
    virtual ~PairType();

    void setEvaluator(SimInfo* info) { evaluator1_(info); evaluator2_(info); }
    void setSeleMan(SimInfo* info) { seleMan1_(info); seleMan2_(info); }
    
  protected:
    
    std::string selectionScript1_;
    SelectionManager seleMan1_;    
    SelectionEvaluator evaluator1_;
    std::string selectionScript2_;
    SelectionManager seleMan2_;    
    SelectionEvaluator evaluator2_;
    
  };
  

  class TripletType : public InteractionType {
    
  public:
    TripletType();
    virtual ~TripletType();
    
    void setEvaluator(SimInfo* info) { evaluator1_(info);
      evaluator2_(info); evaluator3_(info); }
    
    void setSeleMan(SimInfo* info) { seleMan1_(info);
      seleMan2_(info); seleMan3_(info); }
    
    void setSeleScript(const std::string& sele1, const std::string& sele2,
		       const std::string& sele3) {
      selectionScript1_(sele1); selectionScript2_(sele2);
      selectionScript3_(sele3); }


    
  protected:
    std::string selectionScript1_;
    SelectionManager seleMan1_;    
    SelectionEvaluator evaluator1_;
    std::string selectionScript2_;
    SelectionManager seleMan2_;    
    SelectionEvaluator evaluator2_;
    std::string selectionScript3_;
    SelectionManager seleMan3_;    
    SelectionEvaluator evaluator3_;


  };
  

}
#endif
