#ifndef TYPES_ATOMTYPE_HPP
#define TYPES_ATOMTYPE_HPP

#include "utils/PropertyMap.hpp"
#define __C
#include "types/AtomTypeProperties.h"
#include "UseTheForce/DarkSide/atype_interface.h"

namespace oopse {
  /**
   * @class AtomType
   * AtomType is what OOPSE looks to for unchanging data about an atom.
   * Things that belong to AtomType are universal properties (i.e. does
   * this atom have a Charge?  What is it's mass?)  Dynamic properties of
   * an atom are not intended to be properties of an atom type
   */
  class AtomType{
    
  public:
    
    AtomType();
    ~AtomType();
    /**
     * Finishes off the AtomType by communicating the logical portions of the
     * structure to the Fortran atype module
     */
    void    complete();
    
    void    setMass(double m) { mass = m; }
    double  getMass(void) { return mass; }

    void    setIdent(int id) {atp.ident = id;}
    int     getIdent() {return atp.ident;}

    void    setName(char* n) {name = strdup(n);}
    char*   getName() {return name;}
    
    void    setLennardJones() { atp.is_LennardJones = 1; }
    bool    isLennardJones()  { return atp.is_LennardJones; }
    
    void    setElectrostatic() { atp.is_Electrostatic = 1; }
    bool    isElectrostatic()  { return atp.is_Electrostatic; }
             
    void    setEAM() { atp.is_EAM = 1; }
    bool    isEAM()  { return atp.is_EAM; }
        
    void    setCharge() { atp.is_Charge = 1; atp.is_Electrostatic = 1;}
    bool    isCharge()  { return atp.is_Charge; }

    PropertyMap properties;
                
  protected:
    
    AtomTypeProperties atp;
    double mass;
    char* name;
    
  };
}
#endif
