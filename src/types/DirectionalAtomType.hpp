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
    ~DirectionalAtomType();

    /**
     * Finishes off the DirectionalAtomType by communicating the
     * logical portions of the structure to the Fortran atype module
     */
    void    complete();
        
    Mat3x3d getI() {return I;}
    void    setI(Mat3x3d theI) {I = theI;}

    void    setDipole() { atp.is_Dipole = 1; atp.is_Electrostatic = 1; }
    bool    isDipole()  { return atp.is_Dipole; }

    void    setGayBerne() { atp.is_GayBerne = 1; }
    bool    isGayBerne()  { return atp.is_GayBerne; }

    void    setSticky() { atp.is_Sticky = 1; }
    bool    isSticky()  { return atp.is_Sticky; }

    void    setShape() { atp.is_Shape = 1;}
    bool    isShape()  { return atp.is_Shape; }
                
  private:
    
    Mat3x3d I;
    
  };
}
#endif
