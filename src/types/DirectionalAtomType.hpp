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

    void    setDipole() { atp.is_Dipole = 1; atp.is_Electrostatic = 1; }

    void    setGayBerne() { atp.is_GayBerne = 1; }

    void    setSticky() { atp.is_Sticky = 1; }

    void    setShape() { atp.is_Shape = 1;}

                
  private:
    
    Mat3x3d I;
    
  };
}
#endif
