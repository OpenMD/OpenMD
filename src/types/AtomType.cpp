#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "types/AtomType.hpp"
#include "utils/simError.h"
#define __C
#include "UseTheForce/DarkSide/atype_interface.h"

namespace oopse {
  AtomType::AtomType(){
    
    // initialize to an error:
    atp.ident = -1;

    // and massless:
    mass = 0.0;
    
    // atom type is a Tabula Rasa:
    atp.is_LennardJones = 0;
    atp.is_Electrostatic = 0;
    atp.is_Charge = 0;
    atp.is_Directional = 0;
    atp.is_Sticky = 0;
    atp.is_GayBerne = 0;
    atp.is_EAM = 0;
    atp.is_Shape = 0;
    atp.is_FLARB = 0;  
  }
    
  void AtomType::complete() {
    
    int status;

    if (name.empty()) {
      sprintf( painCave.errMsg,
               "Attempting to complete an AtomType without giving "
               "it a name!\n");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    if (atp.ident == -1) {
      sprintf( painCave.errMsg,
               "Attempting to complete AtomType %s without setting the"
               " ident!/n", name.c_str());
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
 
    status = 0;

    makeAtype(&atp, &status);   
    
    if (status != 0) {
      sprintf( painCave.errMsg,
               "Fortran rejected AtomType %s!\n", name.c_str());
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
  }
}
