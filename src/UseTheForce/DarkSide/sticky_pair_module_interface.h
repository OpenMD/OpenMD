/*
 *  sticky_pair_module_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#define __C
#include "config.h"
extern "C"{
  typedef void (set_sticky_params)( double* sticky_w0, 
                                        double* sticky_v0,
                                        double* sticky_v0p,
                                        double* sticky_rl,
                                        double* sticky_ru,
                                        double* sticky_rlp,
                                        double* sticky_rup ){
     
    F90_FUNC(set_sticky_params, SET_STICKY_PARAMS)(sticky_v0
                                                   sticky_v0p
                                                   sticky_rl
                                                   sticky_ru
                                                   sticky_rlp
                                                   sticky_rup
                                                   );
  }
}