/*
 *  sticky_pair_module_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 University of Notre Dame. All rights reserved.
 *
 */
#ifndef USETHEFORCE_DARKSIDE_STICKY_INTERFACE_H
#define USETHEFORCE_DARKSIDE_STICKY_INTERFACE_H

#define __C
#include "config.h"

extern "C" {
  
  void F90_FUNC_(set_sticky_params, SET_STICKY_PARAMS)( double* sticky_w0, 
                                                        double* sticky_v0,
                                                        double* sticky_v0p,
                                                        double* sticky_rl,
                                                        double* sticky_ru,
                                                        double* sticky_rlp,
                                                        double* sticky_rup );
  
}

#endif
