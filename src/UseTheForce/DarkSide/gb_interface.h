/*
 *  gb_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef USETHEFORCE_DARKSIDE_GB_INTERFACE_H
#define USETHEFORCE_DARKSIDE_GB_INTERFACE_H

#define __C

#include "config.h"

#define newGayBerneType F90_FUNC_(set_gb_pair_params, SET_GB_PAIR_PARAMS)

extern "C"{
  void newGayBerneType( double* GB_sigma,
                        double* GB_l2b_ratio,
                        double* GB_eps,
                        double* GB_eps_ratio,
                        double* GB_mu,
                        double* GB_nu );
}
#endif
