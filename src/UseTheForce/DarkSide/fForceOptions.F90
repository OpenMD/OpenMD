!!
!! Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
!!
!! The University of Notre Dame grants you ("Licensee") a
!! non-exclusive, royalty free, license to use, modify and
!! redistribute this software in source and binary code form, provided
!! that the following conditions are met:
!!
!! 1. Acknowledgement of the program authors must be made in any
!!    publication of scientific results based in part on use of the
!!    program.  An acceptable form of acknowledgement is citation of
!!    the article in which the program was described (Matthew
!!    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
!!    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
!!    Parallel Simulation Engine for Molecular Dynamics,"
!!    J. Comput. Chem. 26, pp. 252-271 (2005))
!!
!! 2. Redistributions of source code must retain the above copyright
!!    notice, this list of conditions and the following disclaimer.
!!
!! 3. Redistributions in binary form must reproduce the above copyright
!!    notice, this list of conditions and the following disclaimer in the
!!    documentation and/or other materials provided with the
!!    distribution.
!!
!! This software is provided "AS IS," without a warranty of any
!! kind. All express or implied conditions, representations and
!! warranties, including any implied warranty of merchantability,
!! fitness for a particular purpose or non-infringement, are hereby
!! excluded.  The University of Notre Dame and its licensors shall not
!! be liable for any damages suffered by licensee as a result of
!! using, modifying or distributing the software or its
!! derivatives. In no event will the University of Notre Dame or its
!! licensors be liable for any lost revenue, profit or data, or for
!! direct, indirect, special, consequential, incidental or punitive
!! damages, however caused and regardless of the theory of liability,
!! arising out of the use of or inability to use software, even if the
!! University of Notre Dame has been advised of the possibility of
!! such damages.
!!
!!
!!  fForceOptions.F90
!!  OOPSE-2.0
!!
!!  Created by Charles F. Vardeman II on 12/6/05.
!!
!!  PURPOSE:
!!
!! @author Charles F. Vardeman II 
!! @version $Id: fForceOptions.F90,v 1.1 2005-12-07 19:38:03 chuckv Exp $

!! Handles Mixing options for Fortran.



module  fForceOptions 
   use definitions
   implicit none
   PRIVATE

#define __FORTRAN90
#include "UseTheForce/fForceOptions.h"

type(ForceOptions), save :: fortranForceOptions
logical, save :: haveForceOptions = .false.



public :: getVDW14Scale
public :: getElectrostatic14Scale
public :: getEnergyMixingRule
public :: getDistanceMixingRule
public :: usesGeometricDistanceMixing
public :: usesGeometricEnergyMixing
public :: setForceOptions


contains

subroutine setForceOptions(theseOptions)
type(ForceOptions),intent(in) :: theseOptions
fortranForceOptions = theseOptions
haveForceOptions = .true.
end subroutine setForceOptions


function getVDW14Scale() result(thisScale)
real(kind=dp) :: thisScale
thisScale = fortranForceOptions%vdw14scale
end function getVDW14Scale

function getElectrostatic14Scale() result(thisScale)
real(kind=dp) :: thisScale
thisScale = fortranForceOptions%electrostatic14scale
end function getElectrostatic14Scale


function usesGeometricDistanceMixing() result(doesit)
  logical :: doesit
  doesit = .false.
  if (.not.haveForceOptions) return
  if (fortranForceOptions%DistanceMixingRule == GEOMETRIC_MIXING_RULE) then
     doesit = .true.
  endif
  
end function usesGeometricDistanceMixing


function usesGeometricEnergyMixing() result(doesit)
  logical :: doesit
  doesit = .false.
  if (.not.haveForceOptions) return
  if (fortranForceOptions%EnergyMixingRule == GEOMETRIC_MIXING_RULE) then
     doesit = .true.
  endif
end function usesGeometricEnergyMixing


function getEnergyMixingRule() result(MixingRule)
  integer :: MixingRule
  MixingRule = 0
  if (.not.haveForceOptions) return
  MixingRule = fortranForceOptions%EnergyMixingRule
end function getEnergyMixingRule

function getDistanceMixingRule() result(MixingRule)
  integer :: MixingRule
  MixingRule = 0
  if (.not.haveForceOptions) return
  MixingRule = fortranForceOptions%DistanceMixingRule
end function getDistanceMixingRule



end module fForceOptions
