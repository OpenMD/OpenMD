!!
!! Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
!!
!! The University of Notre Dame grants you ("Licensee") a
!! non-exclusive, royalty free, license to use, modify and
!! redistribute this software in source and binary code form, provided
!! that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright
!!    notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright
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
!! SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
!! research, please cite the appropriate papers when you publish your
!! work.  Good starting points are:
!!                                                                      
!! [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
!! [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
!! [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
!! [4]  Vardeman & Gezelter, in progress (2009).
!!


!! Calculates Long Range forces Lennard-Jones interactions.
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @version $Id$, $Date$, $Name: not supported by cvs2svn $, $Revision$


module lj
  use definitions
  use atype_module
  use vector_class
  use simulation
  use status
  use fForceOptions
  use force_globals

  implicit none
  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"

  logical, save :: useGeometricDistanceMixing = .false.
  logical, save :: haveMixingMap = .false.

  real(kind=DP), save :: defaultCutoff = 0.0_DP
  logical, save :: defaultShiftPot = .false.
  logical, save :: defaultShiftFrc = .false.
  logical, save :: haveDefaultCutoff = .false.

  type, private :: LJtype
     integer       :: atid
     real(kind=dp) :: sigma
     real(kind=dp) :: epsilon
     logical       :: isSoftCore = .false.
  end type LJtype

  type, private :: LJList
     integer               :: Nljtypes = 0
     integer               :: currentLJtype = 0
     type(LJtype), pointer :: LJtypes(:)      => null()
     integer, pointer      :: atidToLJtype(:) => null()
  end type LJList

  type(LJList), save :: LJMap

  type :: MixParameters
     real(kind=DP) :: sigma
     real(kind=DP) :: epsilon
     real(kind=dp) :: sigmai
     real(kind=dp) :: rCut
     logical       :: rCutWasSet = .false.
     logical       :: shiftedPot
     logical       :: isSoftCore = .false.
     logical       :: shiftedFrc
  end type MixParameters

  type(MixParameters), dimension(:,:), allocatable :: MixingMap

  public :: newLJtype
  public :: setLJDefaultCutoff
  public :: getSigma
  public :: getEpsilon
  public :: do_lj_pair
  public :: destroyLJtypes

contains

  subroutine newLJtype(c_ident, sigma, epsilon, isSoftCore, status)
    integer,intent(in) :: c_ident
    real(kind=dp),intent(in) :: sigma
    real(kind=dp),intent(in) :: epsilon
    integer, intent(in) :: isSoftCore
    integer,intent(out) :: status
    integer :: nLJTypes, ntypes, myATID
    integer, pointer :: MatchList(:) => null()
    integer :: current

    status = 0
    ! check to see if this is the first time into this routine...
    if (.not.associated(LJMap%LJtypes)) then

       call getMatchingElementList(atypes, "is_LennardJones", .true., &
            nLJTypes, MatchList)
       
       LJMap%nLJtypes =  nLJTypes

       allocate(LJMap%LJtypes(nLJTypes))

       ntypes = getSize(atypes)

       allocate(LJMap%atidToLJtype(ntypes))
    end if

    LJMap%currentLJtype = LJMap%currentLJtype + 1
    current = LJMap%currentLJtype

    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    LJMap%atidToLJtype(myATID)        = current
    LJMap%LJtypes(current)%atid       = myATID
    LJMap%LJtypes(current)%sigma      = sigma
    LJMap%LJtypes(current)%epsilon    = epsilon 
    if (isSoftCore .eq. 1) then
       LJMap%LJtypes(current)%isSoftCore = .true.
    else
       LJMap%LJtypes(current)%isSoftCore = .false.
    endif
  end subroutine newLJtype

  subroutine setLJDefaultCutoff(thisRcut, shiftedPot, shiftedFrc)
    real(kind=dp), intent(in) :: thisRcut
    logical, intent(in) :: shiftedPot
    logical, intent(in) :: shiftedFrc
    defaultCutoff = thisRcut
    defaultShiftPot = shiftedPot
    defaultShiftFrc = shiftedFrc
    haveDefaultCutoff = .true.

    !we only want to build LJ Mixing map if LJ is being used.
    if(LJMap%nLJTypes /= 0) then
       call createMixingMap()
    end if

  end subroutine setLJDefaultCutoff

  function getSigma(atid) result (s)
    integer, intent(in) :: atid
    integer :: ljt1
    real(kind=dp) :: s

    if (LJMap%currentLJtype == 0) then
       call handleError("LJ", "No members in LJMap")
       return
    end if

    ljt1 = LJMap%atidToLJtype(atid)
    s = LJMap%LJtypes(ljt1)%sigma

  end function getSigma

  function getEpsilon(atid) result (e)
    integer, intent(in) :: atid
    integer :: ljt1
    real(kind=dp) :: e

    if (LJMap%currentLJtype == 0) then
       call handleError("LJ", "No members in LJMap")
       return
    end if

    ljt1 = LJMap%atidToLJtype(atid)
    e = LJMap%LJtypes(ljt1)%epsilon

  end function getEpsilon

  subroutine createMixingMap()
    integer :: nLJtypes, i, j
    real ( kind = dp ) :: s1, s2, e1, e2
    real ( kind = dp ) :: rcut6, tp6, tp12
    logical :: isSoftCore1, isSoftCore2, doShift

    if (LJMap%currentLJtype == 0) then
       call handleError("LJ", "No members in LJMap")
       return
    end if

    nLJtypes = LJMap%nLJtypes

    if (.not. allocated(MixingMap)) then
       allocate(MixingMap(nLJtypes, nLJtypes))
    endif

    useGeometricDistanceMixing = usesGeometricDistanceMixing()
    do i = 1, nLJtypes

       s1 = LJMap%LJtypes(i)%sigma
       e1 = LJMap%LJtypes(i)%epsilon
       isSoftCore1 = LJMap%LJtypes(i)%isSoftCore

       do j = i, nLJtypes
          
          s2 = LJMap%LJtypes(j)%sigma
          e2 = LJMap%LJtypes(j)%epsilon
          isSoftCore2 = LJMap%LJtypes(j)%isSoftCore
          
          MixingMap(i,j)%isSoftCore = isSoftCore1 .or. isSoftCore2

          ! only the distance parameter uses different mixing policies
          if (useGeometricDistanceMixing) then
             MixingMap(i,j)%sigma = sqrt(s1 * s2)
          else
             MixingMap(i,j)%sigma = 0.5_dp * (s1 + s2)
          endif
          
          MixingMap(i,j)%epsilon = sqrt(e1 * e2)

          MixingMap(i,j)%sigmai = 1.0_DP  / (MixingMap(i,j)%sigma)

          if (haveDefaultCutoff) then
             MixingMap(i,j)%shiftedPot = defaultShiftPot
             MixingMap(i,j)%shiftedFrc = defaultShiftFrc
          else
             MixingMap(i,j)%shiftedPot = defaultShiftPot
             MixingMap(i,j)%shiftedFrc = defaultShiftFrc
          endif          

          if (i.ne.j) then
             MixingMap(j,i)%sigma      = MixingMap(i,j)%sigma
             MixingMap(j,i)%epsilon    = MixingMap(i,j)%epsilon
             MixingMap(j,i)%sigmai     = MixingMap(i,j)%sigmai
             MixingMap(j,i)%rCut       = MixingMap(i,j)%rCut
             MixingMap(j,i)%rCutWasSet = MixingMap(i,j)%rCutWasSet
             MixingMap(j,i)%shiftedPot = MixingMap(i,j)%shiftedPot
             MixingMap(j,i)%isSoftCore = MixingMap(i,j)%isSoftCore
             MixingMap(j,i)%shiftedFrc = MixingMap(i,j)%shiftedFrc
          endif

       enddo
    enddo
    
    haveMixingMap = .true.

  end subroutine createMixingMap
          
  subroutine do_lj_pair(atid1, atid2, d, rij, r2, rcut, sw, vdwMult, &
       vpair, fpair, pot, f1)
    
    integer, intent(in) ::  atid1, atid2
    integer :: ljt1, ljt2
    real( kind = dp ), intent(in) :: rij, r2, rcut, vdwMult
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3) :: f1
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair


    ! local Variables
    real( kind = dp ) :: drdx, drdy, drdz
    real( kind = dp ) :: fx, fy, fz
    real( kind = dp ) :: myPot, myPotC, myDeriv, myDerivC, ros, rcos
    real( kind = dp ) :: pot_temp, dudr
    real( kind = dp ) :: sigmai
    real( kind = dp ) :: epsilon
    logical :: isSoftCore, shiftedPot, shiftedFrc
    integer :: id1, id2, localError

    if (.not.haveMixingMap) then
       call createMixingMap()
    endif

    ljt1 = LJMap%atidToLJtype(atid1)
    ljt2 = LJMap%atidToLJtype(atid2)

    sigmai     = MixingMap(ljt1,ljt2)%sigmai
    epsilon    = MixingMap(ljt1,ljt2)%epsilon
    isSoftCore = MixingMap(ljt1,ljt2)%isSoftCore
    shiftedPot = MixingMap(ljt1,ljt2)%shiftedPot
    shiftedFrc = MixingMap(ljt1,ljt2)%shiftedFrc

    ros = rij * sigmai

    if (isSoftCore) then

       call getSoftFunc(ros, myPot, myDeriv)

       if (shiftedPot) then
          rcos = rcut * sigmai
          call getSoftFunc(rcos, myPotC, myDerivC)
          myDerivC = 0.0_dp
       elseif (shiftedFrc) then
          rcos = rcut * sigmai
          call getSoftFunc(rcos, myPotC, myDerivC)
          myPotC = myPotC + myDerivC * (rij - rcut) * sigmai
       else
          myPotC = 0.0_dp
          myDerivC = 0.0_dp
       endif
              
    else

       call getLJfunc(ros, myPot, myDeriv)

       if (shiftedPot) then
          rcos = rcut * sigmai
          call getLJfunc(rcos, myPotC, myDerivC) 
          myDerivC = 0.0_dp
       elseif (shiftedFrc) then
          rcos = rcut * sigmai
          call getLJfunc(rcos, myPotC, myDerivC)
          myPotC = myPotC + myDerivC * (rij - rcut) * sigmai
       else
          myPotC = 0.0_dp
          myDerivC = 0.0_dp
       endif
       
    endif

    pot_temp = vdwMult * epsilon * (myPot - myPotC)
    vpair = vpair + pot_temp
    dudr = sw * vdwMult * epsilon * (myDeriv - myDerivC) * sigmai

    drdx = d(1) / rij
    drdy = d(2) / rij
    drdz = d(3) / rij

    fx = dudr * drdx
    fy = dudr * drdy
    fz = dudr * drdz

    pot = pot + sw*pot_temp

    f1(1) = f1(1) + fx
    f1(2) = f1(2) + fy
    f1(3) = f1(3) + fz

    return    

  end subroutine do_lj_pair

  subroutine destroyLJTypes()

    LJMap%nLJtypes = 0
    LJMap%currentLJtype = 0
    
    if (associated(LJMap%LJtypes)) then
       deallocate(LJMap%LJtypes)
       LJMap%LJtypes => null()
    end if
    
    if (associated(LJMap%atidToLJtype)) then
       deallocate(LJMap%atidToLJtype)
       LJMap%atidToLJtype => null()
    end if
    
    haveMixingMap = .false.

  end subroutine destroyLJTypes

  subroutine getLJfunc(r, myPot, myDeriv)

    real(kind=dp), intent(in) :: r
    real(kind=dp), intent(inout) :: myPot, myDeriv
    real(kind=dp) :: ri, ri2, ri6, ri7, ri12, ri13
    real(kind=dp) :: a, b, c, d, dx
    integer :: j

    ri = 1.0_DP / r
    ri2 = ri*ri
    ri6 = ri2*ri2*ri2
    ri7 = ri6*ri
    ri12 = ri6*ri6
    ri13 = ri12*ri
    
    myPot = 4.0_DP * (ri12 - ri6)
    myDeriv = 24.0_DP * (ri7 - 2.0_DP * ri13)
    
    return
  end subroutine getLJfunc

  subroutine getSoftFunc(r, myPot, myDeriv)
    
    real(kind=dp), intent(in) :: r
    real(kind=dp), intent(inout) :: myPot, myDeriv
    real(kind=dp) :: ri, ri2, ri6, ri7
    real(kind=dp) :: a, b, c, d, dx
    integer :: j
    
    ri = 1.0_DP / r    
    ri2 = ri*ri
    ri6 = ri2*ri2*ri2
    ri7 = ri6*ri
    myPot = 4.0_DP * (ri6)
    myDeriv = - 24.0_DP * ri7 
    
    return
  end subroutine getSoftFunc

end module lj
