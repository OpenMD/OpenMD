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


!! Calculates Long Range forces Lennard-Jones interactions.
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @version $Id: LJ.F90,v 1.19 2005-11-21 22:59:01 gezelter Exp $, $Date: 2005-11-21 22:59:01 $, $Name: not supported by cvs2svn $, $Revision: 1.19 $


module lj
  use atype_module
  use vector_class
  use simulation
  use status
#ifdef IS_MPI
  use mpiSimulation
#endif
  use force_globals

  implicit none
  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"

  integer, parameter :: DP = selected_real_kind(15)

  logical, save :: useGeometricDistanceMixing = .false.
  logical, save :: haveMixingMap = .false.

  real(kind=DP), save :: defaultCutoff = 0.0_DP
  logical, save :: defaultShift = .false.
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
     real(kind=dp) :: sigma6
     real(kind=dp) :: rCut
     real(kind=dp) :: delta
     logical       :: rCutWasSet = .false.
     logical       :: shiftedPot
     logical       :: isSoftCore = .false.
  end type MixParameters

  type(MixParameters), dimension(:,:), allocatable :: MixingMap

  public :: newLJtype
  public :: setLJDefaultCutoff
  public :: getSigma
  public :: getEpsilon
  public :: useGeometricMixing
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

  subroutine setLJDefaultCutoff(thisRcut, shiftedPot)
    real(kind=dp), intent(in) :: thisRcut
    logical, intent(in) :: shiftedPot
    defaultCutoff = thisRcut
    defaultShift = shiftedPot
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

  subroutine useGeometricMixing() 
    useGeometricDistanceMixing = .true.
    haveMixingMap = .false.
    return
  end subroutine useGeometricMixing

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
             MixingMap(i,j)%sigma = dsqrt(s1 * s2)
          else
             MixingMap(i,j)%sigma = 0.5_dp * (s1 + s2)
          endif
          
          MixingMap(i,j)%epsilon = dsqrt(e1 * e2)

          MixingMap(i,j)%sigma6 = (MixingMap(i,j)%sigma)**6

          if (haveDefaultCutoff) then
             rcut6 = defaultCutoff**6
             tp6    = MixingMap(i,j)%sigma6/rcut6
             tp12    = tp6**2          
             MixingMap(i,j)%delta =-4.0_DP*MixingMap(i,j)%epsilon*(tp12 - tp6)
             MixingMap(i,j)%shiftedPot = defaultShift
          else
             MixingMap(i,j)%delta = 0.0_DP
             MixingMap(i,j)%shiftedPot = defaultShift
          endif          

          if (i.ne.j) then
             MixingMap(j,i)%sigma      = MixingMap(i,j)%sigma
             MixingMap(j,i)%epsilon    = MixingMap(i,j)%epsilon
             MixingMap(j,i)%sigma6     = MixingMap(i,j)%sigma6
             MixingMap(j,i)%rCut       = MixingMap(i,j)%rCut
             MixingMap(j,i)%delta      = MixingMap(i,j)%delta
             MixingMap(j,i)%rCutWasSet = MixingMap(i,j)%rCutWasSet
             MixingMap(j,i)%shiftedPot = MixingMap(i,j)%shiftedPot
             MixingMap(j,i)%isSoftCore = MixingMap(i,j)%isSoftCore
          endif

       enddo
    enddo
    
    haveMixingMap = .true.
    
  end subroutine createMixingMap
  
  subroutine do_lj_pair(atom1, atom2, d, rij, r2, rcut, sw, vpair, fpair, &
       pot, f, do_pot)
    
    integer, intent(in) ::  atom1, atom2
    integer :: atid1, atid2, ljt1, ljt2
    real( kind = dp ), intent(in) :: rij, r2, rcut
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), dimension(3,nLocal) :: f    
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    logical, intent(in) :: do_pot

    ! local Variables
    real( kind = dp ) :: drdx, drdy, drdz
    real( kind = dp ) :: fx, fy, fz
    real( kind = dp ) :: pot_temp, dudr
    real( kind = dp ) :: sigma6
    real( kind = dp ) :: epsilon
    real( kind = dp ) :: r6, rc6
    real( kind = dp ) :: t6, tp6
    real( kind = dp ) :: t12, tp12
    real( kind = dp ) :: delta
    logical :: isSoftCore, shiftedPot
    integer :: id1, id2, localError

    if (.not.haveMixingMap) then
       call createMixingMap()
    endif

    ! Look up the correct parameters in the mixing matrix
#ifdef IS_MPI
    atid1 = atid_Row(atom1)
    atid2 = atid_Col(atom2)
#else
    atid1 = atid(atom1)
    atid2 = atid(atom2)
#endif

    ljt1 = LJMap%atidToLJtype(atid1)
    ljt2 = LJMap%atidToLJtype(atid2)

    sigma6     = MixingMap(ljt1,ljt2)%sigma6
    epsilon    = MixingMap(ljt1,ljt2)%epsilon
    delta      = MixingMap(ljt1,ljt2)%delta
    isSoftCore = MixingMap(ljt1,ljt2)%isSoftCore
    shiftedPot = MixingMap(ljt1,ljt2)%shiftedPot

    r6 = r2 * r2 * r2

    t6  = sigma6/ r6
    t12 = t6 * t6     

    if (isSoftCore) then
       
       pot_temp = 4.0E0_DP * epsilon * t6
       if (shiftedPot) then
          rc6 = rcut**6
          tp6 = sigma6 / rc6
          delta =-4.0_DP*epsilon*(tp6)
          pot_temp = pot_temp + delta
       endif
       
       vpair = vpair + pot_temp
       
       dudr = -sw * 24.0E0_DP * epsilon * t6 / rij
       
    else
       pot_temp = 4.0E0_DP * epsilon * (t12 - t6) 
       if (shiftedPot) then
          rc6 = rcut**6
          tp6 = sigma6 / rc6
          tp12 = tp6*tp6
          delta =-4.0_DP*epsilon*(tp12 - tp6)
          pot_temp = pot_temp + delta
       endif
       
       vpair = vpair + pot_temp
       
       dudr = sw * 24.0E0_DP * epsilon * (t6 - 2.0E0_DP*t12) / rij
    endif

    drdx = d(1) / rij
    drdy = d(2) / rij
    drdz = d(3) / rij

    fx = dudr * drdx
    fy = dudr * drdy
    fz = dudr * drdz

#ifdef IS_MPI
    if (do_pot) then
       pot_Row(VDW_POT,atom1) = pot_Row(VDW_POT,atom1) + sw*pot_temp*0.5
       pot_Col(VDW_POT,atom2) = pot_Col(VDW_POT,atom2) + sw*pot_temp*0.5
    endif

    f_Row(1,atom1) = f_Row(1,atom1) + fx 
    f_Row(2,atom1) = f_Row(2,atom1) + fy
    f_Row(3,atom1) = f_Row(3,atom1) + fz

    f_Col(1,atom2) = f_Col(1,atom2) - fx 
    f_Col(2,atom2) = f_Col(2,atom2) - fy
    f_Col(3,atom2) = f_Col(3,atom2) - fz       

#else
    if (do_pot) pot = pot + sw*pot_temp

    f(1,atom1) = f(1,atom1) + fx
    f(2,atom1) = f(2,atom1) + fy
    f(3,atom1) = f(3,atom1) + fz

    f(1,atom2) = f(1,atom2) - fx
    f(2,atom2) = f(2,atom2) - fy
    f(3,atom2) = f(3,atom2) - fz
#endif

#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif

    if (molMembershipList(id1) .ne. molMembershipList(id2)) then

       fpair(1) = fpair(1) + fx
       fpair(2) = fpair(2) + fy
       fpair(3) = fpair(3) + fz

    endif

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

end module lj
