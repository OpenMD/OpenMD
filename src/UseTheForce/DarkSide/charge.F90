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

module charge_charge
  
  use force_globals
  use definitions
  use atype_module
  use vector_class
  use simulation
  use status
#ifdef IS_MPI
  use mpiSimulation
#endif
  implicit none

  PRIVATE

  real(kind=dp), parameter :: pre = 332.06508_DP  
  logical, save :: haveChargeMap = .false.

  public :: newChargeType
  public :: do_charge_pair
  public :: getCharge
  
  type :: ChargeList
     integer :: c_ident
     real(kind=DP) :: charge = 0.0_DP
  end type ChargeList

  type(ChargeList), dimension(:), allocatable :: ChargeMap

contains

  subroutine newChargeType(c_ident, charge, status)
    integer,intent(in) :: c_ident
    real(kind=dp),intent(in) :: charge
    integer,intent(out) :: status
    integer :: nAtypes, myATID

    status = 0
    
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)
    
    !! Be simple-minded and assume that we need a ChargeMap that
    !! is the same size as the total number of atom types

    if (.not.allocated(ChargeMap)) then
       
       nAtypes = getSize(atypes)
    
       if (nAtypes == 0) then
          status = -1
          return
       end if
       
       if (.not. allocated(ChargeMap)) then
          allocate(ChargeMap(nAtypes))
       endif
       
    end if

    if (myATID .gt. size(ChargeMap)) then
       status = -1
       return
    endif
    
    ! set the values for ChargeMap for this atom type:

    ChargeMap(myATID)%c_ident = c_ident
    ChargeMap(myATID)%charge = charge
    
  end subroutine newChargeType
  
  function getCharge(atid) result (c)
    integer, intent(in) :: atid
    integer :: localError
    real(kind=dp) :: c
    
    if (.not.allocated(ChargeMap)) then
       call handleError("charge_charge", "no ChargeMap was present before first call of getCharge!")
       return
    end if
    
    c = ChargeMap(atid)%charge
  end function getCharge
      
  subroutine do_charge_pair(atom1, atom2, d, rij, r2, sw, vpair, fpair, &
       pot, f, do_pot)
    
    logical :: do_pot

    integer atom1, atom2, me1, me2, id1, id2
    integer :: localError
    real(kind=dp) :: rij, r2, q1, q2, sw, vpair
    real(kind=dp) :: drdx, drdy, drdz, dudr, fx, fy, fz
    real(kind=dp) :: vterm

    real( kind = dp ) :: pot
    real( kind = dp ), dimension(3) :: d, fpair
    real( kind = dp ), dimension(3,nLocal) :: f
    

    if (.not.allocated(ChargeMap)) then
       call handleError("charge_charge", "no ChargeMap was present before first call of do_charge_pair!")
       return
    end if

#ifdef IS_MPI
    me1 = atid_Row(atom1)
    me2 = atid_Col(atom2)
#else
    me1 = atid(atom1)
    me2 = atid(atom2)
#endif

    q1 = ChargeMap(me1)%charge
    q2 = ChargeMap(me2)%charge

    vterm = pre * q1 * q2 / rij 

    dudr =  -sw * vterm  / rij

    drdx = d(1) / rij
    drdy = d(2) / rij
    drdz = d(3) / rij
    
    fx = dudr * drdx
    fy = dudr * drdy
    fz = dudr * drdz          

    vpair = vpair + vterm
    
#ifdef IS_MPI
    if (do_pot) then
       pot_Row(atom1) = pot_Row(atom1) + sw*vterm*0.5
       pot_Col(atom2) = pot_Col(atom2) + sw*vterm*0.5
    endif
    
    f_Row(1,atom1) = f_Row(1,atom1) + fx
    f_Row(2,atom1) = f_Row(2,atom1) + fy
    f_Row(3,atom1) = f_Row(3,atom1) + fz
    
    f_Col(1,atom2) = f_Col(1,atom2) - fx
    f_Col(2,atom2) = f_Col(2,atom2) - fy
    f_Col(3,atom2) = f_Col(3,atom2) - fz
    
#else

    if (do_pot) pot = pot + sw*vterm
    
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
  end subroutine do_charge_pair
  
end module charge_charge

subroutine newChargeType(ident, charge, status)

  use charge_charge, ONLY : module_newChargeType => newChargeType

  integer, parameter :: DP = selected_real_kind(15)
  integer,intent(inout) :: ident
  real(kind=dp),intent(inout) :: charge
  integer,intent(inout) :: status
  
  call module_newChargeType(ident, charge, status)
  
end subroutine newChargeType
