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
!! @version $Id: LJ.F90,v 1.10 2005-04-13 20:36:45 chuckv Exp $, $Date: 2005-04-13 20:36:45 $, $Name: not supported by cvs2svn $, $Revision: 1.10 $


module lj
  use atype_module
  use switcheroo
  use vector_class
  use simulation
  use status
#ifdef IS_MPI
  use mpiSimulation
#endif
  use force_globals

  implicit none
  PRIVATE
  
  integer, parameter :: DP = selected_real_kind(15)
  
  type, private :: LjType
     integer :: c_ident
     integer :: atid
     real(kind=dp) :: sigma
     real(kind=dp) :: epsilon
     logical :: soft_pot
  end type LjType
  
  type(LjType), dimension(:), allocatable :: ParameterMap
  
  logical, save :: haveMixingMap = .false.
  
  type :: MixParameters
     real(kind=DP) :: sigma
     real(kind=DP) :: epsilon
     real(kind=dp)  :: sigma6
     real(kind=dp)  :: tp6
     real(kind=dp)  :: tp12
     real(kind=dp)  :: delta 
     logical :: soft_pot     
  end type MixParameters
  
  type(MixParameters), dimension(:,:), allocatable :: MixingMap
  
  real(kind=DP), save :: LJ_rcut
  logical, save :: have_rcut = .false.
  logical, save :: LJ_do_shift = .false.
  logical, save :: useGeometricDistanceMixing = .false.
   
  !! Public methods and data
  
  public :: setCutoffLJ
  public :: useGeometricMixing
  public :: do_lj_pair
  public :: newLJtype  
  public :: getSigma
  public :: getEpsilon
  public :: destroyLJTypes
  
contains

  subroutine newLJtype(c_ident, sigma, epsilon, soft_pot, status)
    integer,intent(in) :: c_ident
    real(kind=dp),intent(in) :: sigma
    real(kind=dp),intent(in) :: epsilon
    integer, intent(in) :: soft_pot
    integer,intent(out) :: status
    integer :: nATypes, myATID
    integer, pointer :: MatchList(:) => null()

    status = 0
    
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    if (.not.allocated(ParameterMap)) then
       
       !call getMatchingElementList(atypes, "is_LennardJones", .true., &
       !     nLJTypes, MatchList)
       nAtypes = getSize(atypes)
       if (nAtypes == 0) then
          status = -1
          return
       end if
        
       if (.not. allocated(ParameterMap)) then
          allocate(ParameterMap(nAtypes))
       endif
      
    end if

    if (myATID .gt. size(ParameterMap)) then
       status = -1
       return
    endif
    
    ! set the values for ParameterMap for this atom type:

    ParameterMap(myATID)%c_ident = c_ident
    ParameterMap(myATID)%atid = myATID
    ParameterMap(myATID)%epsilon = epsilon
    ParameterMap(myATID)%sigma = sigma
    if (soft_pot == 1) then
      ParameterMap(myATID)%soft_pot = .true.
    else
      ParameterMap(myATID)%soft_pot = .false.
    endif
  end subroutine newLJtype

  function getSigma(atid) result (s)
    integer, intent(in) :: atid
    integer :: localError
    real(kind=dp) :: s
    
    if (.not.allocated(ParameterMap)) then
       call handleError("LJ", "no ParameterMap was present before first call of getSigma!")
       return
    end if
    
    s = ParameterMap(atid)%sigma
  end function getSigma

  function getEpsilon(atid) result (e)
    integer, intent(in) :: atid
    integer :: localError
    real(kind=dp) :: e
    
    if (.not.allocated(ParameterMap)) then
       call handleError("LJ", "no ParameterMap was present before first call of getEpsilon!")
       return
    end if
    
    e = ParameterMap(atid)%epsilon
  end function getEpsilon


  subroutine setCutoffLJ(rcut, do_shift, status)
    logical, intent(in):: do_shift
    integer :: status, myStatus
    real(kind=dp) :: rcut

#define __FORTRAN90
#include "UseTheForce/fSwitchingFunction.h"

    status = 0

    LJ_rcut = rcut
    LJ_do_shift = do_shift
    call set_switch(LJ_SWITCH, rcut, rcut)
    have_rcut = .true.
    
    return
  end subroutine setCutoffLJ

  subroutine useGeometricMixing() 
    useGeometricDistanceMixing = .true.
    haveMixingMap = .false.
    return
  end subroutine useGeometricMixing
  
  subroutine createMixingMap(status)
    integer :: nATIDs
    integer :: status
    integer :: i
    integer :: j
    real ( kind = dp ) :: Sigma_i, Sigma_j
    real ( kind = dp ) :: Epsilon_i, Epsilon_j
    real ( kind = dp ) :: rcut6
    logical :: i_is_LJ, j_is_LJ

    status = 0

    if (.not. allocated(ParameterMap)) then
       call handleError("LJ", "no ParameterMap was present before call of createMixingMap!")
       status = -1
       return
    endif
    
    nATIDs = size(ParameterMap)
    
    if (nATIDs == 0) then
       status = -1
       return
    end if

    if (.not. allocated(MixingMap)) then
       allocate(MixingMap(nATIDs, nATIDs))
    endif

    if (.not.have_rcut) then
       status = -1
       return
    endif
        
    rcut6 = LJ_rcut**6
    
    ! This loops through all atypes, even those that don't support LJ forces.
    do i = 1, nATIDs
       call getElementProperty(atypes, i, "is_LennardJones", i_is_LJ)
       if (i_is_LJ) then
          Epsilon_i = ParameterMap(i)%epsilon
          Sigma_i = ParameterMap(i)%sigma
          
          ! do self mixing rule
          MixingMap(i,i)%sigma   = Sigma_i          
          MixingMap(i,i)%sigma6  = Sigma_i ** 6          
          MixingMap(i,i)%tp6     = (MixingMap(i,i)%sigma6)/rcut6          
          MixingMap(i,i)%tp12    = (MixingMap(i,i)%tp6) ** 2
          MixingMap(i,i)%epsilon = Epsilon_i          
          MixingMap(i,i)%delta   = -4.0_DP * MixingMap(i,i)%epsilon * &
               (MixingMap(i,i)%tp12 - MixingMap(i,i)%tp6)
          MixingMap(i,i)%soft_pot = ParameterMap(i)%soft_pot
          
          do j = i + 1, nATIDs
             call getElementProperty(atypes, j, "is_LennardJones", j_is_LJ)
          
             if (j_is_LJ) then
                Epsilon_j = ParameterMap(j)%epsilon
                Sigma_j = ParameterMap(j)%sigma
                
                ! only the distance parameter uses different mixing policies
                if (useGeometricDistanceMixing) then
                   ! only for OPLS as far as we can tell
                   MixingMap(i,j)%sigma = dsqrt(Sigma_i * Sigma_j)
                else
                   ! everyone else
                   MixingMap(i,j)%sigma = 0.5_dp * (Sigma_i + Sigma_j)
                endif
                
                ! energy parameter is always geometric mean:
                MixingMap(i,j)%epsilon = dsqrt(Epsilon_i * Epsilon_j)
                
                MixingMap(i,j)%sigma6 = (MixingMap(i,j)%sigma)**6
                MixingMap(i,j)%tp6    = MixingMap(i,j)%sigma6/rcut6
                MixingMap(i,j)%tp12    = (MixingMap(i,j)%tp6) ** 2
                
                MixingMap(i,j)%delta = -4.0_DP * MixingMap(i,j)%epsilon * &
                     (MixingMap(i,j)%tp12 - MixingMap(i,j)%tp6)
                
                MixingMap(i,j)%soft_pot = ParameterMap(i)%soft_pot .or. ParameterMap(j)%soft_pot
                
                
                MixingMap(j,i)%sigma   = MixingMap(i,j)%sigma
                MixingMap(j,i)%sigma6  = MixingMap(i,j)%sigma6
                MixingMap(j,i)%tp6     = MixingMap(i,j)%tp6
                MixingMap(j,i)%tp12    = MixingMap(i,j)%tp12
                MixingMap(j,i)%epsilon = MixingMap(i,j)%epsilon
                MixingMap(j,i)%delta   = MixingMap(i,j)%delta
                MixingMap(j,i)%soft_pot   = MixingMap(i,j)%soft_pot
             endif
          end do
       endif
    end do
    
    haveMixingMap = .true.

  end subroutine createMixingMap
        
  subroutine do_lj_pair(atom1, atom2, d, rij, r2, sw, vpair, fpair, &
       pot, f, do_pot)

    integer, intent(in) ::  atom1, atom2
    real( kind = dp ), intent(in) :: rij, r2
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
    real( kind = dp ) :: r6
    real( kind = dp ) :: t6
    real( kind = dp ) :: t12
    real( kind = dp ) :: delta
    logical :: soft_pot
    integer :: id1, id2, localError

    if (.not.haveMixingMap) then
       localError = 0
       call createMixingMap(localError)
       if ( localError .ne. 0 ) then
          call handleError("LJ", "MixingMap creation failed!")
          return
       end if
    endif

    ! Look up the correct parameters in the mixing matrix
#ifdef IS_MPI
    sigma6   = MixingMap(atid_Row(atom1),atid_Col(atom2))%sigma6
    epsilon  = MixingMap(atid_Row(atom1),atid_Col(atom2))%epsilon
    delta    = MixingMap(atid_Row(atom1),atid_Col(atom2))%delta
    soft_pot =  MixingMap(atid_Row(atom1),atid_Col(atom2))%soft_pot
#else
    sigma6   = MixingMap(atid(atom1),atid(atom2))%sigma6
    epsilon  = MixingMap(atid(atom1),atid(atom2))%epsilon
    delta    = MixingMap(atid(atom1),atid(atom2))%delta
    soft_pot    = MixingMap(atid(atom1),atid(atom2))%soft_pot
#endif

    r6 = r2 * r2 * r2
    
    t6  = sigma6/ r6
    t12 = t6 * t6     

    if (soft_pot) then
      pot_temp = 4.0E0_DP * epsilon * t12 
      if (LJ_do_shift) then
         pot_temp = pot_temp + delta
      endif
  
      vpair = vpair + pot_temp
         
      dudr = sw * 24.0E0_DP * epsilon * (-2.0E0_DP)*t12 / rij

    else
      pot_temp = 4.0E0_DP * epsilon * (t12 - t6) 
      if (LJ_do_shift) then
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
       pot_Row(atom1) = pot_Row(atom1) + sw*pot_temp*0.5
       pot_Col(atom2) = pot_Col(atom2) + sw*pot_temp*0.5
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
    if(allocated(ParameterMap)) deallocate(ParameterMap)
    if(allocated(MixingMap)) deallocate(MixingMap)
    haveMixingMap = .false.
  end subroutine destroyLJTypes

    
end module lj
