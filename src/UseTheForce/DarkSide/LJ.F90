!! Calculates Long Range forces Lennard-Jones interactions.
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @version $Id: LJ.F90,v 1.4 2004-10-22 20:22:47 gezelter Exp $, $Date: 2004-10-22 20:22:47 $, $Name: not supported by cvs2svn $, $Revision: 1.4 $

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
     integer :: ident
     real(kind=dp) :: sigma
     real(kind=dp) :: epsilon
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
  
contains

  subroutine newLJtype(ident, sigma, epsilon, status)
    integer,intent(in) :: ident
    real(kind=dp),intent(in) :: sigma
    real(kind=dp),intent(in) :: epsilon
    integer,intent(out) :: status
    integer :: nAtypes

    status = 0
    
    !! Be simple-minded and assume that we need a ParameterMap that
    !! is the same size as the total number of atom types

    if (.not.allocated(ParameterMap)) then
       
       nAtypes = getSize(atypes)
    
       if (nAtypes == 0) then
          status = -1
          return
       end if
       
       if (.not. allocated(ParameterMap)) then
          allocate(ParameterMap(nAtypes))
       endif
       
    end if

    if (ident .gt. size(ParameterMap)) then
       status = -1
       return
    endif
    
    ! set the values for ParameterMap for this atom type:

    ParameterMap(ident)%ident = ident
    ParameterMap(ident)%epsilon = epsilon
    ParameterMap(ident)%sigma = sigma
    
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
       call handleError("dipole-dipole", "no ParameterMap was present before first call of getEpsilon!")
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
    integer :: nAtypes
    integer :: status
    integer :: i
    integer :: j
    real ( kind = dp ) :: Sigma_i, Sigma_j
    real ( kind = dp ) :: Epsilon_i, Epsilon_j
    real ( kind = dp ) :: rcut6

    status = 0
    
    nAtypes = size(ParameterMap)
    
    if (nAtypes == 0) then
       status = -1
       return
    end if

    if (.not.have_rcut) then
       status = -1
       return
    endif
    
    if (.not. allocated(MixingMap)) then
       allocate(MixingMap(nAtypes, nAtypes))
    endif
    
    rcut6 = LJ_rcut**6
    
    ! This loops through all atypes, even those that don't support LJ forces.
    do i = 1, nAtypes
       
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
       
       do j = i + 1, nAtypes
          
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
          
          MixingMap(j,i)%sigma   = MixingMap(i,j)%sigma
          MixingMap(j,i)%sigma6  = MixingMap(i,j)%sigma6
          MixingMap(j,i)%tp6     = MixingMap(i,j)%tp6
          MixingMap(j,i)%tp12    = MixingMap(i,j)%tp12
          MixingMap(j,i)%epsilon = MixingMap(i,j)%epsilon
          MixingMap(j,i)%delta   = MixingMap(i,j)%delta
          
       end do
    end do
    
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
#else
    sigma6   = MixingMap(atid(atom1),atid(atom2))%sigma6
    epsilon  = MixingMap(atid(atom1),atid(atom2))%epsilon
    delta    = MixingMap(atid(atom1),atid(atom2))%delta
#endif

    r6 = r2 * r2 * r2
    
    t6  = sigma6/ r6
    t12 = t6 * t6     
   
    pot_temp = 4.0E0_DP * epsilon * (t12 - t6) 
    if (LJ_do_shift) then
       pot_temp = pot_temp + delta
    endif

    vpair = vpair + pot_temp
       
    dudr = sw * 24.0E0_DP * epsilon * (t6 - 2.0E0_DP*t12) / rij
       
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
  
  
  !! Calculates the mixing for sigma or epslon
    
end module lj

subroutine newLJtype(ident, sigma, epsilon, status)
  use lj, ONLY : module_newLJtype => newLJtype
  integer, parameter :: DP = selected_real_kind(15)
  integer,intent(inout) :: ident
  real(kind=dp),intent(inout) :: sigma
  real(kind=dp),intent(inout) :: epsilon
  integer,intent(inout) :: status
  
  call module_newLJtype(ident, sigma, epsilon, status)
  
end subroutine newLJtype

subroutine useGeometricMixing()
  use lj, ONLY: module_useGeometricMixing => useGeometricMixing
  
  call module_useGeometricMixing()
  return
end subroutine useGeometricMixing
