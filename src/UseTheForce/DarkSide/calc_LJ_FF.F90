!! Calculates Long Range forces Lennard-Jones interactions.
!! Corresponds to the force field defined in lj_FF.cpp 
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @version $Id: calc_LJ_FF.F90,v 1.2 2004-09-24 16:27:57 tim Exp $, $Date: 2004-09-24 16:27:57 $, $Name: not supported by cvs2svn $, $Revision: 1.2 $

module lj
  use definitions
  use atype_module
  use switcheroo
  use vector_class
  use simulation
#ifdef IS_MPI
  use mpiSimulation
#endif
  use force_globals

  implicit none
  PRIVATE

#define __FORTRAN90
#include "UseTheForce/fForceField.h"

  integer, save :: LJ_Mixing_Policy
  real(kind=DP), save :: LJ_rcut
  logical, save :: havePolicy = .false.
  logical, save :: haveCut = .false.
  logical, save :: LJ_do_shift = .false.
  
  !! Logical has lj force field module been initialized?
  
  logical, save :: LJ_FF_initialized = .false.
  
  !! Public methods and data
  public :: init_LJ_FF
  public :: setCutoffLJ
  public :: do_lj_pair
  
  type :: lj_mixed_params
     !! Lennard-Jones epsilon
     real ( kind = dp )  :: epsilon = 0.0_dp
     !! Lennard-Jones Sigma
     real ( kind = dp )  :: sigma = 0.0_dp
     !! Lennard-Jones Sigma to sixth
     real ( kind = dp )  :: sigma6 = 0.0_dp
     !! 
     real ( kind = dp )  :: tp6
     real ( kind = dp )  :: tp12
     real ( kind = dp )  :: delta  = 0.0_dp
  end type lj_mixed_params
  
  type (lj_mixed_params), dimension(:,:), pointer :: ljMixed
  
contains

  subroutine init_LJ_FF(mix_Policy, status)
    integer, intent(in) :: mix_Policy
    integer, intent(out) :: status
    integer :: myStatus
    
    if (mix_Policy == LB_MIXING_RULE) then
       LJ_Mixing_Policy = LB_MIXING_RULE
    else
       if (mix_Policy == EXPLICIT_MIXING_RULE) then
          LJ_Mixing_Policy = EXPLICIT_MIXING_RULE
       else
          write(*,*) 'Unknown Mixing Policy!'
          status = -1
          return
       endif
    endif

    havePolicy = .true.

    if (haveCut) then
       status = 0
       call createMixingList(myStatus)
       if (myStatus /= 0) then
          status = -1
          return
       end if
       
       LJ_FF_initialized = .true.
    end if
  
  end subroutine init_LJ_FF
  
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
    haveCut = .true.

    if (havePolicy) then
       status = 0
       call createMixingList(myStatus)
       if (myStatus /= 0) then
          status = -1
          return
       end if
       
       LJ_FF_initialized = .true.
    end if    
    
    return
  end subroutine setCutoffLJ
  
  subroutine createMixingList(status)
    integer :: nAtypes
    integer :: status
    integer :: i
    integer :: j
    real ( kind = dp ) :: mySigma_i,mySigma_j
    real ( kind = dp ) :: myEpsilon_i,myEpsilon_j
    real ( kind = dp ) :: rcut6
    logical :: I_isLJ, J_isLJ
    status = 0
    
    nAtypes = getSize(atypes)
    if (nAtypes == 0) then
       status = -1
       return
    end if
        
    if (.not. associated(ljMixed)) then
       allocate(ljMixed(nAtypes, nAtypes))
    endif

    rcut6 = LJ_rcut**6

! This loops through all atypes, even those that don't support LJ forces.
    do i = 1, nAtypes

       call getElementProperty(atypes, i, "is_LJ", I_isLJ)

       if (I_isLJ) then

          call getElementProperty(atypes, i, "lj_epsilon", myEpsilon_i)
          call getElementProperty(atypes, i, "lj_sigma",   mySigma_i)
          ! do self mixing rule
          ljMixed(i,i)%sigma   = mySigma_i
          
          ljMixed(i,i)%sigma6  = (ljMixed(i,i)%sigma) ** 6
          
          ljMixed(i,i)%tp6     = (ljMixed(i,i)%sigma6)/rcut6
          
          ljMixed(i,i)%tp12    = (ljMixed(i,i)%tp6) ** 2
          
          
          ljMixed(i,i)%epsilon = myEpsilon_i
          
          ljMixed(i,i)%delta = -4.0_DP * ljMixed(i,i)%epsilon * &
            (ljMixed(i,i)%tp12 - ljMixed(i,i)%tp6)
          
          do j = i + 1, nAtypes

             call getElementProperty(atypes, j, "is_LJ", J_isLJ)
             
             if (J_isLJ) then                
             
                call getElementProperty(atypes,j,"lj_epsilon",myEpsilon_j)
                call getElementProperty(atypes,j,"lj_sigma",  mySigma_j)
                
                ljMixed(i,j)%sigma  =  &
                     calcLJMix("sigma",mySigma_i, &
                     mySigma_j)
                
                ljMixed(i,j)%sigma6 = &
                     (ljMixed(i,j)%sigma)**6
                
                
                ljMixed(i,j)%tp6     = ljMixed(i,j)%sigma6/rcut6
                
                ljMixed(i,j)%tp12    = (ljMixed(i,j)%tp6) ** 2
                
                
                ljMixed(i,j)%epsilon = &
                     calcLJMix("epsilon",myEpsilon_i, &
                     myEpsilon_j)
                
                ljMixed(i,j)%delta = -4.0_DP * ljMixed(i,j)%epsilon * &
                     (ljMixed(i,j)%tp12 - ljMixed(i,j)%tp6)
                
                
                ljMixed(j,i)%sigma   = ljMixed(i,j)%sigma
                ljMixed(j,i)%sigma6  = ljMixed(i,j)%sigma6
                ljMixed(j,i)%tp6     = ljMixed(i,j)%tp6
                ljMixed(j,i)%tp12    = ljMixed(i,j)%tp12
                ljMixed(j,i)%epsilon = ljMixed(i,j)%epsilon
                ljMixed(j,i)%delta   = ljMixed(i,j)%delta
             endif                
          end do
       endif
    end do
    
  end subroutine createMixingList
  
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
    integer :: id1, id2

    ! Look up the correct parameters in the mixing matrix
#ifdef IS_MPI
    sigma6   = ljMixed(atid_Row(atom1),atid_Col(atom2))%sigma6
    epsilon  = ljMixed(atid_Row(atom1),atid_Col(atom2))%epsilon
    delta    = ljMixed(atid_Row(atom1),atid_Col(atom2))%delta
#else
    sigma6   = ljMixed(atid(atom1),atid(atom2))%sigma6
    epsilon  = ljMixed(atid(atom1),atid(atom2))%epsilon
    delta    = ljMixed(atid(atom1),atid(atom2))%delta
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
  
  function calcLJMix(thisParam,param1,param2,status) result(myMixParam)
    character(len=*) :: thisParam
    real(kind = dp)  :: param1
    real(kind = dp)  :: param2
    real(kind = dp ) :: myMixParam

    integer, optional :: status   

    myMixParam = 0.0_dp
    
    if (present(status)) status = 0
    select case (LJ_Mixing_Policy)
    case (1)
       select case (thisParam)
       case ("sigma")
          myMixParam = 0.5_dp * (param1 + param2)
       case ("epsilon")
          myMixParam = sqrt(param1 * param2)
       case default
          status = -1
       end select
    case default
       status = -1
    end select
  end function calcLJMix
  
end module lj
