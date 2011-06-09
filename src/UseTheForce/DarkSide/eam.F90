!!
!! Copyright (c) 2005, 2009 The University of Notre Dame. All Rights Reserved.
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
module eam
  use definitions
  use simulation
  use force_globals
  use status
  use atype_module
  use vector_class
  use interpolation
  implicit none
  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"

  logical, save :: EAM_FF_initialized = .false.
  integer, save :: EAM_Mixing_Policy
  real(kind = dp), save :: EAM_rcut
  logical, save :: haveRcut = .false.

  character(len = statusMsgSize) :: errMesg
  integer :: eam_err

  character(len = 200) :: errMsg
  character(len=*), parameter :: RoutineName =  "EAM MODULE"
  !! Logical that determines if eam arrays should be zeroed
  logical :: cleanme = .true.

  type, private :: EAMtype
     integer           :: eam_atype       
     real( kind = DP ) :: eam_lattice        
     real( kind = DP ) :: eam_rcut     
     integer           :: eam_atype_map
     type(cubicSpline) :: rho 
     type(cubicSpline) :: Z
     type(cubicSpline) :: F
     type(cubicSpline) :: phi
  end type EAMtype


  type, private :: EAMTypeList
     integer           :: n_eam_types = 0
     integer           :: currentAddition = 0

     type (EAMtype), pointer  :: EAMParams(:) => null()
     integer, pointer         :: atidToEAMType(:) => null()
  end type EAMTypeList

  type (eamTypeList), save :: EAMList

  !! standard eam stuff  


  public :: setCutoffEAM
  public :: do_eam_pair
  public :: newEAMtype
  public :: calc_eam_prepair_rho
  public :: calc_eam_preforce_Frho
  public :: destroyEAMTypes
  public :: getEAMCut
  public :: lookupEAMSpline
  public :: lookupEAMSpline1d

contains

  subroutine newEAMtype(lattice_constant,eam_nrho,eam_drho,eam_nr,&
       eam_dr,rcut,eam_Z_r,eam_rho_r,eam_F_rho, c_ident, status)
    real (kind = dp )                      :: lattice_constant
    integer                                :: eam_nrho
    real (kind = dp )                      :: eam_drho
    integer                                :: eam_nr
    real (kind = dp )                      :: eam_dr
    real (kind = dp )                      :: rcut
    real (kind = dp ), dimension(eam_nr)   :: eam_Z_r, rvals
    real (kind = dp ), dimension(eam_nr)   :: eam_rho_r, eam_phi_r
    real (kind = dp ), dimension(eam_nrho) :: eam_F_rho, rhovals
    integer                                :: c_ident
    integer                                :: status

    integer                                :: nAtypes,nEAMTypes,myATID
    integer                                :: maxVals
    integer                                :: alloc_stat
    integer                                :: current, j
    integer,pointer                        :: Matchlist(:) => null()

    status = 0

    !! Assume that atypes has already been set and get the total number of types in atypes
    !! Also assume that every member of atypes is a EAM model.

    ! check to see if this is the first time into 
    if (.not.associated(EAMList%EAMParams)) then
       call getMatchingElementList(atypes, "is_EAM", .true., nEAMtypes, MatchList)
       EAMList%n_eam_types = nEAMtypes
       allocate(EAMList%EAMParams(nEAMTypes))
       nAtypes = getSize(atypes)
       allocate(EAMList%atidToEAMType(nAtypes))
       EAMList%atidToEAMType = -1
    end if

    EAMList%currentAddition = EAMList%currentAddition + 1
    current = EAMList%currentAddition

    myATID =  getFirstMatchingElement(atypes, "c_ident", c_ident)
    EAMList%atidToEAMType(myATID) = current

    EAMList%EAMParams(current)%eam_atype    = c_ident
    EAMList%EAMParams(current)%eam_lattice  = lattice_constant
    EAMList%EAMParams(current)%eam_rcut     = rcut

    ! Build array of r values
    do j = 1, eam_nr
       rvals(j) = real(j-1,kind=dp) * eam_dr
    end do
    
    ! Build array of rho values
    do j = 1, eam_nrho
       rhovals(j) = real(j-1,kind=dp) * eam_drho
    end do

    ! convert from eV to kcal / mol:
    do j = 1, eam_nrho
       eam_F_rho(j) = eam_F_rho(j) * 23.06054E0_DP
    end do
    
    ! precompute the pair potential and get it into kcal / mol:
    eam_phi_r(1) = 0.0E0_DP
    do j = 2, eam_nr
       eam_phi_r(j) = 331.999296E0_DP * (eam_Z_r(j)**2) / rvals(j)
    enddo

    call newSpline(EAMList%EAMParams(current)%rho, rvals, eam_rho_r, .true.)
    call newSpline(EAMList%EAMParams(current)%Z,   rvals, eam_Z_r, .true.)
    call newSpline(EAMList%EAMParams(current)%F, rhovals, eam_F_rho, .true.)
    call newSpline(EAMList%EAMParams(current)%phi, rvals, eam_phi_r, .true.)
  end subroutine newEAMtype


  ! kills all eam types entered and sets simulation to uninitalized
  subroutine destroyEAMtypes()
    integer :: i
    type(EAMType), pointer :: tempEAMType=>null()

    do i = 1, EAMList%n_eam_types
       tempEAMType => eamList%EAMParams(i)
       call deallocate_EAMType(tempEAMType)
    end do
    if(associated( eamList%EAMParams)) deallocate( eamList%EAMParams)
    eamList%EAMParams => null()

    eamList%n_eam_types = 0
    eamList%currentAddition = 0
  end subroutine destroyEAMtypes

  function getEAMCut(atomID) result(cutValue)
    integer, intent(in) :: atomID
    integer :: eamID
    real(kind=dp) :: cutValue
    
    eamID = EAMList%atidToEAMType(atomID)
    cutValue = EAMList%EAMParams(eamID)%eam_rcut
  end function getEAMCut


  subroutine setCutoffEAM(rcut)
    real(kind=dp) :: rcut
    EAM_rcut = rcut
  end subroutine setCutoffEAM


  subroutine deallocate_EAMType(thisEAMType)
    type (EAMtype), pointer :: thisEAMType

    call deleteSpline(thisEAMType%F)
    call deleteSpline(thisEAMType%rho)
    call deleteSpline(thisEAMType%phi)
    call deleteSpline(thisEAMType%Z)

  end subroutine deallocate_EAMType

  !! Calculates rho_r
  subroutine calc_eam_prepair_rho(atid1, atid2, d, r, rijsq, rho_i_at_j, rho_j_at_i)
    integer :: Atid1, Atid2
    real(kind = dp), dimension(3) :: d
    real(kind = dp), intent(inout)               :: r
    real(kind = dp), intent(inout)               :: rijsq
    ! value of electron density rho do to atom i at atom j
    real(kind = dp), intent(inout) :: rho_i_at_j
    ! value of electron density rho do to atom j at atom i
    real(kind = dp), intent(inout) :: rho_j_at_i
    real(kind = dp) :: rci, rcj
    integer :: eam_err
   
    integer :: myid_atom1 ! EAM atid
    integer :: myid_atom2 

    ! check to see if we need to be cleaned at the start of a force loop

    Myid_atom1 = Eamlist%atidtoeamtype(Atid1)
    Myid_atom2 = Eamlist%atidtoeamtype(Atid2)
    
    if (r.lt.EAMList%EAMParams(myid_atom1)%eam_rcut) then
       
       ! get cutoff for atom 1
       rci = EAMList%EAMParams(myid_atom1)%eam_rcut
       ! get type specific cutoff for atom 2
       rcj = EAMList%EAMParams(myid_atom2)%eam_rcut
       
       if (r.lt.rci) then
          call lookupEAMSpline(EAMList%EAMParams(myid_atom1)%rho, r, &
               rho_i_at_j)
       endif

       if (r.lt.rcj) then
          call lookupEAMSpline(EAMList%EAMParams(myid_atom2)%rho, r, &
               rho_j_at_i)
       endif
    endif
  end subroutine calc_eam_prepair_rho


  !! Calculate the functional F(rho) for all local atoms
  subroutine calc_eam_preforce_Frho(nlocal, pot, particle_pot)
    integer :: nlocal
    real(kind=dp) :: pot
    integer :: i, j
    integer :: atom
    real(kind=dp) :: U,U1
    real( kind = dp ), dimension(nlocal) :: particle_pot
    integer :: atype1
    integer :: me, atid1

!    cleanme = .true.
    !! Calculate F(rho) and derivative 
    do atom = 1, nlocal
       atid1 = atid(atom)
       me = EAMList%atidtoEAMtype(atid1)
       ! me is set to -1 for non EAM atoms.
       ! Punt if we are a non-EAM atom type.
       if (me == -1) then
          frho(atom) = 0.0_dp
          dfrhodrho(atom) = 0.0_dp
       else          
          
          call lookupEAMSpline1d(EAMList%EAMParams(me)%F, rho(atom), &
               u, u1)
       
          frho(atom) = u
          dfrhodrho(atom) = u1
       endif
       pot = pot + u
       particle_pot(atom) = particle_pot(atom) + u

    enddo

  end subroutine calc_eam_preforce_Frho
  
  !! Does EAM pairwise Force calculation.  
  subroutine do_eam_pair(atid1, atid2, d, rij, r2, sw, vpair, &
       fpair, pot, f1, rho_i, rho_j, dfrhodrho_i, dfrhodrho_j, &
       fshift_i, fshift_j)
    !Arguments    
    integer, intent(in) ::  atid1, atid2
    real( kind = dp ), intent(in) :: rij, r2
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), dimension(3) :: f1
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    real( kind = dp ), intent(inout) :: dfrhodrho_i, dfrhodrho_j
    real( kind = dp ), intent(inout) :: rho_i, rho_j
    real( kind = dp ), intent(inout):: fshift_i, fshift_j

    real( kind = dp ) :: drdx, drdy, drdz
    real( kind = dp ) :: phab, pha, dvpdr
    real( kind = dp ) :: rha, drha, dpha
    real( kind = dp ) :: rhb, drhb, dphb
    real( kind = dp ) :: dudr
    real( kind = dp ) :: rci, rcj
    real( kind = dp ) :: drhoidr, drhojdr
    real( kind = dp ) :: Fx, Fy, Fz
    real( kind = dp ) :: r, phb
    real( kind = dp ) :: u1, u2

    integer :: id1, id2
    integer  :: mytype_atom1
    integer  :: mytype_atom2

    phab = 0.0E0_DP
    dvpdr = 0.0E0_DP
    drha = 0.0
    drhb = 0.0

    if (rij .lt. EAM_rcut) then

       mytype_atom1 = EAMList%atidToEAMType(atid1)
       mytype_atom2 = EAMList%atidTOEAMType(atid2)


       ! get cutoff for atom 1
       rci = EAMList%EAMParams(mytype_atom1)%eam_rcut
       ! get type specific cutoff for atom 2
       rcj = EAMList%EAMParams(mytype_atom2)%eam_rcut

       drdx = d(1)/rij
       drdy = d(2)/rij
       drdz = d(3)/rij

       if (rij.lt.rci) then

          ! Calculate rho and drho for atom1

          call lookupEAMSpline1d(EAMList%EAMParams(mytype_atom1)%rho, &
               rij, rha, drha)
          
          ! Calculate Phi(r) for atom1.
          
          call lookupEAMSpline1d(EAMList%EAMParams(mytype_atom1)%phi, &
               rij, pha, dpha)

       endif

       if (rij.lt.rcj) then      

          ! Calculate rho and drho for atom2

          call lookupEAMSpline1d(EAMList%EAMParams(mytype_atom2)%rho, &
               rij, rhb, drhb)

          ! Calculate Phi(r) for atom2.

          call lookupEAMSpline1d(EAMList%EAMParams(mytype_atom2)%phi, &
               rij, phb, dphb)

       endif

       if (rij.lt.rci) then 
          phab = phab + 0.5E0_DP*(rhb/rha)*pha
          dvpdr = dvpdr + 0.5E0_DP*((rhb/rha)*dpha + &
               pha*((drhb/rha) - (rhb*drha/rha/rha)))
       endif

       if (rij.lt.rcj) then
          phab = phab + 0.5E0_DP*(rha/rhb)*phb
          dvpdr = dvpdr + 0.5E0_DP*((rha/rhb)*dphb + &
               phb*((drha/rhb) - (rha*drhb/rhb/rhb)))
       endif

       drhoidr = drha
       drhojdr = drhb

       dudr = drhojdr*dfrhodrho_i + drhoidr*dfrhodrho_j + dvpdr        

       fx = dudr * drdx
       fy = dudr * drdy
       fz = dudr * drdz

       ! particle_pot is the difference between the full potential 
       ! and the full potential without the presence of a particular
       ! particle (atom1).
       !
       ! This reduces the density at other particle locations, so
       ! we need to recompute the density at atom2 assuming atom1
       ! didn't contribute.  This then requires recomputing the
       ! density functional for atom2 as well.
       !
       ! Most of the particle_pot heavy lifting comes from the
       ! pair interaction, and will be handled by vpair.
       
       call lookupEAMSpline1d(EAMList%EAMParams(mytype_atom1)%F, &
            rho_i-rhb, fshift_i, u1)
       call lookupEAMSpline1d(EAMList%EAMParams(mytype_atom2)%F, &
            rho_j-rha, fshift_j, u2)            

       pot = pot + phab

       f1(1) = f1(1) + fx
       f1(2) = f1(2) + fy
       f1(3) = f1(3) + fz

       vpair = vpair + phab

    endif
  end subroutine do_eam_pair

  subroutine lookupEAMSpline(cs, xval, yval)
    
    implicit none

    type (cubicSpline), intent(in) :: cs
    real( kind = DP ), intent(in)  :: xval
    real( kind = DP ), intent(out) :: yval
    real( kind = DP ) ::  dx
    integer :: i, j
    !
    !  Find the interval J = [ cs%x(J), cs%x(J+1) ] that contains 
    !  or is nearest to xval.
    
    j = MAX(1, MIN(cs%n-1, int((xval-cs%x(1)) * cs%dx_i) + 1))
    
    dx = xval - cs%x(j)
    yval = cs%y(j) + dx*(cs%b(j) + dx*(cs%c(j) + dx*cs%d(j)))
    
    return
  end subroutine lookupEAMSpline

  subroutine lookupEAMSpline1d(cs, xval, yval, dydx)
    
    implicit none

    type (cubicSpline), intent(in) :: cs
    real( kind = DP ), intent(in)  :: xval
    real( kind = DP ), intent(out) :: yval, dydx
    real( kind = DP ) :: dx
    integer :: i, j
    
    !  Find the interval J = [ cs%x(J), cs%x(J+1) ] that contains 
    !  or is nearest to xval.


    j = MAX(1, MIN(cs%n-1, int((xval-cs%x(1)) * cs%dx_i) + 1))
    
    dx = xval - cs%x(j)
    yval = cs%y(j) + dx*(cs%b(j) + dx*(cs%c(j) + dx*cs%d(j)))

    dydx = cs%b(j) + dx*(2.0d0 * cs%c(j) + 3.0d0 * dx * cs%d(j))
       
    return
  end subroutine lookupEAMSpline1d

end module eam
