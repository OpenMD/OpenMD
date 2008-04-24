!!
!! Copyright (c) 2007 The University of Notre Dame. All Rights Reserved.
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


!! Calculates Metal-Non Metal interactions.
!! @author Charles F. Vardeman II 
!! @version $Id: MetalNonMetal.F90,v 1.11 2008-04-24 21:04:27 gezelter Exp $, $Date: 2008-04-24 21:04:27 $, $Name: not supported by cvs2svn $, $Revision: 1.11 $


module MetalNonMetal
  use definitions
  use atype_module
  use vector_class
  use simulation
  use status
  use fForceOptions
#ifdef IS_MPI
  use mpiSimulation
#endif
  use force_globals

  implicit none
  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"
#include "UseTheForce/DarkSide/fMnMInteractions.h"

  logical, save :: useGeometricDistanceMixing = .false.
  logical, save :: haveInteractionLookup = .false.

  real(kind=DP), save :: defaultCutoff = 0.0_DP
  
  logical, save :: defaultShiftPot = .false.
  logical, save :: defaultShiftFrc = .false.
  logical, save :: haveDefaultCutoff = .false.

  type :: MnMinteraction
     integer :: metal_atid
     integer :: nonmetal_atid
     integer :: interaction_type
     real(kind=dp) :: R0
     real(kind=dp) :: D0
     real(kind=dp) :: beta0
     real(kind=dp) :: betaH
     real(kind=dp) :: ca1
     real(kind=dp) :: cb1
     real(kind=dp) :: sigma
     real(kind=dp) :: epsilon
     real(kind=dp) :: rCut = 0.0_dp
     logical       :: rCutWasSet = .false.
     logical       :: shiftedPot = .false.
     logical       :: shiftedFrc = .false.
  end type MnMinteraction

  type :: MnMinteractionMap
     PRIVATE
     integer :: initialCapacity = 10
     integer :: capacityIncrement = 0
     integer :: interactionCount = 0
     type(MnMinteraction), pointer :: interactions(:) => null()
  end type MnMinteractionMap

  type (MnMInteractionMap), pointer :: MnM_Map

  integer,  allocatable, dimension(:,:) :: MnMinteractionLookup

  public :: setMnMDefaultCutoff
  public :: addInteraction
  public :: deleteInteractions
  public :: MNMtype
  public :: do_mnm_pair

contains


  subroutine do_mnm_pair(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot, A, F,t, Do_pot)
    integer, intent(in) ::  atom1, atom2
    integer :: atid1, atid2, ljt1, ljt2
    real( kind = dp ), intent(in) :: rij, r2, rcut
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3,nLocal) :: f 
    real (kind=dp), intent(in), dimension(9,nLocal) :: A
    real (kind=dp), intent(inout), dimension(3,nLocal) :: t   
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    logical, intent(in) :: do_pot

    integer :: interaction_id
    integer :: interaction_type

#ifdef IS_MPI
    atid1 = atid_Row(atom1)
    atid2 = atid_Col(atom2)
#else
    atid1 = atid(atom1)
    atid2 = atid(atom2)
#endif

    if(.not.haveInteractionLookup)  then
      call createInteractionLookup(MnM_MAP)
      haveInteractionLookup =.true.
    end if
    
    interaction_id = MnMinteractionLookup(atid1, atid2)
    interaction_type = MnM_Map%interactions(interaction_id)%interaction_type
    

    
    select case (interaction_type)
    case (MNM_LENNARDJONES)
       call calc_mnm_lennardjones(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, &
            Fpair, Pot, F, Do_pot, interaction_id)
    case(MNM_REPULSIVEMORSE, MNM_SHIFTEDMORSE)
       call calc_mnm_morse(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
            Pot, F, Do_pot, interaction_id, interaction_type)
    case(MNM_MAW)
       call calc_mnm_maw(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
            Pot,A, F,t, Do_pot, interaction_id)
    case default
    call handleError("MetalNonMetal","Unknown interaction type")      
    end select

  end subroutine do_mnm_pair

  subroutine calc_mnm_lennardjones(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, &
    Fpair, Pot, F, Do_pot, interaction_id)
    
    integer, intent(in) ::  atom1, atom2
    real( kind = dp ), intent(in) :: rij, r2, rcut
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3,nLocal) :: f    
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    logical, intent(in) :: do_pot
    integer, intent(in) :: interaction_id

    ! local Variables
    real( kind = dp ) :: drdx, drdy, drdz
    real( kind = dp ) :: fx, fy, fz
    real( kind = dp ) :: myPot, myPotC, myDeriv, myDerivC, ros, rcos
    real( kind = dp ) :: pot_temp, dudr
    real( kind = dp ) :: sigmai
    real( kind = dp ) :: epsilon
    logical :: isSoftCore, shiftedPot, shiftedFrc
    integer :: id1, id2, localError

    sigmai     = 1.0_dp / MnM_Map%interactions(interaction_id)%sigma
    epsilon    = MnM_Map%interactions(interaction_id)%epsilon
    shiftedPot = MnM_Map%interactions(interaction_id)%shiftedPot
    shiftedFrc = MnM_Map%interactions(interaction_id)%shiftedFrc

    ros = rij * sigmai

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

    pot_temp = epsilon * (myPot - myPotC)
    vpair = vpair + pot_temp
    dudr = sw * epsilon * (myDeriv - myDerivC) * sigmai

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
  end subroutine calc_mnm_lennardjones

  subroutine calc_mnm_morse(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot, f, Do_pot, interaction_id, interaction_type)
    integer, intent(in) ::  atom1, atom2
    real( kind = dp ), intent(in) :: rij, r2, rcut
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3,nLocal) :: f    
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    logical, intent(in) :: do_pot
    integer, intent(in) :: interaction_id, interaction_type
    logical :: shiftedPot, shiftedFrc

    ! Local Variables
    real(kind=dp) :: Beta0
    real(kind=dp) :: R0
    real(kind=dp) :: D0
    real(kind=dp) :: expt
    real(kind=dp) :: expt2
    real(kind=dp) :: expfnc
    real(kind=dp) :: expfnc2
    real(kind=dp) :: D_expt
    real(kind=dp) :: D_expt2
    real(kind=dp) :: rcos
    real(kind=dp) :: exptC
    real(kind=dp) :: expt2C
    real(kind=dp) :: expfncC
    real(kind=dp) :: expfnc2C
    real(kind=dp) :: D_expt2C
    real(kind=dp) :: D_exptC
    
    real(kind=dp) :: myPot
    real(kind=dp) :: myPotC
    real(kind=dp) :: myDeriv
    real(kind=dp) :: myDerivC
    real(kind=dp) :: pot_temp
    real(kind=dp) :: fx,fy,fz
    real(kind=dp) :: drdx,drdy,drdz
    real(kind=dp) :: dudr
    integer :: id1,id2
    

    D0 = MnM_Map%interactions(interaction_id)%D0
    R0 = MnM_Map%interactions(interaction_id)%r0
    Beta0 = MnM_Map%interactions(interaction_id)%Beta0
    shiftedPot = MnM_Map%interactions(interaction_id)%shiftedPot
    shiftedFrc = MnM_Map%interactions(interaction_id)%shiftedFrc
        
    ! V(r) = D_e exp(-a(r-re)(exp(-a(r-re))-2)
    
    expt     = -beta0*(rij - R0) 
    expfnc   = exp(expt)
    expfnc2  = expfnc*expfnc

    if (shiftedPot .or. shiftedFrc) then
       exptC     = -beta0*(rcut - R0) 
       expfncC   = exp(exptC)
       expfnc2C  = expfncC*expfncC
    endif
       
    select case (interaction_type)
    case (MNM_SHIFTEDMORSE)

       myPot  = D0 * (expfnc2  - 2.0_dp * expfnc)
       myDeriv   = 2.0*D0*beta0*(expfnc - expfnc2)

       if (shiftedPot) then
          myPotC = D0 * (expfnc2C - 2.0_dp * expfncC)
          myDerivC = 0.0_dp
       elseif (shiftedFrc) then
          myPotC = D0 * (expfnc2C - 2.0_dp * expfncC)
          myDerivC  = 2.0*D0*beta0*(expfnc2C - expfnc2C)
          myPotC = myPotC + myDerivC * (rij - rcut) 
       else
          myPotC = 0.0_dp
          myDerivC = 0.0_dp
       endif

    case (MNM_REPULSIVEMORSE)

       myPot  = D0 * expfnc2
       myDeriv  = -2.0_dp * D0 * beta0 * expfnc2

       if (shiftedPot) then
          myPotC = D0 * expfnc2C
          myDerivC = 0.0_dp
       elseif (shiftedFrc) then
          myPotC = D0 * expfnc2C
          myDerivC = -2.0_dp * D0* beta0 * expfnc2C
          myPotC = myPotC + myDerivC * (rij - rcut) 
       else
          myPotC = 0.0_dp
          myDerivC = 0.0_dp
       endif
    end select

    pot_temp = (myPot - myPotC)
    vpair = vpair + pot_temp
    dudr = sw * (myDeriv - myDerivC)

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
  end subroutine calc_mnm_morse
  
  subroutine calc_mnm_maw(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot, A, F,t, Do_pot, interaction_id)
    integer, intent(in) ::  atom1, atom2
    real( kind = dp ), intent(in) :: rij, r2, rcut
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3,nLocal) :: f    
    real (kind=dp),intent(in), dimension(9,nLocal) :: A
    real (kind=dp),intent(inout), dimension(3,nLocal) :: t   

    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    logical, intent(in) :: do_pot

    integer, intent(in) :: interaction_id

    real(kind = dp) :: D0, R0, beta0, ca1, cb1, pot_temp
    real(kind = dp) :: expt0, expfnc0, expfnc02
    real(kind = dp) :: exptH, expfncH, expfncH2
    real(kind = dp) :: x, y, z, x2, y2, z2, r3, r4
    real(kind = dp) :: drdx, drdy, drdz
    real(kind = dp) :: dvdx, dvdy, dvdz
    real(kind = dp) :: Vang, dVangdx, dVangdy, dVangdz
    real(kind = dp) :: dVangdux, dVangduy, dVangduz
    real(kind = dp) :: dVmorsedr
    real(kind = dp) :: Vmorse, dVmorsedx, dVmorsedy, dVmorsedz
    real(kind = dp) :: a1, b1, s
    real(kind = dp) :: da1dx, da1dy, da1dz, da1dux, da1duy, da1duz
    real(kind = dp) :: db1dx, db1dy, db1dz, db1dux, db1duy, db1duz
    real(kind = dp) :: fx, fy, fz, tx, ty, tz, fxl, fyl, fzl
    real(kind = dp), parameter :: st = sqrt(3.0_dp)
    integer :: atid1, atid2, id1, id2
    logical :: shiftedPot, shiftedFrc
    
#ifdef IS_MPI
    atid1 = atid_Row(atom1)
    atid2 = atid_Col(atom2)

    if (atid2.eq.MnM_Map%interactions(interaction_id)%metal_atid) then
       ! rotate the inter-particle separation into the two different 
       ! body-fixed coordinate systems:

       x = A_row(1,atom1)*d(1) + A_row(2,atom1)*d(2) + A_row(3,atom1)*d(3)
       y = A_row(4,atom1)*d(1) + A_row(5,atom1)*d(2) + A_row(6,atom1)*d(3)
       z = A_row(7,atom1)*d(1) + A_row(8,atom1)*d(2) + A_row(9,atom1)*d(3)
    else
       ! negative sign because this is the vector from j to i:
       
       x = -(A_Col(1,atom2)*d(1) + A_Col(2,atom2)*d(2) + A_Col(3,atom2)*d(3))
       y = -(A_Col(4,atom2)*d(1) + A_Col(5,atom2)*d(2) + A_Col(6,atom2)*d(3))
       z = -(A_Col(7,atom2)*d(1) + A_Col(8,atom2)*d(2) + A_Col(9,atom2)*d(3))
    endif
       
#else
    atid1 = atid(atom1)
    atid2 = atid(atom2)
    
    if (atid2.eq.MnM_Map%interactions(interaction_id)%metal_atid) then
       ! rotate the inter-particle separation into the two different 
       ! body-fixed coordinate systems:
       
       x = a(1,atom1)*d(1) + a(2,atom1)*d(2) + a(3,atom1)*d(3)
       y = a(4,atom1)*d(1) + a(5,atom1)*d(2) + a(6,atom1)*d(3)
       z = a(7,atom1)*d(1) + a(8,atom1)*d(2) + a(9,atom1)*d(3)
    else
       ! negative sign because this is the vector from j to i:

       x = -(a(1,atom2)*d(1) + a(2,atom2)*d(2) + a(3,atom2)*d(3))
       y = -(a(4,atom2)*d(1) + a(5,atom2)*d(2) + a(6,atom2)*d(3))
       z = -(a(7,atom2)*d(1) + a(8,atom2)*d(2) + a(9,atom2)*d(3))
    endif

#endif
  
    D0 = MnM_Map%interactions(interaction_id)%D0
    R0 = MnM_Map%interactions(interaction_id)%r0
    beta0 = MnM_Map%interactions(interaction_id)%beta0    
    ca1 = MnM_Map%interactions(interaction_id)%ca1
    cb1 = MnM_Map%interactions(interaction_id)%cb1
    
    shiftedPot = MnM_Map%interactions(interaction_id)%shiftedPot
    shiftedFrc = MnM_Map%interactions(interaction_id)%shiftedFrc   
    
    expt0     = -beta0*(rij - R0) 
    expfnc0   = exp(expt0)
    expfnc02  = expfnc0*expfnc0
  
!!$    if (shiftedPot .or. shiftedFrc) then
!!$       exptC0     = -beta0*(rcut - R0) 
!!$       expfncC0   = exp(exptC0)
!!$       expfncC02  = expfncC0*expfncC0
!!$       exptCH     = -betaH*(rcut - R0) 
!!$       expfncCH   = exp(exptCH)
!!$       expfncCH2  = expfncCH*expfncCH
!!$    endif

    drdx = x / rij
    drdy = y / rij
    drdz = z / rij
    
    x2 = x*x
    y2 = y*y
    z2 = z*z
    r3 = r2*rij
    r4 = r2*r2

    Vmorse = D0 * (expfnc02  - 2.0_dp * expfnc0)

    ! angular modulation of morse part of potential to approximate 
    ! the squares of the two HOMO lone pair orbitals in water:
    !
    ! s = 1 / (4 pi)
    ! a1 = (s + pz)^2 = (1 - sqrt(3)*cos(theta))^2 / (4 pi)
    ! b1 = px^2 = 3 * (sin(theta)*cos(phi))^2 / (4 pi)   

    ! we'll leave out the 4 pi for now

    s = 1.0_dp
    a1 = (1.0_dp + st * z / rij )**2
    b1 = 3.0_dp * (x2 - y2)/ r2

    Vang = s + ca1 * a1 + cb1 * b1
    
    pot_temp = Vmorse*Vang 
         
    vpair = vpair + pot_temp
    
    if (do_pot) then
#ifdef IS_MPI
       pot_row(VDW_POT,atom1) = pot_row(VDW_POT,atom1) + 0.5_dp*pot_temp*sw
       pot_col(VDW_POT,atom2) = pot_col(VDW_POT,atom2) + 0.5_dp*pot_temp*sw
#else
       pot = pot + pot_temp*sw
#endif
    endif

    dVmorsedr = 2.0_dp*D0*beta0*(expfnc0 - expfnc02)

    dVmorsedx = dVmorsedr * drdx
    dVmorsedy = dVmorsedr * drdy
    dVmorsedz = dVmorsedr * drdz
    
    da1dx = 2.0_dp * st * x * z / r3 - 6.0_dp * x * z2 / r4
    da1dy = 2.0_dp * st * y * z / r3 - 6.0_dp * y * z2 / r4
    da1dz = 2.0_dp * st * (x2 + y2) * (st * z - rij ) / r4

    db1dx = 6.0_dp * x * z2 / r4
    db1dy = 6.0_dp * y * z2 / r4
    db1dz = -6.0_dp * (x2 + y2) * z / r4

    dVangdx = ca1 * da1dx + cb1 * db1dx
    dVangdy = ca1 * da1dy + cb1 * db1dy
    dVangdz = ca1 * da1dz + cb1 * db1dz
    
    ! chain rule to put these back on x, y, z
    dvdx = Vang * dVmorsedx + Vmorse * dVangdx
    dvdy = Vang * dVmorsedy + Vmorse * dVangdy
    dvdz = Vang * dVmorsedz + Vmorse * dVangdz
    
    ! Torques for Vang using method of Price:
    ! S. L. Price, A. J. Stone, and M. Alderton, Mol. Phys. 52, 987 (1984).

    da1dux =   6.0_dp * y * z / r2 - 2.0_dp * st * y / rij
    da1duy =  -6.0_dp * x * z / r2 + 2.0_dp * st * x / rij
    da1duz =   0.0_dp

    db1dux =   6.0_dp * y * z / r2
    db1duy =   6.0_dp * x * z / r2
    db1duz = -12.0_dp * y * x / r2

    dVangdux = ca1 * da1dux + cb1 * db1dux
    dVangduy = ca1 * da1duy + cb1 * db1duy
    dVangduz = ca1 * da1duz + cb1 * db1duz
    
    ! do the torques first since they are easy:
    ! remember that these are still in the body fixed axes    
    
    tx = Vmorse * dVangdux * sw
    ty = Vmorse * dVangduy * sw
    tz = Vmorse * dVangduz * sw

    ! go back to lab frame using transpose of rotation matrix:

#ifdef IS_MPI
    if (atid2.eq.MnM_Map%interactions(interaction_id)%metal_atid) then
       t_Row(1,atom1) = t_Row(1,atom1) + a_Row(1,atom1)*tx + &
            a_Row(4,atom1)*ty + a_Row(7,atom1)*tz
       t_Row(2,atom1) = t_Row(2,atom1) + a_Row(2,atom1)*tx + &
            a_Row(5,atom1)*ty + a_Row(8,atom1)*tz
       t_Row(3,atom1) = t_Row(3,atom1) + a_Row(3,atom1)*tx + &
            a_Row(6,atom1)*ty + a_Row(9,atom1)*tz
    else
       t_Col(1,atom2) = t_Col(1,atom2) + a_Col(1,atom2)*tx + &
            a_Col(4,atom2)*ty + a_Col(7,atom2)*tz
       t_Col(2,atom2) = t_Col(2,atom2) + a_Col(2,atom2)*tx + &
            a_Col(5,atom2)*ty + a_Col(8,atom2)*tz
       t_Col(3,atom2) = t_Col(3,atom2) + a_Col(3,atom2)*tx + &
            a_Col(6,atom2)*ty + a_Col(9,atom2)*tz
    endif
#else
    if (atid2.eq.MnM_Map%interactions(interaction_id)%metal_atid) then
       t(1,atom1) = t(1,atom1) + a(1,atom1)*tx + a(4,atom1)*ty + &
            a(7,atom1)*tz
       t(2,atom1) = t(2,atom1) + a(2,atom1)*tx + a(5,atom1)*ty + &
            a(8,atom1)*tz
       t(3,atom1) = t(3,atom1) + a(3,atom1)*tx + a(6,atom1)*ty + &
            a(9,atom1)*tz       
    else
       t(1,atom2) = t(1,atom2) + a(1,atom2)*tx + a(4,atom2)*ty + &
            a(7,atom2)*tz
       t(2,atom2) = t(2,atom2) + a(2,atom2)*tx + a(5,atom2)*ty + &
            a(8,atom2)*tz
       t(3,atom2) = t(3,atom2) + a(3,atom2)*tx + a(6,atom2)*ty + &
            a(9,atom2)*tz
    endif
#endif
    ! Now, on to the forces (still in body frame of water)

    fx = dvdx * sw
    fy = dvdy * sw
    fz = dvdz * sw

    ! rotate the terms back into the lab frame:

#ifdef IS_MPI
    if (atid2.eq.MnM_Map%interactions(interaction_id)%metal_atid) then
       fxl = a_Row(1,atom1)*fx + a_Row(4,atom1)*fy + a_Row(7,atom1)*fz
       fyl = a_Row(2,atom1)*fx + a_Row(5,atom1)*fy + a_Row(8,atom1)*fz
       fzl = a_Row(3,atom1)*fx + a_Row(6,atom1)*fy + a_Row(9,atom1)*fz
    else
	   ! negative sign because this is the vector from j to i:
       fxl = -(a_Col(1,atom2)*fx + a_Col(4,atom2)*fy + a_Col(7,atom2)*fz)
       fyl = -(a_Col(2,atom2)*fx + a_Col(5,atom2)*fy + a_Col(8,atom2)*fz)
       fzl = -(a_Col(3,atom2)*fx + a_Col(6,atom2)*fy + a_Col(9,atom2)*fz)
    endif
    f_Row(1,atom1) = f_Row(1,atom1) + fxl
    f_Row(2,atom1) = f_Row(2,atom1) + fyl
    f_Row(3,atom1) = f_Row(3,atom1) + fzl
    
    f_Col(1,atom2) = f_Col(1,atom2) - fxl
    f_Col(2,atom2) = f_Col(2,atom2) - fyl
    f_Col(3,atom2) = f_Col(3,atom2) - fzl
#else
    if (atid2.eq.MnM_Map%interactions(interaction_id)%metal_atid) then
       fxl = a(1,atom1)*fx + a(4,atom1)*fy + a(7,atom1)*fz
       fyl = a(2,atom1)*fx + a(5,atom1)*fy + a(8,atom1)*fz
       fzl = a(3,atom1)*fx + a(6,atom1)*fy + a(9,atom1)*fz
    else
	   ! negative sign because this is the vector from j to i:
       fxl = -(a(1,atom2)*fx + a(4,atom2)*fy + a(7,atom2)*fz)
       fyl = -(a(2,atom2)*fx + a(5,atom2)*fy + a(8,atom2)*fz)
       fzl = -(a(3,atom2)*fx + a(6,atom2)*fy + a(9,atom2)*fz)
    endif
    f(1,atom1) = f(1,atom1) + fxl
    f(2,atom1) = f(2,atom1) + fyl
    f(3,atom1) = f(3,atom1) + fzl

    f(1,atom2) = f(1,atom2) - fxl
    f(2,atom2) = f(2,atom2) - fyl
    f(3,atom2) = f(3,atom2) - fzl       
#endif

#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif
    
    if (molMembershipList(id1) .ne. molMembershipList(id2)) then
       
       fpair(1) = fpair(1) + fxl
       fpair(2) = fpair(2) + fyl
       fpair(3) = fpair(3) + fzl
       
    endif
    
    return
  end subroutine calc_mnm_maw


  subroutine  setMnMDefaultCutoff(thisRcut, shiftedPot, shiftedFrc)
    real(kind=dp), intent(in) :: thisRcut
    logical, intent(in) :: shiftedPot
    logical, intent(in) :: shiftedFrc
    integer i, nInteractions
    defaultCutoff = thisRcut
    defaultShiftPot = shiftedPot
    defaultShiftFrc = shiftedFrc

    if (associated(MnM_Map)) then
       if(MnM_Map%interactionCount /= 0) then
          nInteractions = MnM_Map%interactionCount

          do i = 1, nInteractions
             MnM_Map%interactions(i)%shiftedPot = shiftedPot
             MnM_Map%interactions(i)%shiftedFrc = shiftedFrc
             MnM_Map%interactions(i)%rCut = thisRcut
             MnM_Map%interactions(i)%rCutWasSet = .true.
          enddo
       end if
  end if

  end subroutine setMnMDefaultCutoff

  subroutine copyAllData(v1, v2)
    type(MnMinteractionMap), pointer  :: v1
    type(MnMinteractionMap), pointer  :: v2
    integer :: i, j

    do i = 1, v1%interactionCount
       v2%interactions(i) = v1%interactions(i)
    enddo

    v2%interactionCount = v1%interactionCount
    return
  end subroutine copyAllData

  subroutine addInteraction(myInteraction)
    type(MNMtype) :: myInteraction
    type(MnMinteraction) :: nt
    integer :: id
 
    nt%interaction_type = myInteraction%MNMInteractionType
    nt%metal_atid = &
        getFirstMatchingElement(atypes, "c_ident", myInteraction%metal_atid)
    nt%nonmetal_atid = &
        getFirstMatchingElement(atypes, "c_ident", myInteraction%nonmetal_atid)
       
    select case (nt%interaction_type)
    case (MNM_LENNARDJONES)
       nt%sigma = myInteraction%sigma
       nt%epsilon = myInteraction%epsilon
    case(MNM_REPULSIVEMORSE, MNM_SHIFTEDMORSE)
       nt%R0 = myInteraction%R0
       nt%D0 = myInteraction%D0
       nt%beta0 = myInteraction%beta0
    case(MNM_MAW)
       nt%R0 = myInteraction%R0
       nt%D0 = myInteraction%D0
       nt%beta0 = myInteraction%beta0
       nt%ca1 = myInteraction%ca1
       nt%cb1 = myInteraction%cb1
    case default
       call handleError("MNM", "Unknown Interaction type")
    end select
    
    if (.not. associated(MnM_Map)) then
       call ensureCapacityHelper(MnM_Map, 1)
    else
       call ensureCapacityHelper(MnM_Map, MnM_Map%interactionCount + 1)
    end if
    
    MnM_Map%interactionCount = MnM_Map%interactionCount + 1
    id = MnM_Map%interactionCount
    MnM_Map%interactions(id) = nt
  end subroutine addInteraction

  subroutine ensureCapacityHelper(this, minCapacity)
    type(MnMinteractionMap), pointer :: this, that
    integer, intent(in) :: minCapacity
    integer :: oldCapacity 
    integer :: newCapacity
    logical :: resizeFlag 

    resizeFlag = .false.

    !  first time: allocate a new vector with default size

    if (.not. associated(this)) then
       this => MnMinitialize(minCapacity, 0)
    endif

    oldCapacity = size(this%interactions)

    if (minCapacity > oldCapacity) then
       if (this%capacityIncrement .gt. 0) then
          newCapacity = oldCapacity + this%capacityIncrement
       else
          newCapacity = oldCapacity * 2
       endif
       if (newCapacity .lt. minCapacity) then
          newCapacity = minCapacity
       endif
       resizeFlag = .true.
    else
       newCapacity = oldCapacity
    endif

    if (resizeFlag) then
       that => MnMinitialize(newCapacity, this%capacityIncrement)
       call copyAllData(this, that)
       this => MnMdestroy(this)
       this => that
    endif
  end subroutine ensureCapacityHelper

  function MnMinitialize(cap, capinc) result(this)
    integer, intent(in) :: cap, capinc
    integer :: error
    type(MnMinteractionMap), pointer :: this 

    nullify(this)

    if (cap < 0) then
       write(*,*) 'Bogus Capacity:', cap
       return
    endif
    allocate(this,stat=error)
    if ( error /= 0 ) then
       write(*,*) 'Could not allocate MnMinteractionMap!'
       return
    end if

    this%initialCapacity = cap
    this%capacityIncrement = capinc

    allocate(this%interactions(this%initialCapacity), stat=error)
    if(error /= 0) write(*,*) 'Could not allocate MnMinteraction!'

  end function MnMinitialize

  subroutine createInteractionLookup(this) 
    type (MnMInteractionMap),pointer :: this
    integer :: biggestAtid, i, metal_atid, nonmetal_atid, error

    biggestAtid=-1
    do i = 1, this%interactionCount
       metal_atid = this%interactions(i)%metal_atid
       nonmetal_atid = this%interactions(i)%nonmetal_atid

       if (metal_atid .gt. biggestAtid) biggestAtid = metal_atid
       if (nonmetal_atid .gt. biggestAtid) biggestAtid = nonmetal_atid
    enddo
    
    allocate(MnMinteractionLookup(biggestAtid,biggestAtid), stat=error)
    if (error /= 0) write(*,*) 'Could not allocate MnMinteractionLookup'

    do i = 1, this%interactionCount
       metal_atid = this%interactions(i)%metal_atid
       nonmetal_atid = this%interactions(i)%nonmetal_atid
    
       MnMinteractionLookup(metal_atid, nonmetal_atid) = i
       MnMinteractionLookup(nonmetal_atid, metal_atid) = i
    enddo
  end subroutine createInteractionLookup
    
  function MnMdestroy(this) result(null_this)
    logical :: done
    type(MnMinteractionMap), pointer :: this 
    type(MnMinteractionMap), pointer :: null_this 

    if (.not. associated(this)) then
       null_this => null()
       return
    end if

    !! Walk down the list and deallocate each of the map's components
    if(associated(this%interactions)) then
       deallocate(this%interactions)
       this%interactions=>null()
    endif
    deallocate(this)
    this => null()
    null_this => null()
  end function MnMdestroy

  subroutine deleteInteractions()    
    MnM_Map => MnMdestroy(MnM_Map)
    return
  end subroutine deleteInteractions

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
end module MetalNonMetal
