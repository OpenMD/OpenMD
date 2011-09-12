!!
!! Copyright (c) 2007 The University of Notre Dame. All Rights Reserved.
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


!! Calculates Metal-Non Metal interactions.
!! @author Charles F. Vardeman II 
!! @version $Id$, $Date$, $Name: not supported by cvs2svn $, $Revision$


module MetalNonMetal
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
     integer       :: nRep
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


  subroutine do_mnm_pair(atid1, atid2, D, Rij, R2, Rcut, Sw, vdwMult, &
       Vpair, Fpair, Pot, A1, A2, f1, t1, t2)
    integer, intent(in) ::  atid1, atid2
    integer :: ljt1, ljt2
    real( kind = dp ), intent(in) :: rij, r2, rcut, vdwMult
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3) :: f1 
    real (kind=dp), intent(inout), dimension(9) :: A1, A2
    real (kind=dp), intent(inout), dimension(3) :: t1, t2
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair

    integer :: interaction_id
    integer :: interaction_type

    if(.not.haveInteractionLookup)  then
      call createInteractionLookup(MnM_MAP)
      haveInteractionLookup =.true.
    end if
    
    interaction_id = MnMinteractionLookup(atid1, atid2)
    interaction_type = MnM_Map%interactions(interaction_id)%interaction_type
    
    select case (interaction_type)    
    case (MNM_LENNARDJONES)
       call calc_mnm_lennardjones(D, Rij, R2, Rcut, Sw, &
            vdwMult, Vpair, Fpair, Pot, f1, interaction_id)
    case (MNM_REPULSIVEPOWER)
       call calc_mnm_repulsivepower(D, Rij, R2, Rcut, Sw, &
            vdwMult, Vpair, Fpair, Pot, f1, interaction_id)
    case(MNM_REPULSIVEMORSE, MNM_SHIFTEDMORSE)
       call calc_mnm_morse(D, Rij, R2, Rcut, Sw, vdwMult, &
            Vpair, Fpair, Pot, f1, interaction_id, interaction_type)
    case(MNM_MAW)
       call calc_mnm_maw(atid1, atid2, D, Rij, R2, Rcut, Sw, vdwMult, &
            Vpair, Fpair, Pot, A1, A2, f1, t1, t2, interaction_id)
    case default
    call handleError("MetalNonMetal","Unknown interaction type")      
    end select

  end subroutine do_mnm_pair

  subroutine calc_mnm_lennardjones(D, Rij, R2, Rcut, Sw, &
       vdwMult,Vpair, Fpair, Pot, f1, interaction_id)
    
    real( kind = dp ), intent(in) :: rij, r2, rcut, vdwMult
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3) :: f1
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
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
  end subroutine calc_mnm_lennardjones


  subroutine calc_mnm_repulsivepower(D, Rij, R2, Rcut, Sw, &
       vdwMult,Vpair, Fpair, Pot, f1, interaction_id)
    
    real( kind = dp ), intent(in) :: rij, r2, rcut, vdwMult
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3) :: f1
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    integer, intent(in) :: interaction_id

    ! local Variables
    real( kind = dp ) :: drdx, drdy, drdz
    real( kind = dp ) :: fx, fy, fz
    real( kind = dp ) :: myPot, myPotC, myDeriv, myDerivC, ros, rcos
    real( kind = dp ) :: pot_temp, dudr
    real( kind = dp ) :: sigmai
    real( kind = dp ) :: epsilon
    logical :: isSoftCore, shiftedPot, shiftedFrc
    integer :: id1, id2, localError, n

    sigmai     = 1.0_dp / MnM_Map%interactions(interaction_id)%sigma
    epsilon    = MnM_Map%interactions(interaction_id)%epsilon
    n          = MnM_Map%interactions(interaction_id)%nRep
    shiftedPot = MnM_Map%interactions(interaction_id)%shiftedPot
    shiftedFrc = MnM_Map%interactions(interaction_id)%shiftedFrc

    ros = rij * sigmai

    call getNRepulsionFunc(ros, n, myPot, myDeriv)
    
    if (shiftedPot) then
       rcos = rcut * sigmai
       call getNRepulsionFunc(rcos, n, myPotC, myDerivC) 
       myDerivC = 0.0_dp
    elseif (shiftedFrc) then
       rcos = rcut * sigmai
       call getNRepulsionFunc(rcos, n, myPotC, myDerivC)
       myPotC = myPotC + myDerivC * (rij - rcut) * sigmai
    else
       myPotC = 0.0_dp
       myDerivC = 0.0_dp
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
  end subroutine calc_mnm_repulsivepower


  subroutine calc_mnm_morse(D, Rij, R2, Rcut, Sw, vdwMult, &
       Vpair, Fpair, Pot, f1, interaction_id, interaction_type)
    real( kind = dp ), intent(in) :: rij, r2, rcut, vdwMult
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3) :: f1
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
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

    pot_temp = vdwMult * (myPot - myPotC)
    vpair = vpair + pot_temp
    dudr = sw * vdwMult * (myDeriv - myDerivC)

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
  end subroutine calc_mnm_morse
  
  subroutine calc_mnm_maw(atid1, atid2, D, Rij, R2, Rcut, Sw, vdwMult, &
       Vpair, Fpair, Pot, A1, A2, f1, t1, t2, interaction_id)
    real( kind = dp ), intent(in) :: rij, r2, rcut, vdwMult
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3) :: f1  
    real (kind=dp),intent(in), dimension(9) :: A1, A2
    real (kind=dp),intent(inout), dimension(3) :: t1, t2

    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair

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
    real(kind = dp) :: ta1, b1, s
    real(kind = dp) :: da1dx, da1dy, da1dz, da1dux, da1duy, da1duz
    real(kind = dp) :: db1dx, db1dy, db1dz, db1dux, db1duy, db1duz
    real(kind = dp) :: fx, fy, fz, tx, ty, tz, fxl, fyl, fzl
!    real(kind = dp), parameter :: st = sqrt(3.0_dp)
    real(kind = dp), parameter :: st = 1.732050807568877
    integer :: atid1, atid2, id1, id2
    logical :: shiftedPot, shiftedFrc
   
    if (atid2.eq.MnM_Map%interactions(interaction_id)%metal_atid) then
       ! rotate the inter-particle separation into the two different 
       ! body-fixed coordinate systems:
       
       x = A1(1)*d(1) + A1(2)*d(2) + A1(3)*d(3)
       y = A1(4)*d(1) + A1(5)*d(2) + A1(6)*d(3)
       z = A1(7)*d(1) + A1(8)*d(2) + A1(9)*d(3)
    else
       ! negative sign because this is the vector from j to i:
       
       x = -(A2(1)*d(1) + A2(2)*d(2) + A2(3)*d(3))
       y = -(A2(4)*d(1) + A2(5)*d(2) + A2(6)*d(3))
       z = -(A2(7)*d(1) + A2(8)*d(2) + A2(9)*d(3))
    endif
      
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
    !********************** old form*************************
    ! s = 1 / (4 pi)
    ! ta1 = (s - pz)^2 = (1 - sqrt(3)*cos(theta))^2 / (4 pi)
    ! b1 = px^2 = 3 * (sin(theta)*cos(phi))^2 / (4 pi)   
    !********************** old form*************************
    ! we'll leave out the 4 pi for now
    
    ! new functional form just using the p orbitals.
    ! Vmorse(r)*[a*p_x + b p_z + (1-a-b)]
    ! which is 
    ! Vmorse(r)*[a sin^2(theta) cos^2(phi) + b cos(theta) + (1-a-b)]
    ! Vmorse(r)*[a*x2/r2 + b*z/r + (1-a-b)]



    s = 1.0_dp
!    ta1 = (1.0_dp - st * z / rij )**2
!    b1 = 3.0_dp * x2 / r2

!    Vang = s + ca1 * ta1 + cb1 * b1
 
    Vang = ca1 * x2/r2 + cb1 * z/rij + (0.8_dp-ca1/3.0_dp)

    pot_temp = vdwMult * Vmorse*Vang 
         
    vpair = vpair + pot_temp
    pot = pot + pot_temp*sw
    
    dVmorsedr = 2.0_dp*D0*beta0*(expfnc0 - expfnc02)

    dVmorsedx = dVmorsedr * drdx
    dVmorsedy = dVmorsedr * drdy
    dVmorsedz = dVmorsedr * drdz
    
    !da1dx = 2.0_dp * st * x * z / r3 - 6.0_dp * x * z2 / r4
    !da1dy = 2.0_dp * st * y * z / r3 - 6.0_dp * y * z2 / r4
    !da1dz = 2.0_dp * st * (x2 + y2) * (st * z - rij ) / r4

    !db1dx = 6.0_dp * x * (1.0_dp - x2 / r2) / r2
    !db1dy = -6.0_dp * x2 * y / r4
    !db1dz = -6.0_dp * x2 * z / r4

    !dVangdx = ca1 * da1dx + cb1 * db1dx
    !dVangdy = ca1 * da1dy + cb1 * db1dy
    !dVangdz = ca1 * da1dz + cb1 * db1dz
    
    dVangdx = -2.0*ca1*x2*x/r4 + 2.0*ca1*x/r2 - cb1*x*z/r3
    dVangdy = -2.0*ca1*x2*y/r4                - cb1*z*y/r3
    dVangdz = -2.0*ca1*x2*z/r4 + cb1/rij      - cb1*z2 /r3

    ! chain rule to put these back on x, y, z
    dvdx = Vang * dVmorsedx + Vmorse * dVangdx
    dvdy = Vang * dVmorsedy + Vmorse * dVangdy
    dvdz = Vang * dVmorsedz + Vmorse * dVangdz
    
    ! Torques for Vang using method of Price:
    ! S. L. Price, A. J. Stone, and M. Alderton, Mol. Phys. 52, 987 (1984).

    !da1dux =   6.0_dp * y * z / r2 - 2.0_dp * st * y / rij
    !da1duy =  -6.0_dp * x * z / r2 + 2.0_dp * st * x / rij
    !da1duz =   0.0_dp

    !db1dux =   0.0_dp
    !db1duy =   6.0_dp * x * z / r2
    !db1duz =  -6.0_dp * y * x / r2

    !dVangdux = ca1 * da1dux + cb1 * db1dux
    !dVangduy = ca1 * da1duy + cb1 * db1duy
    !dVangduz = ca1 * da1duz + cb1 * db1duz
    
    dVangdux = cb1*y/rij
    dVangduy = 2.0*ca1*x*z/r2 - cb1*x/rij
    dVangduz = -2.0*ca1*y*x/r2

    ! do the torques first since they are easy:
    ! remember that these are still in the body fixed axes    
    
    tx = vdwMult * Vmorse * dVangdux * sw
    ty = vdwMult * Vmorse * dVangduy * sw
    tz = vdwMult * Vmorse * dVangduz * sw

    ! go back to lab frame using transpose of rotation matrix:

    if (atid2.eq.MnM_Map%interactions(interaction_id)%metal_atid) then
       t1(1) = t1(1) + a1(1)*tx + a1(4)*ty + a1(7)*tz
       t1(2) = t1(2) + a1(2)*tx + a1(5)*ty + a1(8)*tz
       t1(3) = t1(3) + a1(3)*tx + a1(6)*ty + a1(9)*tz
    else
       t2(1) = t2(1) + a2(1)*tx + a2(4)*ty + a2(7)*tz
       t2(2) = t2(2) + a2(2)*tx + a2(5)*ty + a2(8)*tz
       t2(3) = t2(3) + a2(3)*tx + a2(6)*ty + a2(9)*tz
    endif

    ! Now, on to the forces (still in body frame of water)

    fx = vdwMult * dvdx * sw
    fy = vdwMult * dvdy * sw
    fz = vdwMult * dvdz * sw

    ! rotate the terms back into the lab frame:

    if (atid2.eq.MnM_Map%interactions(interaction_id)%metal_atid) then
       fxl = a1(1)*fx + a1(4)*fy + a1(7)*fz
       fyl = a1(2)*fx + a1(5)*fy + a1(8)*fz
       fzl = a1(3)*fx + a1(6)*fy + a1(9)*fz
    else
       ! negative sign because this is the vector from j to i:
       fxl = -(a2(1)*fx + a2(4)*fy + a2(7)*fz)
       fyl = -(a2(2)*fx + a2(5)*fy + a2(8)*fz)
       fzl = -(a2(3)*fx + a2(6)*fy + a2(9)*fz)
    endif

    f1(1) = f1(1) + fxl
    f1(2) = f1(2) + fyl
    f1(3) = f1(3) + fzl

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
    case (MNM_REPULSIVEPOWER)
       nt%sigma = myInteraction%sigma
       nt%epsilon = myInteraction%epsilon
       nt%nRep = myInteraction%nRep
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

  subroutine getNRepulsionFunc(r, n, myPot, myDeriv)
    
    real(kind=dp), intent(in) :: r
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: myPot, myDeriv
    real(kind=dp) :: ri, rin, rin1
    
    ri = 1.0_DP / r   

    rin = ri**n
    rin1 = rin * ri

    myPot = rin
    myDeriv = -n * rin1
    
    return
  end subroutine getNRepulsionFunc
end module MetalNonMetal
