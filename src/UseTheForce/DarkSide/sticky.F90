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

!! This Module Calculates forces due to SSD potential and VDW interactions
!! [Chandra and Ichiye, J. Chem. Phys. 111, 2701 (1999)].

!! This module contains the Public procedures:


!! Corresponds to the force field defined in ssd_FF.cpp 
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @author Christopher Fennell
!! @author J. Daniel Gezelter
!! @version $Id$, $Date$, $Name: not supported by cvs2svn $, $Revision$

module sticky

  use force_globals
  use definitions
  use atype_module
  use vector_class
  use simulation
  use status
  use interpolation
  implicit none

  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"

  public :: newStickyType
  public :: do_sticky_pair
  public :: destroyStickyTypes
  public :: do_sticky_power_pair
  public :: getStickyCut
  public :: getStickyPowerCut

  type :: StickyList
     integer :: c_ident
     real( kind = dp ) :: w0 = 0.0_dp
     real( kind = dp ) :: v0 = 0.0_dp
     real( kind = dp ) :: v0p = 0.0_dp
     real( kind = dp ) :: rl = 0.0_dp
     real( kind = dp ) :: ru = 0.0_dp
     real( kind = dp ) :: rlp = 0.0_dp
     real( kind = dp ) :: rup = 0.0_dp
     real( kind = dp ) :: rbig = 0.0_dp
     type(cubicSpline) :: stickySpline
     type(cubicSpline) :: stickySplineP
  end type StickyList

  type(StickyList), dimension(:),allocatable :: StickyMap
  logical, save :: hasStickyMap = .false.

contains

  subroutine newStickyType(c_ident, w0, v0, v0p, rl, ru, rlp, rup, isError)

    integer, intent(in) :: c_ident
    integer, intent(inout) :: isError
    real( kind = dp ), intent(in) :: w0, v0, v0p
    real( kind = dp ), intent(in) :: rl, ru
    real( kind = dp ), intent(in) :: rlp, rup
    real( kind = dp ), dimension(2) :: rCubVals, sCubVals, rpCubVals, spCubVals
    integer :: nATypes, myATID


    isError = 0
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    !! Be simple-minded and assume that we need a StickyMap that
    !! is the same size as the total number of atom types

    if (.not.allocated(StickyMap)) then

       nAtypes = getSize(atypes)

       if (nAtypes == 0) then
          isError = -1
          return
       end if

       if (.not. allocated(StickyMap)) then
          allocate(StickyMap(nAtypes))
       endif

    end if

    if (myATID .gt. size(StickyMap)) then
       isError = -1
       return
    endif

    ! set the values for StickyMap for this atom type:

    StickyMap(myATID)%c_ident = c_ident

    ! we could pass all 5 parameters if we felt like it...

    StickyMap(myATID)%w0 = w0
    StickyMap(myATID)%v0 = v0
    StickyMap(myATID)%v0p = v0p
    StickyMap(myATID)%rl = rl
    StickyMap(myATID)%ru = ru
    StickyMap(myATID)%rlp = rlp
    StickyMap(myATID)%rup = rup

    if (StickyMap(myATID)%ru .gt. StickyMap(myATID)%rup) then
       StickyMap(myATID)%rbig = StickyMap(myATID)%ru
    else
       StickyMap(myATID)%rbig = StickyMap(myATID)%rup
    endif

    ! build the 2 cubic splines for the sticky switching functions

    rCubVals(1) = rl
    rCubVals(2) = ru
    sCubVals(1) = 1.0_dp
    sCubVals(2) = 0.0_dp      
    call newSpline(StickyMap(myATID)%stickySpline, rCubVals, sCubVals, .true.)
    rpCubVals(1) = rlp
    rpCubVals(2) = rup
    spCubVals(1) = 1.0_dp
    spCubVals(2) = 0.0_dp      
    call newSpline(StickyMap(myATID)%stickySplineP,rpCubVals,spCubVals,.true.)

    hasStickyMap = .true.

    return
  end subroutine newStickyType

  function getStickyCut(atomID) result(cutValue)
    integer, intent(in) :: atomID
    real(kind=dp) :: cutValue

    cutValue = StickyMap(atomID)%rbig
  end function getStickyCut

  function getStickyPowerCut(atomID) result(cutValue)
    integer, intent(in) :: atomID
    real(kind=dp) :: cutValue

    cutValue = StickyMap(atomID)%rbig
  end function getStickyPowerCut

  subroutine do_sticky_pair(me1, me2, d, rij, r2, sw, vpair, fpair, &
       pot, A1, A2, f1, t1, t2)

    !! This routine does only the sticky portion of the SSD potential
    !! [Chandra and Ichiye, J. Chem. Phys. 111, 2701 (1999)].
    !! The Lennard-Jones and dipolar interaction must be handled separately.

    !! We assume that the rotation matrices have already been calculated
    !! and placed in the A array.

    !! i and j are pointers to the two SSD atoms

    integer, intent(in) :: me1, me2
    real (kind=dp), intent(inout) :: rij, r2
    real (kind=dp), dimension(3), intent(in) :: d
    real (kind=dp), dimension(3), intent(inout) :: fpair
    real (kind=dp) :: pot, vpair, sw
    real (kind=dp), dimension(9) :: A1, A2
    real (kind=dp), dimension(3) :: f1
    real (kind=dp), dimension(3) :: t1, t2

    real (kind=dp) :: xi, yi, zi, xj, yj, zj, xi2, yi2, zi2, xj2, yj2, zj2
    real (kind=dp) :: r3, r5, r6, s, sp, dsdr, dspdr
    real (kind=dp) :: wi, wj, w, wip, wjp, wp
    real (kind=dp) :: dwidx, dwidy, dwidz, dwjdx, dwjdy, dwjdz
    real (kind=dp) :: dwipdx, dwipdy, dwipdz, dwjpdx, dwjpdy, dwjpdz
    real (kind=dp) :: dwidux, dwiduy, dwiduz, dwjdux, dwjduy, dwjduz
    real (kind=dp) :: dwipdux, dwipduy, dwipduz, dwjpdux, dwjpduy, dwjpduz
    real (kind=dp) :: zif, zis, zjf, zjs, uglyi, uglyj
    real (kind=dp) :: drdx, drdy, drdz
    real (kind=dp) :: txi, tyi, tzi, txj, tyj, tzj
    real (kind=dp) :: fxii, fyii, fzii, fxjj, fyjj, fzjj
    real (kind=dp) :: fxij, fyij, fzij, fxji, fyji, fzji       
    real (kind=dp) :: fxradial, fyradial, fzradial
    real (kind=dp) :: rijtest, rjitest
    real (kind=dp) :: radcomxi, radcomyi, radcomzi
    real (kind=dp) :: radcomxj, radcomyj, radcomzj
    integer :: id1, id2

    real (kind=dp) :: w0, v0, v0p, rl, ru, rlp, rup, rbig, dx

    if (me1.eq.me2) then
       w0  = StickyMap(me1)%w0 
       v0  = StickyMap(me1)%v0 
       v0p = StickyMap(me1)%v0p
       rl  = StickyMap(me1)%rl 
       ru  = StickyMap(me1)%ru 
       rlp = StickyMap(me1)%rlp
       rup = StickyMap(me1)%rup
       rbig = StickyMap(me1)%rbig
    else
       ! This is silly, but if you want 2 sticky types in your 
       ! simulation, we'll let you do it with the Lorentz-
       ! Berthelot mixing rules.
       ! (Warning: you'll be SLLLLLLLLLLLLLLLOOOOOOOOOOWWWWWWWWWWW)
       rl   = 0.5_dp * ( StickyMap(me1)%rl + StickyMap(me2)%rl )
       ru   = 0.5_dp * ( StickyMap(me1)%ru + StickyMap(me2)%ru )
       rlp  = 0.5_dp * ( StickyMap(me1)%rlp + StickyMap(me2)%rlp )
       rup  = 0.5_dp * ( StickyMap(me1)%rup + StickyMap(me2)%rup )
       rbig = max(ru, rup)
       w0  = sqrt( StickyMap(me1)%w0   * StickyMap(me2)%w0  )
       v0  = sqrt( StickyMap(me1)%v0   * StickyMap(me2)%v0  )
       v0p = sqrt( StickyMap(me1)%v0p  * StickyMap(me2)%v0p )
    endif

    if ( rij .LE. rbig ) then

       r3 = r2*rij
       r5 = r3*r2

       drdx = d(1) / rij
       drdy = d(2) / rij
       drdz = d(3) / rij

       ! rotate the inter-particle separation into the two different 
       ! body-fixed coordinate systems:

       xi = A1(1)*d(1) + A1(2)*d(2) + A1(3)*d(3)
       yi = A1(4)*d(1) + A1(5)*d(2) + A1(6)*d(3)
       zi = A1(7)*d(1) + A1(8)*d(2) + A1(9)*d(3)

       ! negative sign because this is the vector from j to i:

       xj = -(A2(1)*d(1) + A2(2)*d(2) + A2(3)*d(3))
       yj = -(A2(4)*d(1) + A2(5)*d(2) + A2(6)*d(3))
       zj = -(A2(7)*d(1) + A2(8)*d(2) + A2(9)*d(3))

       xi2 = xi*xi
       yi2 = yi*yi
       zi2 = zi*zi

       xj2 = xj*xj
       yj2 = yj*yj
       zj2 = zj*zj

       ! calculate the switching info. from the splines
       if (me1.eq.me2) then
          s = 0.0_dp
          dsdr = 0.0_dp
          sp = 0.0_dp
          dspdr = 0.0_dp
          
          if (rij.lt.ru) then
             if (rij.lt.rl) then
                s = 1.0_dp
                dsdr = 0.0_dp
             else         
                ! we are in the switching region 
                dx = rij - rl
                s = StickyMap(me1)%stickySpline%y(1) + &
                     dx*(dx*(StickyMap(me1)%stickySpline%c(1) + &
                     dx*StickyMap(me1)%stickySpline%d(1)))
                dsdr = dx*(2.0_dp * StickyMap(me1)%stickySpline%c(1) + &
                     3.0_dp * dx * StickyMap(me1)%stickySpline%d(1))
             endif
          endif
          if (rij.lt.rup) then
             if (rij.lt.rlp) then
                sp = 1.0_dp
                dspdr = 0.0_dp
             else
                ! we are in the switching region 
                dx = rij - rlp
                sp = StickyMap(me1)%stickySplineP%y(1) + &
                     dx*(dx*(StickyMap(me1)%stickySplineP%c(1) + &
                     dx*StickyMap(me1)%stickySplineP%d(1)))
                dspdr = dx*(2.0_dp * StickyMap(me1)%stickySplineP%c(1) + &
                     3.0_dp * dx * StickyMap(me1)%stickySplineP%d(1))
             endif
          endif
       else
          ! calculate the switching function explicitly rather than from 
          ! the splines with mixed sticky maps
          call calc_sw_fnc(rij, rl, ru, rlp, rup, s, sp, dsdr, dspdr)
       endif

       wi = 2.0_dp*(xi2-yi2)*zi / r3
       wj = 2.0_dp*(xj2-yj2)*zj / r3
       w = wi+wj

       zif = zi/rij - 0.6_dp
       zis = zi/rij + 0.8_dp

       zjf = zj/rij - 0.6_dp
       zjs = zj/rij + 0.8_dp

       wip = zif*zif*zis*zis - w0
       wjp = zjf*zjf*zjs*zjs - w0
       wp = wip + wjp

       vpair = vpair + 0.5_dp*(v0*s*w + v0p*sp*wp)

       pot = pot + 0.5_dp*(v0*s*w + v0p*sp*wp)*sw

       dwidx =   4.0_dp*xi*zi/r3  - 6.0_dp*xi*zi*(xi2-yi2)/r5
       dwidy = - 4.0_dp*yi*zi/r3  - 6.0_dp*yi*zi*(xi2-yi2)/r5
       dwidz =   2.0_dp*(xi2-yi2)/r3  - 6.0_dp*zi2*(xi2-yi2)/r5

       dwjdx =   4.0_dp*xj*zj/r3  - 6.0_dp*xj*zj*(xj2-yj2)/r5
       dwjdy = - 4.0_dp*yj*zj/r3  - 6.0_dp*yj*zj*(xj2-yj2)/r5
       dwjdz =   2.0_dp*(xj2-yj2)/r3  - 6.0_dp*zj2*(xj2-yj2)/r5

       uglyi = zif*zif*zis + zif*zis*zis
       uglyj = zjf*zjf*zjs + zjf*zjs*zjs

       dwipdx = -2.0_dp*xi*zi*uglyi/r3
       dwipdy = -2.0_dp*yi*zi*uglyi/r3
       dwipdz = 2.0_dp*(1.0_dp/rij - zi2/r3)*uglyi

       dwjpdx = -2.0_dp*xj*zj*uglyj/r3
       dwjpdy = -2.0_dp*yj*zj*uglyj/r3
       dwjpdz = 2.0_dp*(1.0_dp/rij - zj2/r3)*uglyj

       dwidux = 4.0_dp*(yi*zi2 + 0.5_dp*yi*(xi2-yi2))/r3
       dwiduy = 4.0_dp*(xi*zi2 - 0.5_dp*xi*(xi2-yi2))/r3
       dwiduz = - 8.0_dp*xi*yi*zi/r3

       dwjdux = 4.0_dp*(yj*zj2 + 0.5_dp*yj*(xj2-yj2))/r3
       dwjduy = 4.0_dp*(xj*zj2 - 0.5_dp*xj*(xj2-yj2))/r3
       dwjduz = - 8.0_dp*xj*yj*zj/r3

       dwipdux =  2.0_dp*yi*uglyi/rij
       dwipduy = -2.0_dp*xi*uglyi/rij
       dwipduz =  0.0_dp

       dwjpdux =  2.0_dp*yj*uglyj/rij
       dwjpduy = -2.0_dp*xj*uglyj/rij
       dwjpduz =  0.0_dp

       ! do the torques first since they are easy:
       ! remember that these are still in the body fixed axes

       txi = 0.5_dp*(v0*s*dwidux + v0p*sp*dwipdux)*sw
       tyi = 0.5_dp*(v0*s*dwiduy + v0p*sp*dwipduy)*sw
       tzi = 0.5_dp*(v0*s*dwiduz + v0p*sp*dwipduz)*sw

       txj = 0.5_dp*(v0*s*dwjdux + v0p*sp*dwjpdux)*sw
       tyj = 0.5_dp*(v0*s*dwjduy + v0p*sp*dwjpduy)*sw
       tzj = 0.5_dp*(v0*s*dwjduz + v0p*sp*dwjpduz)*sw

       ! go back to lab frame using transpose of rotation matrix:

       t1(1) = t1(1) + a1(1)*txi + a1(4)*tyi + a1(7)*tzi
       t1(2) = t1(2) + a1(2)*txi + a1(5)*tyi + a1(8)*tzi
       t1(3) = t1(3) + a1(3)*txi + a1(6)*tyi + a1(9)*tzi

       t2(1) = t2(1) + a2(1)*txj + a2(4)*tyj + a2(7)*tzj
       t2(2) = t2(2) + a2(2)*txj + a2(5)*tyj + a2(8)*tzj
       t2(3) = t2(3) + a2(3)*txj + a2(6)*tyj + a2(9)*tzj

       ! Now, on to the forces:

       ! first rotate the i terms back into the lab frame:

       radcomxi = (v0*s*dwidx+v0p*sp*dwipdx)*sw
       radcomyi = (v0*s*dwidy+v0p*sp*dwipdy)*sw
       radcomzi = (v0*s*dwidz+v0p*sp*dwipdz)*sw

       radcomxj = (v0*s*dwjdx+v0p*sp*dwjpdx)*sw
       radcomyj = (v0*s*dwjdy+v0p*sp*dwjpdy)*sw
       radcomzj = (v0*s*dwjdz+v0p*sp*dwjpdz)*sw

       fxii = a1(1)*(radcomxi) + a1(4)*(radcomyi) + a1(7)*(radcomzi)
       fyii = a1(2)*(radcomxi) + a1(5)*(radcomyi) + a1(8)*(radcomzi)
       fzii = a1(3)*(radcomxi) + a1(6)*(radcomyi) + a1(9)*(radcomzi)

       fxjj = a2(1)*(radcomxj) + a2(4)*(radcomyj) + a2(7)*(radcomzj)
       fyjj = a2(2)*(radcomxj) + a2(5)*(radcomyj) + a2(8)*(radcomzj)
       fzjj = a2(3)*(radcomxj) + a2(6)*(radcomyj) + a2(9)*(radcomzj)

       fxij = -fxii
       fyij = -fyii
       fzij = -fzii

       fxji = -fxjj
       fyji = -fyjj
       fzji = -fzjj

       ! now assemble these with the radial-only terms:

       fxradial = 0.5_dp*(v0*dsdr*drdx*w + v0p*dspdr*drdx*wp + fxii + fxji)
       fyradial = 0.5_dp*(v0*dsdr*drdy*w + v0p*dspdr*drdy*wp + fyii + fyji)
       fzradial = 0.5_dp*(v0*dsdr*drdz*w + v0p*dspdr*drdz*wp + fzii + fzji)

       f1(1) = f1(1) + fxradial
       f1(2) = f1(2) + fyradial
       f1(3) = f1(3) + fzradial

    endif
  end subroutine do_sticky_pair

  !! calculates the switching functions and their derivatives for a given
  subroutine calc_sw_fnc(r, rl, ru, rlp, rup, s, sp, dsdr, dspdr)

    real (kind=dp), intent(in) :: r, rl, ru, rlp, rup
    real (kind=dp), intent(inout) :: s, sp, dsdr, dspdr

    ! distances must be in angstroms
    s = 0.0_dp
    dsdr = 0.0_dp
    sp = 0.0_dp
    dspdr = 0.0_dp
    
    if (r.lt.ru) then
       if (r.lt.rl) then
          s = 1.0_dp
          dsdr = 0.0_dp
       else
          s = ((ru + 2.0_dp*r - 3.0_dp*rl) * (ru-r)**2) / &
               ((ru - rl)**3)
          dsdr = 6.0_dp*(r-ru)*(r-rl)/((ru - rl)**3)
       endif
    endif

    if (r.lt.rup) then
       if (r.lt.rlp) then
          sp = 1.0_dp       
          dspdr = 0.0_dp
       else
          sp = ((rup + 2.0_dp*r - 3.0_dp*rlp) * (rup-r)**2) / &
               ((rup - rlp)**3)
          dspdr = 6.0_dp*(r-rup)*(r-rlp)/((rup - rlp)**3)       
       endif
    endif

    return
  end subroutine calc_sw_fnc

  subroutine destroyStickyTypes()   
    if(allocated(StickyMap)) deallocate(StickyMap)
  end subroutine destroyStickyTypes
  
  subroutine do_sticky_power_pair(me1, me2, d, rij, r2, sw, vpair, fpair, &
       pot, A1, A2, f1, t1, t2)
    
    !! i and j are pointers to the two SSD atoms
    
    real (kind=dp), intent(inout) :: rij, r2
    real (kind=dp), dimension(3), intent(in) :: d
    real (kind=dp), dimension(3), intent(inout) :: fpair
    real (kind=dp) :: pot, vpair, sw
    real (kind=dp), dimension(9) :: A1, A2
    real (kind=dp), dimension(3) :: f1
    real (kind=dp), dimension(3) :: t1, t2


    real (kind=dp) :: xi, yi, zi, xj, yj, zj, xi2, yi2, zi2, xj2, yj2, zj2
    real (kind=dp) :: xihat, yihat, zihat, xjhat, yjhat, zjhat
    real (kind=dp) :: rI, rI2, rI3, rI4, rI5, rI6, rI7, s, sp, dsdr, dspdr
    real (kind=dp) :: wi, wj, w, wi2, wj2, eScale, v0scale
    real (kind=dp) :: dwidx, dwidy, dwidz, dwjdx, dwjdy, dwjdz
    real (kind=dp) :: dwidux, dwiduy, dwiduz, dwjdux, dwjduy, dwjduz
    real (kind=dp) :: drdx, drdy, drdz
    real (kind=dp) :: txi, tyi, tzi, txj, tyj, tzj
    real (kind=dp) :: fxii, fyii, fzii, fxjj, fyjj, fzjj
    real (kind=dp) :: fxij, fyij, fzij, fxji, fyji, fzji       
    real (kind=dp) :: fxradial, fyradial, fzradial
    real (kind=dp) :: rijtest, rjitest
    real (kind=dp) :: radcomxi, radcomyi, radcomzi
    real (kind=dp) :: radcomxj, radcomyj, radcomzj
    integer :: id1, id2
    integer :: me1, me2
    real (kind=dp) :: w0, v0, v0p, rl, ru, rlp, rup, rbig
    real (kind=dp) :: zi3, zi4, zi5, zj3, zj4, zj5
    real (kind=dp) :: frac1, frac2
          
    if (.not.allocated(StickyMap)) then
       call handleError("sticky", "no StickyMap was present before first call of do_sticky_power_pair!")
       return
    end if
 
    if (me1.eq.me2) then
       w0  = StickyMap(me1)%w0 
       v0  = StickyMap(me1)%v0 
       v0p = StickyMap(me1)%v0p
       rl  = StickyMap(me1)%rl 
       ru  = StickyMap(me1)%ru 
       rlp = StickyMap(me1)%rlp
       rup = StickyMap(me1)%rup
       rbig = StickyMap(me1)%rbig
    else
       ! This is silly, but if you want 2 sticky types in your 
       ! simulation, we'll let you do it with the Lorentz-
       ! Berthelot mixing rules.
       ! (Warning: you'll be SLLLLLLLLLLLLLLLOOOOOOOOOOWWWWWWWWWWW)
       rl   = 0.5_dp * ( StickyMap(me1)%rl + StickyMap(me2)%rl )
       ru   = 0.5_dp * ( StickyMap(me1)%ru + StickyMap(me2)%ru )
       rlp  = 0.5_dp * ( StickyMap(me1)%rlp + StickyMap(me2)%rlp )
       rup  = 0.5_dp * ( StickyMap(me1)%rup + StickyMap(me2)%rup )
       rbig = max(ru, rup)
       w0  = sqrt( StickyMap(me1)%w0   * StickyMap(me2)%w0  )
       v0  = sqrt( StickyMap(me1)%v0   * StickyMap(me2)%v0  )
       v0p = sqrt( StickyMap(me1)%v0p  * StickyMap(me2)%v0p )
    endif

    if ( rij .LE. rbig ) then

       rI = 1.0_dp/rij
       rI2 = rI*rI
       rI3 = rI2*rI
       rI4 = rI2*rI2
       rI5 = rI3*rI2
       rI6 = rI3*rI3
       rI7 = rI4*rI3
              
       drdx = d(1) * rI
       drdy = d(2) * rI
       drdz = d(3) * rI

       ! rotate the inter-particle separation into the two different 
       ! body-fixed coordinate systems:

       xi = A1(1)*d(1) + A1(2)*d(2) + A1(3)*d(3)
       yi = A1(4)*d(1) + A1(5)*d(2) + A1(6)*d(3)
       zi = A1(7)*d(1) + A1(8)*d(2) + A1(9)*d(3)

       ! negative sign because this is the vector from j to i:

       xj = -(A2(1)*d(1) + A2(2)*d(2) + A2(3)*d(3))
       yj = -(A2(4)*d(1) + A2(5)*d(2) + A2(6)*d(3))
       zj = -(A2(7)*d(1) + A2(8)*d(2) + A2(9)*d(3))

       xi2 = xi*xi
       yi2 = yi*yi
       zi2 = zi*zi
       zi3 = zi2*zi
       zi4 = zi2*zi2
       zi5 = zi3*zi2
       xihat = xi*rI
       yihat = yi*rI
       zihat = zi*rI
       
       xj2 = xj*xj
       yj2 = yj*yj
       zj2 = zj*zj
       zj3 = zj2*zj
       zj4 = zj2*zj2
       zj5 = zj3*zj2
       xjhat = xj*rI
       yjhat = yj*rI
       zjhat = zj*rI
       
       call calc_sw_fnc(rij, rl, ru, rlp, rup, s, sp, dsdr, dspdr)
           
       frac1 = 0.25_dp
       frac2 = 0.75_dp
       
       wi = 2.0_dp*(xi2-yi2)*zi*rI3
       wj = 2.0_dp*(xj2-yj2)*zj*rI3
       
       wi2 = wi*wi
       wj2 = wj*wj

       w = frac1*wi*wi2 + frac2*wi + frac1*wj*wj2 + frac2*wj + v0p

       vpair = vpair + 0.5_dp*(v0*s*w)
       
       pot = pot + 0.5_dp*(v0*s*w)*sw

       dwidx = ( 4.0_dp*xi*zi*rI3 - 6.0_dp*xi*zi*(xi2-yi2)*rI5 )
       dwidy = ( -4.0_dp*yi*zi*rI3 - 6.0_dp*yi*zi*(xi2-yi2)*rI5 )
       dwidz = ( 2.0_dp*(xi2-yi2)*rI3 - 6.0_dp*zi2*(xi2-yi2)*rI5 )
       
       dwidx = frac1*3.0_dp*wi2*dwidx + frac2*dwidx
       dwidy = frac1*3.0_dp*wi2*dwidy + frac2*dwidy
       dwidz = frac1*3.0_dp*wi2*dwidz + frac2*dwidz

       dwjdx = ( 4.0_dp*xj*zj*rI3  - 6.0_dp*xj*zj*(xj2-yj2)*rI5 )
       dwjdy = ( -4.0_dp*yj*zj*rI3  - 6.0_dp*yj*zj*(xj2-yj2)*rI5 )
       dwjdz = ( 2.0_dp*(xj2-yj2)*rI3  - 6.0_dp*zj2*(xj2-yj2)*rI5 )

       dwjdx = frac1*3.0_dp*wj2*dwjdx + frac2*dwjdx
       dwjdy = frac1*3.0_dp*wj2*dwjdy + frac2*dwjdy
       dwjdz = frac1*3.0_dp*wj2*dwjdz + frac2*dwjdz
       
       dwidux = ( 4.0_dp*(yi*zi2 + 0.5_dp*yi*(xi2-yi2))*rI3 )
       dwiduy = ( 4.0_dp*(xi*zi2 - 0.5_dp*xi*(xi2-yi2))*rI3 )
       dwiduz = ( -8.0_dp*xi*yi*zi*rI3 )

       dwidux = frac1*3.0_dp*wi2*dwidux + frac2*dwidux
       dwiduy = frac1*3.0_dp*wi2*dwiduy + frac2*dwiduy
       dwiduz = frac1*3.0_dp*wi2*dwiduz + frac2*dwiduz

       dwjdux = ( 4.0_dp*(yj*zj2 + 0.5_dp*yj*(xj2-yj2))*rI3 )
       dwjduy = ( 4.0_dp*(xj*zj2 - 0.5_dp*xj*(xj2-yj2))*rI3 )
       dwjduz = ( -8.0_dp*xj*yj*zj*rI3 )

       dwjdux = frac1*3.0_dp*wj2*dwjdux + frac2*dwjdux
       dwjduy = frac1*3.0_dp*wj2*dwjduy + frac2*dwjduy
       dwjduz = frac1*3.0_dp*wj2*dwjduz + frac2*dwjduz

       ! do the torques first since they are easy:
       ! remember that these are still in the body fixed axes

       txi = 0.5_dp*(v0*s*dwidux)*sw
       tyi = 0.5_dp*(v0*s*dwiduy)*sw
       tzi = 0.5_dp*(v0*s*dwiduz)*sw

       txj = 0.5_dp*(v0*s*dwjdux)*sw
       tyj = 0.5_dp*(v0*s*dwjduy)*sw
       tzj = 0.5_dp*(v0*s*dwjduz)*sw
 
       ! go back to lab frame using transpose of rotation matrix:

       t1(1) = t1(1) + a1(1)*txi + a1(4)*tyi + a1(7)*tzi
       t1(2) = t1(2) + a1(2)*txi + a1(5)*tyi + a1(8)*tzi
       t1(3) = t1(3) + a1(3)*txi + a1(6)*tyi + a1(9)*tzi

       t2(1) = t2(1) + a2(1)*txj + a2(4)*tyj + a2(7)*tzj
       t2(2) = t2(2) + a2(2)*txj + a2(5)*tyj + a2(8)*tzj
       t2(3) = t2(3) + a2(3)*txj + a2(6)*tyj + a2(9)*tzj

       ! Now, on to the forces:

       ! first rotate the i terms back into the lab frame:

       radcomxi = (v0*s*dwidx)*sw
       radcomyi = (v0*s*dwidy)*sw
       radcomzi = (v0*s*dwidz)*sw

       radcomxj = (v0*s*dwjdx)*sw
       radcomyj = (v0*s*dwjdy)*sw
       radcomzj = (v0*s*dwjdz)*sw

       fxii = a1(1)*(radcomxi) + a1(4)*(radcomyi) + a1(7)*(radcomzi)
       fyii = a1(2)*(radcomxi) + a1(5)*(radcomyi) + a1(8)*(radcomzi)
       fzii = a1(3)*(radcomxi) + a1(6)*(radcomyi) + a1(9)*(radcomzi)

       fxjj = a2(1)*(radcomxj) + a2(4)*(radcomyj) + a2(7)*(radcomzj)
       fyjj = a2(2)*(radcomxj) + a2(5)*(radcomyj) + a2(8)*(radcomzj)
       fzjj = a2(3)*(radcomxj) + a2(6)*(radcomyj) + a2(9)*(radcomzj)

       fxij = -fxii
       fyij = -fyii
       fzij = -fzii

       fxji = -fxjj
       fyji = -fyjj
       fzji = -fzjj

       ! now assemble these with the radial-only terms:

       fxradial = 0.5_dp*(v0*dsdr*w*drdx + fxii + fxji)
       fyradial = 0.5_dp*(v0*dsdr*w*drdy + fyii + fyji)
       fzradial = 0.5_dp*(v0*dsdr*w*drdz + fzii + fzji)

       f1(1) = f1(1) + fxradial
       f1(2) = f1(2) + fyradial
       f1(3) = f1(3) + fzradial

    endif
  end subroutine do_sticky_power_pair

end module sticky
