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

!! This Module Calculates forces due to SSD potential and VDW interactions
!! [Chandra and Ichiye, J. Chem. Phys. 111, 2701 (1999)].

!! This module contains the Public procedures:


!! Corresponds to the force field defined in ssd_FF.cpp 
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @author Christopher Fennell
!! @author J. Daniel Gezelter
!! @version $Id: sticky.F90,v 1.8 2005-05-05 14:47:35 chrisfen Exp $, $Date: 2005-05-05 14:47:35 $, $Name: not supported by cvs2svn $, $Revision: 1.8 $

module sticky

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

  public :: newStickyType
  public :: do_sticky_pair
  public :: destroyStickyTypes
  public :: do_sticky_power_pair


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
  end type StickyList

  type(StickyList), dimension(:),allocatable :: StickyMap

contains

  subroutine newStickyType(c_ident, w0, v0, v0p, rl, ru, rlp, rup, isError)

    integer, intent(in) :: c_ident
    integer, intent(inout) :: isError
    real( kind = dp ), intent(in) :: w0, v0, v0p
    real( kind = dp ), intent(in) :: rl, ru
    real( kind = dp ), intent(in) :: rlp, rup
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

    return
  end subroutine newStickyType

  subroutine do_sticky_pair(atom1, atom2, d, rij, r2, sw, vpair, fpair, &
       pot, A, f, t, do_pot)

    !! This routine does only the sticky portion of the SSD potential
    !! [Chandra and Ichiye, J. Chem. Phys. 111, 2701 (1999)].
    !! The Lennard-Jones and dipolar interaction must be handled separately.

    !! We assume that the rotation matrices have already been calculated
    !! and placed in the A array.

    !! i and j are pointers to the two SSD atoms

    integer, intent(in) :: atom1, atom2
    real (kind=dp), intent(inout) :: rij, r2
    real (kind=dp), dimension(3), intent(in) :: d
    real (kind=dp), dimension(3), intent(inout) :: fpair
    real (kind=dp) :: pot, vpair, sw
    real (kind=dp), dimension(9,nLocal) :: A
    real (kind=dp), dimension(3,nLocal) :: f
    real (kind=dp), dimension(3,nLocal) :: t
    logical, intent(in) :: do_pot

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
    integer :: me1, me2
    real (kind=dp) :: w0, v0, v0p, rl, ru, rlp, rup, rbig

    if (.not.allocated(StickyMap)) then
       call handleError("sticky", "no StickyMap was present before first call of do_sticky_pair!")
       return
    end if

#ifdef IS_MPI
    me1 = atid_Row(atom1)
    me2 = atid_Col(atom2)
#else
    me1 = atid(atom1)
    me2 = atid(atom2)
#endif

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

#ifdef IS_MPI
       ! rotate the inter-particle separation into the two different 
       ! body-fixed coordinate systems:

       xi = A_row(1,atom1)*d(1) + A_row(2,atom1)*d(2) + A_row(3,atom1)*d(3)
       yi = A_row(4,atom1)*d(1) + A_row(5,atom1)*d(2) + A_row(6,atom1)*d(3)
       zi = A_row(7,atom1)*d(1) + A_row(8,atom1)*d(2) + A_row(9,atom1)*d(3)

       ! negative sign because this is the vector from j to i:

       xj = -(A_Col(1,atom2)*d(1) + A_Col(2,atom2)*d(2) + A_Col(3,atom2)*d(3))
       yj = -(A_Col(4,atom2)*d(1) + A_Col(5,atom2)*d(2) + A_Col(6,atom2)*d(3))
       zj = -(A_Col(7,atom2)*d(1) + A_Col(8,atom2)*d(2) + A_Col(9,atom2)*d(3))
#else
       ! rotate the inter-particle separation into the two different 
       ! body-fixed coordinate systems:

       xi = a(1,atom1)*d(1) + a(2,atom1)*d(2) + a(3,atom1)*d(3)
       yi = a(4,atom1)*d(1) + a(5,atom1)*d(2) + a(6,atom1)*d(3)
       zi = a(7,atom1)*d(1) + a(8,atom1)*d(2) + a(9,atom1)*d(3)

       ! negative sign because this is the vector from j to i:

       xj = -(a(1,atom2)*d(1) + a(2,atom2)*d(2) + a(3,atom2)*d(3))
       yj = -(a(4,atom2)*d(1) + a(5,atom2)*d(2) + a(6,atom2)*d(3))
       zj = -(a(7,atom2)*d(1) + a(8,atom2)*d(2) + a(9,atom2)*d(3))
#endif

       xi2 = xi*xi
       yi2 = yi*yi
       zi2 = zi*zi

       xj2 = xj*xj
       yj2 = yj*yj
       zj2 = zj*zj

       call calc_sw_fnc(rij, rl, ru, rlp, rup, s, sp, dsdr, dspdr)

       wi = 2.0d0*(xi2-yi2)*zi / r3
       wj = 2.0d0*(xj2-yj2)*zj / r3
       w = wi+wj

       zif = zi/rij - 0.6d0
       zis = zi/rij + 0.8d0

       zjf = zj/rij - 0.6d0
       zjs = zj/rij + 0.8d0

       wip = zif*zif*zis*zis - w0
       wjp = zjf*zjf*zjs*zjs - w0
       wp = wip + wjp

       vpair = vpair + 0.5d0*(v0*s*w + v0p*sp*wp)
       if (do_pot) then
#ifdef IS_MPI 
          pot_row(atom1) = pot_row(atom1) + 0.25d0*(v0*s*w + v0p*sp*wp)*sw
          pot_col(atom2) = pot_col(atom2) + 0.25d0*(v0*s*w + v0p*sp*wp)*sw
#else
          pot = pot + 0.5d0*(v0*s*w + v0p*sp*wp)*sw
#endif  
       endif

       dwidx =   4.0d0*xi*zi/r3  - 6.0d0*xi*zi*(xi2-yi2)/r5
       dwidy = - 4.0d0*yi*zi/r3  - 6.0d0*yi*zi*(xi2-yi2)/r5
       dwidz =   2.0d0*(xi2-yi2)/r3  - 6.0d0*zi2*(xi2-yi2)/r5

       dwjdx =   4.0d0*xj*zj/r3  - 6.0d0*xj*zj*(xj2-yj2)/r5
       dwjdy = - 4.0d0*yj*zj/r3  - 6.0d0*yj*zj*(xj2-yj2)/r5
       dwjdz =   2.0d0*(xj2-yj2)/r3  - 6.0d0*zj2*(xj2-yj2)/r5

       uglyi = zif*zif*zis + zif*zis*zis
       uglyj = zjf*zjf*zjs + zjf*zjs*zjs

       dwipdx = -2.0d0*xi*zi*uglyi/r3
       dwipdy = -2.0d0*yi*zi*uglyi/r3
       dwipdz = 2.0d0*(1.0d0/rij - zi2/r3)*uglyi

       dwjpdx = -2.0d0*xj*zj*uglyj/r3
       dwjpdy = -2.0d0*yj*zj*uglyj/r3
       dwjpdz = 2.0d0*(1.0d0/rij - zj2/r3)*uglyj

       dwidux = 4.0d0*(yi*zi2 + 0.5d0*yi*(xi2-yi2))/r3
       dwiduy = 4.0d0*(xi*zi2 - 0.5d0*xi*(xi2-yi2))/r3
       dwiduz = - 8.0d0*xi*yi*zi/r3

       dwjdux = 4.0d0*(yj*zj2 + 0.5d0*yj*(xj2-yj2))/r3
       dwjduy = 4.0d0*(xj*zj2 - 0.5d0*xj*(xj2-yj2))/r3
       dwjduz = - 8.0d0*xj*yj*zj/r3

       dwipdux =  2.0d0*yi*uglyi/rij
       dwipduy = -2.0d0*xi*uglyi/rij
       dwipduz =  0.0d0

       dwjpdux =  2.0d0*yj*uglyj/rij
       dwjpduy = -2.0d0*xj*uglyj/rij
       dwjpduz =  0.0d0

       ! do the torques first since they are easy:
       ! remember that these are still in the body fixed axes

       txi = 0.5d0*(v0*s*dwidux + v0p*sp*dwipdux)*sw
       tyi = 0.5d0*(v0*s*dwiduy + v0p*sp*dwipduy)*sw
       tzi = 0.5d0*(v0*s*dwiduz + v0p*sp*dwipduz)*sw

       txj = 0.5d0*(v0*s*dwjdux + v0p*sp*dwjpdux)*sw
       tyj = 0.5d0*(v0*s*dwjduy + v0p*sp*dwjpduy)*sw
       tzj = 0.5d0*(v0*s*dwjduz + v0p*sp*dwjpduz)*sw

       ! go back to lab frame using transpose of rotation matrix:

#ifdef IS_MPI
       t_Row(1,atom1) = t_Row(1,atom1) + a_Row(1,atom1)*txi + &
            a_Row(4,atom1)*tyi + a_Row(7,atom1)*tzi
       t_Row(2,atom1) = t_Row(2,atom1) + a_Row(2,atom1)*txi + &
            a_Row(5,atom1)*tyi + a_Row(8,atom1)*tzi
       t_Row(3,atom1) = t_Row(3,atom1) + a_Row(3,atom1)*txi + &
            a_Row(6,atom1)*tyi + a_Row(9,atom1)*tzi

       t_Col(1,atom2) = t_Col(1,atom2) + a_Col(1,atom2)*txj + &
            a_Col(4,atom2)*tyj + a_Col(7,atom2)*tzj
       t_Col(2,atom2) = t_Col(2,atom2) + a_Col(2,atom2)*txj + &
            a_Col(5,atom2)*tyj + a_Col(8,atom2)*tzj
       t_Col(3,atom2) = t_Col(3,atom2) + a_Col(3,atom2)*txj + &
            a_Col(6,atom2)*tyj + a_Col(9,atom2)*tzj
#else
       t(1,atom1) = t(1,atom1) + a(1,atom1)*txi + a(4,atom1)*tyi + a(7,atom1)*tzi
       t(2,atom1) = t(2,atom1) + a(2,atom1)*txi + a(5,atom1)*tyi + a(8,atom1)*tzi
       t(3,atom1) = t(3,atom1) + a(3,atom1)*txi + a(6,atom1)*tyi + a(9,atom1)*tzi

       t(1,atom2) = t(1,atom2) + a(1,atom2)*txj + a(4,atom2)*tyj + a(7,atom2)*tzj
       t(2,atom2) = t(2,atom2) + a(2,atom2)*txj + a(5,atom2)*tyj + a(8,atom2)*tzj
       t(3,atom2) = t(3,atom2) + a(3,atom2)*txj + a(6,atom2)*tyj + a(9,atom2)*tzj
#endif    
       ! Now, on to the forces:

       ! first rotate the i terms back into the lab frame:

       radcomxi = (v0*s*dwidx+v0p*sp*dwipdx)*sw
       radcomyi = (v0*s*dwidy+v0p*sp*dwipdy)*sw
       radcomzi = (v0*s*dwidz+v0p*sp*dwipdz)*sw

       radcomxj = (v0*s*dwjdx+v0p*sp*dwjpdx)*sw
       radcomyj = (v0*s*dwjdy+v0p*sp*dwjpdy)*sw
       radcomzj = (v0*s*dwjdz+v0p*sp*dwjpdz)*sw

#ifdef IS_MPI    
       fxii = a_Row(1,atom1)*(radcomxi) + &
            a_Row(4,atom1)*(radcomyi) + &
            a_Row(7,atom1)*(radcomzi)
       fyii = a_Row(2,atom1)*(radcomxi) + &
            a_Row(5,atom1)*(radcomyi) + &
            a_Row(8,atom1)*(radcomzi)
       fzii = a_Row(3,atom1)*(radcomxi) + &
            a_Row(6,atom1)*(radcomyi) + &
            a_Row(9,atom1)*(radcomzi)

       fxjj = a_Col(1,atom2)*(radcomxj) + &
            a_Col(4,atom2)*(radcomyj) + &
            a_Col(7,atom2)*(radcomzj)
       fyjj = a_Col(2,atom2)*(radcomxj) + &
            a_Col(5,atom2)*(radcomyj) + &
            a_Col(8,atom2)*(radcomzj)
       fzjj = a_Col(3,atom2)*(radcomxj)+ &
            a_Col(6,atom2)*(radcomyj) + &
            a_Col(9,atom2)*(radcomzj)
#else
       fxii = a(1,atom1)*(radcomxi) + &
            a(4,atom1)*(radcomyi) + &
            a(7,atom1)*(radcomzi)
       fyii = a(2,atom1)*(radcomxi) + &
            a(5,atom1)*(radcomyi) + &
            a(8,atom1)*(radcomzi)
       fzii = a(3,atom1)*(radcomxi) + &
            a(6,atom1)*(radcomyi) + &
            a(9,atom1)*(radcomzi)

       fxjj = a(1,atom2)*(radcomxj) + &
            a(4,atom2)*(radcomyj) + &
            a(7,atom2)*(radcomzj)
       fyjj = a(2,atom2)*(radcomxj) + &
            a(5,atom2)*(radcomyj) + &
            a(8,atom2)*(radcomzj)
       fzjj = a(3,atom2)*(radcomxj)+ &
            a(6,atom2)*(radcomyj) + &
            a(9,atom2)*(radcomzj)
#endif

       fxij = -fxii
       fyij = -fyii
       fzij = -fzii

       fxji = -fxjj
       fyji = -fyjj
       fzji = -fzjj

       ! now assemble these with the radial-only terms:

       fxradial = 0.5d0*(v0*dsdr*drdx*w + v0p*dspdr*drdx*wp + fxii + fxji)
       fyradial = 0.5d0*(v0*dsdr*drdy*w + v0p*dspdr*drdy*wp + fyii + fyji)
       fzradial = 0.5d0*(v0*dsdr*drdz*w + v0p*dspdr*drdz*wp + fzii + fzji)

#ifdef IS_MPI
       f_Row(1,atom1) = f_Row(1,atom1) + fxradial
       f_Row(2,atom1) = f_Row(2,atom1) + fyradial
       f_Row(3,atom1) = f_Row(3,atom1) + fzradial

       f_Col(1,atom2) = f_Col(1,atom2) - fxradial
       f_Col(2,atom2) = f_Col(2,atom2) - fyradial
       f_Col(3,atom2) = f_Col(3,atom2) - fzradial
#else
       f(1,atom1) = f(1,atom1) + fxradial
       f(2,atom1) = f(2,atom1) + fyradial
       f(3,atom1) = f(3,atom1) + fzradial

       f(1,atom2) = f(1,atom2) - fxradial
       f(2,atom2) = f(2,atom2) - fyradial
       f(3,atom2) = f(3,atom2) - fzradial
#endif

#ifdef IS_MPI
       id1 = AtomRowToGlobal(atom1)
       id2 = AtomColToGlobal(atom2)
#else
       id1 = atom1
       id2 = atom2
#endif

       if (molMembershipList(id1) .ne. molMembershipList(id2)) then

          fpair(1) = fpair(1) + fxradial
          fpair(2) = fpair(2) + fyradial
          fpair(3) = fpair(3) + fzradial

       endif
    endif
  end subroutine do_sticky_pair

  !! calculates the switching functions and their derivatives for a given
  subroutine calc_sw_fnc(r, rl, ru, rlp, rup, s, sp, dsdr, dspdr)

    real (kind=dp), intent(in) :: r, rl, ru, rlp, rup
    real (kind=dp), intent(inout) :: s, sp, dsdr, dspdr

    ! distances must be in angstroms

    if (r.lt.rl) then
       s = 1.0d0
       dsdr = 0.0d0
    elseif (r.gt.ru) then
       s = 0.0d0
       dsdr = 0.0d0
    else
       s = ((ru + 2.0d0*r - 3.0d0*rl) * (ru-r)**2) / &
            ((ru - rl)**3)
       dsdr = 6.0d0*(r-ru)*(r-rl)/((ru - rl)**3)
    endif

    if (r.lt.rlp) then
       sp = 1.0d0       
       dspdr = 0.0d0
    elseif (r.gt.rup) then
       sp = 0.0d0
       dspdr = 0.0d0
    else
       sp = ((rup + 2.0d0*r - 3.0d0*rlp) * (rup-r)**2) / &
            ((rup - rlp)**3)
       dspdr = 6.0d0*(r-rup)*(r-rlp)/((rup - rlp)**3)       
    endif

    return
  end subroutine calc_sw_fnc

  subroutine destroyStickyTypes()   
    if(allocated(StickyMap)) deallocate(StickyMap)
  end subroutine destroyStickyTypes
  
    subroutine do_sticky_power_pair(atom1, atom2, d, rij, r2, sw, vpair, fpair, &
       pot, A, f, t, do_pot)
    !! We assume that the rotation matrices have already been calculated
    !! and placed in the A array.

    !! i and j are pointers to the two SSD atoms

    integer, intent(in) :: atom1, atom2
    real (kind=dp), intent(inout) :: rij, r2
    real (kind=dp), dimension(3), intent(in) :: d
    real (kind=dp), dimension(3), intent(inout) :: fpair
    real (kind=dp) :: pot, vpair, sw
    real (kind=dp), dimension(9,nLocal) :: A
    real (kind=dp), dimension(3,nLocal) :: f
    real (kind=dp), dimension(3,nLocal) :: t
    logical, intent(in) :: do_pot

    real (kind=dp) :: xi, yi, zi, xj, yj, zj, xi2, yi2, zi2, xj2, yj2, zj2
    real (kind=dp) :: r3, r5, r6, s, sp, dsdr, dspdr
    real (kind=dp) :: wi, wj, w, wip, wjp, wp, wi2, wj2
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
    integer :: me1, me2
    real (kind=dp) :: w0, v0, v0p, rl, ru, rlp, rup, rbig

    if (.not.allocated(StickyMap)) then
       call handleError("sticky", "no StickyMap was present before first call of do_sticky_power_pair!")
       return
    end if

#ifdef IS_MPI
    me1 = atid_Row(atom1)
    me2 = atid_Col(atom2)
#else
    me1 = atid(atom1)
    me2 = atid(atom2)
#endif
 
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

#ifdef IS_MPI
       ! rotate the inter-particle separation into the two different 
       ! body-fixed coordinate systems:

       xi = A_row(1,atom1)*d(1) + A_row(2,atom1)*d(2) + A_row(3,atom1)*d(3)
       yi = A_row(4,atom1)*d(1) + A_row(5,atom1)*d(2) + A_row(6,atom1)*d(3)
       zi = A_row(7,atom1)*d(1) + A_row(8,atom1)*d(2) + A_row(9,atom1)*d(3)

       ! negative sign because this is the vector from j to i:

       xj = -(A_Col(1,atom2)*d(1) + A_Col(2,atom2)*d(2) + A_Col(3,atom2)*d(3))
       yj = -(A_Col(4,atom2)*d(1) + A_Col(5,atom2)*d(2) + A_Col(6,atom2)*d(3))
       zj = -(A_Col(7,atom2)*d(1) + A_Col(8,atom2)*d(2) + A_Col(9,atom2)*d(3))
#else
       ! rotate the inter-particle separation into the two different 
       ! body-fixed coordinate systems:

       xi = a(1,atom1)*d(1) + a(2,atom1)*d(2) + a(3,atom1)*d(3)
       yi = a(4,atom1)*d(1) + a(5,atom1)*d(2) + a(6,atom1)*d(3)
       zi = a(7,atom1)*d(1) + a(8,atom1)*d(2) + a(9,atom1)*d(3)

       ! negative sign because this is the vector from j to i:

       xj = -(a(1,atom2)*d(1) + a(2,atom2)*d(2) + a(3,atom2)*d(3))
       yj = -(a(4,atom2)*d(1) + a(5,atom2)*d(2) + a(6,atom2)*d(3))
       zj = -(a(7,atom2)*d(1) + a(8,atom2)*d(2) + a(9,atom2)*d(3))
#endif

       xi2 = xi*xi
       yi2 = yi*yi
       zi2 = zi*zi

       xj2 = xj*xj
       yj2 = yj*yj
       zj2 = zj*zj

       call calc_sw_fnc(rij, rl, ru, rlp, rup, s, sp, dsdr, dspdr)

       wi = 2.0d0*(xi2-yi2)*zi / r3
       wj = 2.0d0*(xj2-yj2)*zj / r3
       !rootwi = sqrt(abs(wi))
       !rootwj = sqrt(abs(wj))
       wi2 = wi*wi
       wj2 = wj*wj

       
       w = wi*wi2+wj*wj2

       zif = zi/rij - 0.6d0
       zis = zi/rij + 0.8d0

       zjf = zj/rij - 0.6d0
       zjs = zj/rij + 0.8d0

       wip = zif*zif*zis*zis - w0
       wjp = zjf*zjf*zjs*zjs - w0
       wp = wip + wjp

       vpair = vpair + 0.5d0*(v0*s*w + v0p*sp*wp)
       if (do_pot) then
#ifdef IS_MPI 
          pot_row(atom1) = pot_row(atom1) + 0.25d0*(v0*s*w + v0p*sp*wp)*sw
          pot_col(atom2) = pot_col(atom2) + 0.25d0*(v0*s*w + v0p*sp*wp)*sw
#else
          pot = pot + 0.5d0*(v0*s*w + v0p*sp*wp)*sw
#endif  
       endif

!       dwidx = 1.5d0*rootwi*( 4.0d0*xi*zi/r3 - 6.0d0*xi*zi*(xi2-yi2)/r5 )
!       dwidy = 1.5d0*rootwi*( -4.0d0*yi*zi/r3 - 6.0d0*yi*zi*(xi2-yi2)/r5 )
!       dwidz = 1.5d0*rootwi*( 2.0d0*(xi2-yi2)/r3 - 6.0d0*zi2*(xi2-yi2)/r5 )

!       dwjdx = 1.5d0*rootwj*( 4.0d0*xj*zj/r3  - 6.0d0*xj*zj*(xj2-yj2)/r5 )
!       dwjdy = 1.5d0*rootwj*( -4.0d0*yj*zj/r3  - 6.0d0*yj*zj*(xj2-yj2)/r5 )
!       dwjdz = 1.5d0*rootwj*( 2.0d0*(xj2-yj2)/r3  - 6.0d0*zj2*(xj2-yj2)/r5 )
       
       dwidx = 3.0d0*wi2*( 4.0d0*xi*zi/r3 - 6.0d0*xi*zi*(xi2-yi2)/r5 )
       dwidy = 3.0d0*wi2*( -4.0d0*yi*zi/r3 - 6.0d0*yi*zi*(xi2-yi2)/r5 )
       dwidz = 3.0d0*wi2*( 2.0d0*(xi2-yi2)/r3 - 6.0d0*zi2*(xi2-yi2)/r5 )

       dwjdx = 3.0d0*wj2*( 4.0d0*xj*zj/r3  - 6.0d0*xj*zj*(xj2-yj2)/r5 )
       dwjdy = 3.0d0*wj2*( -4.0d0*yj*zj/r3  - 6.0d0*yj*zj*(xj2-yj2)/r5 )
       dwjdz = 3.0d0*wj2*( 2.0d0*(xj2-yj2)/r3  - 6.0d0*zj2*(xj2-yj2)/r5 )

       uglyi = zif*zif*zis + zif*zis*zis
       uglyj = zjf*zjf*zjs + zjf*zjs*zjs

       dwipdx = -2.0d0*xi*zi*uglyi/r3
       dwipdy = -2.0d0*yi*zi*uglyi/r3
       dwipdz = 2.0d0*(1.0d0/rij - zi2/r3)*uglyi

       dwjpdx = -2.0d0*xj*zj*uglyj/r3
       dwjpdy = -2.0d0*yj*zj*uglyj/r3
       dwjpdz = 2.0d0*(1.0d0/rij - zj2/r3)*uglyj

!       dwidux = 1.5d0*rootwi*( 4.0d0*(yi*zi2 + 0.5d0*yi*(xi2-yi2))/r3 )
!       dwiduy = 1.5d0*rootwi*( 4.0d0*(xi*zi2 - 0.5d0*xi*(xi2-yi2))/r3 )
!       dwiduz = 1.5d0*rootwi*( -8.0d0*xi*yi*zi/r3 )

!       dwjdux = 1.5d0*rootwj*( 4.0d0*(yj*zj2 + 0.5d0*yj*(xj2-yj2))/r3 )
!       dwjduy = 1.5d0*rootwj*( 4.0d0*(xj*zj2 - 0.5d0*xj*(xj2-yj2))/r3 )
!       dwjduz = 1.5d0*rootwj*( -8.0d0*xj*yj*zj/r3 )
       
       dwidux = 3.0d0*wi2*( 4.0d0*(yi*zi2 + 0.5d0*yi*(xi2-yi2))/r3 )
       dwiduy = 3.0d0*wi2*( 4.0d0*(xi*zi2 - 0.5d0*xi*(xi2-yi2))/r3 )
       dwiduz = 3.0d0*wi2*( -8.0d0*xi*yi*zi/r3 )

       dwjdux = 3.0d0*wj2*( 4.0d0*(yj*zj2 + 0.5d0*yj*(xj2-yj2))/r3 )
       dwjduy = 3.0d0*wj2*( 4.0d0*(xj*zj2 - 0.5d0*xj*(xj2-yj2))/r3 )
       dwjduz = 3.0d0*wj2*( -8.0d0*xj*yj*zj/r3 )

       dwipdux =  2.0d0*yi*uglyi/rij
       dwipduy = -2.0d0*xi*uglyi/rij
       dwipduz =  0.0d0

       dwjpdux =  2.0d0*yj*uglyj/rij
       dwjpduy = -2.0d0*xj*uglyj/rij
       dwjpduz =  0.0d0

       ! do the torques first since they are easy:
       ! remember that these are still in the body fixed axes

       txi = 0.5d0*(v0*s*dwidux + v0p*sp*dwipdux)*sw
       tyi = 0.5d0*(v0*s*dwiduy + v0p*sp*dwipduy)*sw
       tzi = 0.5d0*(v0*s*dwiduz + v0p*sp*dwipduz)*sw

       txj = 0.5d0*(v0*s*dwjdux + v0p*sp*dwjpdux)*sw
       tyj = 0.5d0*(v0*s*dwjduy + v0p*sp*dwjpduy)*sw
       tzj = 0.5d0*(v0*s*dwjduz + v0p*sp*dwjpduz)*sw
 
       ! go back to lab frame using transpose of rotation matrix:

#ifdef IS_MPI
       t_Row(1,atom1) = t_Row(1,atom1) + a_Row(1,atom1)*txi + &
            a_Row(4,atom1)*tyi + a_Row(7,atom1)*tzi
       t_Row(2,atom1) = t_Row(2,atom1) + a_Row(2,atom1)*txi + &
            a_Row(5,atom1)*tyi + a_Row(8,atom1)*tzi
       t_Row(3,atom1) = t_Row(3,atom1) + a_Row(3,atom1)*txi + &
            a_Row(6,atom1)*tyi + a_Row(9,atom1)*tzi

       t_Col(1,atom2) = t_Col(1,atom2) + a_Col(1,atom2)*txj + &
            a_Col(4,atom2)*tyj + a_Col(7,atom2)*tzj
       t_Col(2,atom2) = t_Col(2,atom2) + a_Col(2,atom2)*txj + &
            a_Col(5,atom2)*tyj + a_Col(8,atom2)*tzj
       t_Col(3,atom2) = t_Col(3,atom2) + a_Col(3,atom2)*txj + &
            a_Col(6,atom2)*tyj + a_Col(9,atom2)*tzj
#else
       t(1,atom1) = t(1,atom1) + a(1,atom1)*txi + a(4,atom1)*tyi + a(7,atom1)*tzi
       t(2,atom1) = t(2,atom1) + a(2,atom1)*txi + a(5,atom1)*tyi + a(8,atom1)*tzi
       t(3,atom1) = t(3,atom1) + a(3,atom1)*txi + a(6,atom1)*tyi + a(9,atom1)*tzi

       t(1,atom2) = t(1,atom2) + a(1,atom2)*txj + a(4,atom2)*tyj + a(7,atom2)*tzj
       t(2,atom2) = t(2,atom2) + a(2,atom2)*txj + a(5,atom2)*tyj + a(8,atom2)*tzj
       t(3,atom2) = t(3,atom2) + a(3,atom2)*txj + a(6,atom2)*tyj + a(9,atom2)*tzj
#endif    
       ! Now, on to the forces:

       ! first rotate the i terms back into the lab frame:

       radcomxi = (v0*s*dwidx+v0p*sp*dwipdx)*sw
       radcomyi = (v0*s*dwidy+v0p*sp*dwipdy)*sw
       radcomzi = (v0*s*dwidz+v0p*sp*dwipdz)*sw

       radcomxj = (v0*s*dwjdx+v0p*sp*dwjpdx)*sw
       radcomyj = (v0*s*dwjdy+v0p*sp*dwjpdy)*sw
       radcomzj = (v0*s*dwjdz+v0p*sp*dwjpdz)*sw

#ifdef IS_MPI    
       fxii = a_Row(1,atom1)*(radcomxi) + &
            a_Row(4,atom1)*(radcomyi) + &
            a_Row(7,atom1)*(radcomzi)
       fyii = a_Row(2,atom1)*(radcomxi) + &
            a_Row(5,atom1)*(radcomyi) + &
            a_Row(8,atom1)*(radcomzi)
       fzii = a_Row(3,atom1)*(radcomxi) + &
            a_Row(6,atom1)*(radcomyi) + &
            a_Row(9,atom1)*(radcomzi)

       fxjj = a_Col(1,atom2)*(radcomxj) + &
            a_Col(4,atom2)*(radcomyj) + &
            a_Col(7,atom2)*(radcomzj)
       fyjj = a_Col(2,atom2)*(radcomxj) + &
            a_Col(5,atom2)*(radcomyj) + &
            a_Col(8,atom2)*(radcomzj)
       fzjj = a_Col(3,atom2)*(radcomxj)+ &
            a_Col(6,atom2)*(radcomyj) + &
            a_Col(9,atom2)*(radcomzj)
#else
       fxii = a(1,atom1)*(radcomxi) + &
            a(4,atom1)*(radcomyi) + &
            a(7,atom1)*(radcomzi)
       fyii = a(2,atom1)*(radcomxi) + &
            a(5,atom1)*(radcomyi) + &
            a(8,atom1)*(radcomzi)
       fzii = a(3,atom1)*(radcomxi) + &
            a(6,atom1)*(radcomyi) + &
            a(9,atom1)*(radcomzi)

       fxjj = a(1,atom2)*(radcomxj) + &
            a(4,atom2)*(radcomyj) + &
            a(7,atom2)*(radcomzj)
       fyjj = a(2,atom2)*(radcomxj) + &
            a(5,atom2)*(radcomyj) + &
            a(8,atom2)*(radcomzj)
       fzjj = a(3,atom2)*(radcomxj)+ &
            a(6,atom2)*(radcomyj) + &
            a(9,atom2)*(radcomzj)
#endif

       fxij = -fxii
       fyij = -fyii
       fzij = -fzii

       fxji = -fxjj
       fyji = -fyjj
       fzji = -fzjj

       ! now assemble these with the radial-only terms:

       fxradial = 0.5d0*(v0*dsdr*drdx*w + v0p*dspdr*drdx*wp + fxii + fxji)
       fyradial = 0.5d0*(v0*dsdr*drdy*w + v0p*dspdr*drdy*wp + fyii + fyji)
       fzradial = 0.5d0*(v0*dsdr*drdz*w + v0p*dspdr*drdz*wp + fzii + fzji)

#ifdef IS_MPI
       f_Row(1,atom1) = f_Row(1,atom1) + fxradial
       f_Row(2,atom1) = f_Row(2,atom1) + fyradial
       f_Row(3,atom1) = f_Row(3,atom1) + fzradial

       f_Col(1,atom2) = f_Col(1,atom2) - fxradial
       f_Col(2,atom2) = f_Col(2,atom2) - fyradial
       f_Col(3,atom2) = f_Col(3,atom2) - fzradial
#else
       f(1,atom1) = f(1,atom1) + fxradial
       f(2,atom1) = f(2,atom1) + fyradial
       f(3,atom1) = f(3,atom1) + fzradial

       f(1,atom2) = f(1,atom2) - fxradial
       f(2,atom2) = f(2,atom2) - fyradial
       f(3,atom2) = f(3,atom2) - fzradial
#endif

#ifdef IS_MPI
       id1 = AtomRowToGlobal(atom1)
       id2 = AtomColToGlobal(atom2)
#else
       id1 = atom1
       id2 = atom2
#endif

       if (molMembershipList(id1) .ne. molMembershipList(id2)) then

          fpair(1) = fpair(1) + fxradial
          fpair(2) = fpair(2) + fyradial
          fpair(3) = fpair(3) + fzradial

       endif
    endif
  end subroutine do_sticky_power_pair

end module sticky
