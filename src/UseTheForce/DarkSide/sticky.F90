!! This Module Calculates forces due to SSD potential and VDW interactions
!! [Chandra and Ichiye, J. Chem. Phys. 111, 2701 (1999)].

!! This module contains the Public procedures:


!! Corresponds to the force field defined in ssd_FF.cpp 
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @author Christopher Fennel
!! @author J. Daniel Gezelter
!! @version $Id: sticky.F90,v 1.2 2004-10-20 21:52:20 gezelter Exp $, $Date: 2004-10-20 21:52:20 $, $Name: not supported by cvs2svn $, $Revision: 1.2 $

module sticky_pair

  use force_globals
  use definitions
  use simulation
#ifdef IS_MPI
  use mpiSimulation
#endif

  implicit none

  PRIVATE

  logical, save :: sticky_initialized = .false.
  real( kind = dp ), save :: SSD_w0 = 0.0_dp
  real( kind = dp ), save :: SSD_v0 = 0.0_dp
  real( kind = dp ), save :: SSD_v0p = 0.0_dp
  real( kind = dp ), save :: SSD_rl = 0.0_dp
  real( kind = dp ), save :: SSD_ru = 0.0_dp
  real( kind = dp ), save :: SSD_rlp = 0.0_dp
  real( kind = dp ), save :: SSD_rup = 0.0_dp
  real( kind = dp ), save :: SSD_rbig = 0.0_dp

  public :: check_sticky_FF
  public :: set_sticky_params
  public :: do_sticky_pair

contains

  subroutine check_sticky_FF(status)
    integer :: status 
    status = -1
    if (sticky_initialized) status = 0
    return
  end subroutine check_sticky_FF

  subroutine set_sticky_params(sticky_w0, sticky_v0, sticky_v0p, &
       sticky_rl, sticky_ru, sticky_rlp, sticky_rup)

    real( kind = dp ), intent(in) :: sticky_w0, sticky_v0, sticky_v0p
    real( kind = dp ), intent(in) :: sticky_rl, sticky_ru
    real( kind = dp ), intent(in) :: sticky_rlp, sticky_rup
    
    ! we could pass all 5 parameters if we felt like it...
    
    SSD_w0 = sticky_w0
    SSD_v0 = sticky_v0
    SSD_v0p = sticky_v0p
    SSD_rl = sticky_rl
    SSD_ru = sticky_ru
    SSD_rlp = sticky_rlp
    SSD_rup = sticky_rup

    if (SSD_ru .gt. SSD_rup) then
       SSD_rbig = SSD_ru
    else
       SSD_rbig = SSD_rup
    endif
   
    sticky_initialized = .true.
    return
  end subroutine set_sticky_params

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

    if (.not.sticky_initialized) then
       write(*,*) 'Sticky forces not initialized!'
       return
    endif


    if ( rij .LE. SSD_rbig ) then

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

       call calc_sw_fnc(rij, s, sp, dsdr, dspdr)

       wi = 2.0d0*(xi2-yi2)*zi / r3
       wj = 2.0d0*(xj2-yj2)*zj / r3
       w = wi+wj

       zif = zi/rij - 0.6d0
       zis = zi/rij + 0.8d0

       zjf = zj/rij - 0.6d0
       zjs = zj/rij + 0.8d0

       wip = zif*zif*zis*zis - SSD_w0
       wjp = zjf*zjf*zjs*zjs - SSD_w0
       wp = wip + wjp

       vpair = vpair + 0.5d0*(SSD_v0*s*w + SSD_v0p*sp*wp)
       if (do_pot) then
#ifdef IS_MPI 
          pot_row(atom1) = pot_row(atom1) + 0.25d0*(SSD_v0*s*w + SSD_v0p*sp*wp)*sw
          pot_col(atom2) = pot_col(atom2) + 0.25d0*(SSD_v0*s*w + SSD_v0p*sp*wp)*sw
#else
          pot = pot + 0.5d0*(SSD_v0*s*w + SSD_v0p*sp*wp)*sw
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

       txi = 0.5d0*(SSD_v0*s*dwidux + SSD_v0p*sp*dwipdux)*sw
       tyi = 0.5d0*(SSD_v0*s*dwiduy + SSD_v0p*sp*dwipduy)*sw
       tzi = 0.5d0*(SSD_v0*s*dwiduz + SSD_v0p*sp*dwipduz)*sw

       txj = 0.5d0*(SSD_v0*s*dwjdux + SSD_v0p*sp*dwjpdux)*sw
       tyj = 0.5d0*(SSD_v0*s*dwjduy + SSD_v0p*sp*dwjpduy)*sw
       tzj = 0.5d0*(SSD_v0*s*dwjduz + SSD_v0p*sp*dwjpduz)*sw

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

       radcomxi = (SSD_v0*s*dwidx+SSD_v0p*sp*dwipdx)*sw
       radcomyi = (SSD_v0*s*dwidy+SSD_v0p*sp*dwipdy)*sw
       radcomzi = (SSD_v0*s*dwidz+SSD_v0p*sp*dwipdz)*sw

       radcomxj = (SSD_v0*s*dwjdx+SSD_v0p*sp*dwjpdx)*sw
       radcomyj = (SSD_v0*s*dwjdy+SSD_v0p*sp*dwjpdy)*sw
       radcomzj = (SSD_v0*s*dwjdz+SSD_v0p*sp*dwjpdz)*sw

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

       fxradial = 0.5d0*(SSD_v0*dsdr*drdx*w + SSD_v0p*dspdr*drdx*wp + fxii + fxji)
       fyradial = 0.5d0*(SSD_v0*dsdr*drdy*w + SSD_v0p*dspdr*drdy*wp + fyii + fyji)
       fzradial = 0.5d0*(SSD_v0*dsdr*drdz*w + SSD_v0p*dspdr*drdz*wp + fzii + fzji)

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
  subroutine calc_sw_fnc(r, s, sp, dsdr, dspdr)
    
    real (kind=dp), intent(in) :: r
    real (kind=dp), intent(inout) :: s, sp, dsdr, dspdr
    
    ! distances must be in angstroms
    
    if (r.lt.SSD_rl) then
       s = 1.0d0
       dsdr = 0.0d0
    elseif (r.gt.SSD_ru) then
       s = 0.0d0
       dsdr = 0.0d0
    else
       s = ((SSD_ru + 2.0d0*r - 3.0d0*SSD_rl) * (SSD_ru-r)**2) / &
            ((SSD_ru - SSD_rl)**3)
       dsdr = 6.0d0*(r-SSD_ru)*(r-SSD_rl)/((SSD_ru - SSD_rl)**3)
    endif

    if (r.lt.SSD_rlp) then
       sp = 1.0d0       
       dspdr = 0.0d0
    elseif (r.gt.SSD_rup) then
       sp = 0.0d0
       dspdr = 0.0d0
    else
       sp = ((SSD_rup + 2.0d0*r - 3.0d0*SSD_rlp) * (SSD_rup-r)**2) / &
            ((SSD_rup - SSD_rlp)**3)
       dspdr = 6.0d0*(r-SSD_rup)*(r-SSD_rlp)/((SSD_rup - SSD_rlp)**3)       
    endif
    
    return
  end subroutine calc_sw_fnc
end module sticky_pair

  subroutine makeStickyType(sticky_w0, sticky_v0, sticky_v0p, &
       sticky_rl, sticky_ru, sticky_rlp, sticky_rup)
    use definitions, ONLY : dp   
    use sticky_pair, ONLY : set_sticky_params
    real( kind = dp ), intent(inout) :: sticky_w0, sticky_v0, sticky_v0p
    real( kind = dp ), intent(inout) :: sticky_rl, sticky_ru
    real( kind = dp ), intent(inout) :: sticky_rlp, sticky_rup
    
    call set_sticky_params(sticky_w0, sticky_v0, sticky_v0p, &
       sticky_rl, sticky_ru, sticky_rlp, sticky_rup)
       
  end subroutine makeStickyType
