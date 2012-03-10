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


module gayberne
  use force_globals
  use definitions
  use simulation
  use atype_module
  use vector_class
  use linearalgebra
  use status
  use lj
  use fForceOptions
  
  implicit none

  private

#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"

  logical, save :: haveGBMap = .false.
  logical, save :: haveMixingMap = .false.
  real(kind=dp), save :: mu = 2.0_dp
  real(kind=dp), save :: nu = 1.0_dp


  public :: newGBtype
  public :: complete_GB_FF
  public :: do_gb_pair
  public :: getGayBerneCut
  public :: destroyGBtypes

  type :: GBtype
     integer          :: atid
     real(kind = dp ) :: d
     real(kind = dp ) :: l
     real(kind = dp ) :: eps
     real(kind = dp ) :: eps_ratio 
     real(kind = dp ) :: dw
     logical          :: isLJ
  end type GBtype
  
  type, private :: GBList
     integer               :: nGBtypes = 0
     integer               :: currentGBtype = 0
     type(GBtype), pointer :: GBtypes(:)      => null()
     integer, pointer      :: atidToGBtype(:) => null()
  end type GBList
  
  type(GBList), save :: GBMap
  
  type :: GBMixParameters
     real(kind=DP) :: sigma0
     real(kind=DP) :: eps0
     real(kind=DP) :: dw
     real(kind=DP) :: x2
     real(kind=DP) :: xa2
     real(kind=DP) :: xai2
     real(kind=DP) :: xp2
     real(kind=DP) :: xpap2
     real(kind=DP) :: xpapi2
  end type GBMixParameters
  
  type(GBMixParameters), dimension(:,:), allocatable :: GBMixingMap
  
contains
  
  subroutine newGBtype(c_ident, d, l, eps, eps_ratio, dw, status)
    
    integer, intent(in) :: c_ident
    real( kind = dp ), intent(in) :: d, l, eps, eps_ratio, dw
    integer, intent(out) :: status
    
    integer :: nGBTypes, nLJTypes, ntypes, myATID
    integer, pointer :: MatchList(:) => null()
    integer :: current, i
    status = 0
    
    if (.not.associated(GBMap%GBtypes)) then
       
       call getMatchingElementList(atypes, "is_GayBerne", .true., &
            nGBtypes, MatchList)
       
       call getMatchingElementList(atypes, "is_LennardJones", .true., &
            nLJTypes, MatchList)
       
       GBMap%nGBtypes = nGBtypes + nLJTypes
       
       allocate(GBMap%GBtypes(nGBtypes + nLJTypes))
       
       ntypes = getSize(atypes)
       
       allocate(GBMap%atidToGBtype(ntypes))       
    endif
    
    GBMap%currentGBtype = GBMap%currentGBtype + 1
    current = GBMap%currentGBtype
    
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)
    
    GBMap%atidToGBtype(myATID)       = current
    GBMap%GBtypes(current)%atid      = myATID
    GBMap%GBtypes(current)%d         = d
    GBMap%GBtypes(current)%l         = l
    GBMap%GBtypes(current)%eps       = eps
    GBMap%GBtypes(current)%eps_ratio = eps_ratio
    GBMap%GBtypes(current)%dw        = dw
    GBMap%GBtypes(current)%isLJ      = .false.
    
    return
  end subroutine newGBtype
  
  subroutine complete_GB_FF(status)
    integer :: status
    integer :: i, j, l, m, lm, function_type
    real(kind=dp) :: thisDP, sigma
    integer :: alloc_stat, iTheta, iPhi, nSteps, nAtypes, myATID, current
    logical :: thisProperty
    
    status = 0
    if (GBMap%currentGBtype == 0) then
       call handleError("complete_GB_FF", "No members in GBMap")
       status = -1
       return
    end if
    
    nAtypes = getSize(atypes)
    
    if (nAtypes == 0) then
       status = -1
       return
    end if
    
    ! atypes comes from c side
    do i = 1, nAtypes
       
       myATID = getFirstMatchingElement(atypes, 'c_ident', i)
       call getElementProperty(atypes, myATID, "is_LennardJones", thisProperty)
       
       if (thisProperty) then
          GBMap%currentGBtype = GBMap%currentGBtype + 1
          current = GBMap%currentGBtype
          
          GBMap%atidToGBtype(myATID) = current
          GBMap%GBtypes(current)%atid      = myATID       
          GBMap%GBtypes(current)%isLJ      = .true.          
          GBMap%GBtypes(current)%d         = getSigma(myATID) / sqrt(2.0_dp)
          GBMap%GBtypes(current)%l         = GBMap%GBtypes(current)%d
          GBMap%GBtypes(current)%eps       = getEpsilon(myATID)
          GBMap%GBtypes(current)%eps_ratio = 1.0_dp
          GBMap%GBtypes(current)%dw        = 1.0_dp
          
       endif
       
    end do
    
    haveGBMap = .true.

    
  end subroutine complete_GB_FF

  subroutine createGBMixingMap()
    integer :: nGBtypes, i, j
    real (kind = dp) :: d1, l1, e1, er1, dw1
    real (kind = dp) :: d2, l2, e2, er2, dw2
    real (kind = dp) :: er, ermu, xp, ap2

    if (GBMap%currentGBtype == 0) then
       call handleError("GB", "No members in GBMap")
       return
    end if
    
    nGBtypes = GBMap%nGBtypes

    if (.not. allocated(GBMixingMap)) then
       allocate(GBMixingMap(nGBtypes, nGBtypes))
    endif

    do i = 1, nGBtypes

       d1 = GBMap%GBtypes(i)%d
       l1 = GBMap%GBtypes(i)%l
       e1 = GBMap%GBtypes(i)%eps
       er1 = GBMap%GBtypes(i)%eps_ratio
       dw1 = GBMap%GBtypes(i)%dw

       do j = 1, nGBtypes

          d2 = GBMap%GBtypes(j)%d
          l2 = GBMap%GBtypes(j)%l
          e2 = GBMap%GBtypes(j)%eps
          er2 = GBMap%GBtypes(j)%eps_ratio
          dw2 = GBMap%GBtypes(j)%dw

!  Cleaver paper uses sqrt of squares to get sigma0 for
!  mixed interactions.
            
          GBMixingMap(i,j)%sigma0 = sqrt(d1*d1 + d2*d2)
          GBMixingMap(i,j)%xa2 = (l1*l1 - d1*d1)/(l1*l1 + d2*d2)
          GBMixingMap(i,j)%xai2 = (l2*l2 - d2*d2)/(l2*l2 + d1*d1)
          GBMixingMap(i,j)%x2 = (l1*l1 - d1*d1) * (l2*l2 - d2*d2) / &
               ((l2*l2 + d1*d1) * (l1*l1 + d2*d2))

          ! assumed LB mixing rules for now:

          GBMixingMap(i,j)%dw = 0.5_dp * (dw1 + dw2)
          GBMixingMap(i,j)%eps0 = sqrt(e1 * e2)

          er = sqrt(er1 * er2)
          ermu = er**(1.0_dp / mu)
          xp = (1.0_dp - ermu) / (1.0_dp + ermu)
          ap2 = 1.0_dp / (1.0_dp + ermu)

          GBMixingMap(i,j)%xp2 = xp*xp
          GBMixingMap(i,j)%xpap2 = xp*ap2
          GBMixingMap(i,j)%xpapi2 = xp/ap2
       enddo
    enddo
    haveMixingMap = .true.
    mu = getGayBerneMu()
    nu = getGayBerneNu()    
  end subroutine createGBMixingMap
  

  !! gay berne cutoff should be a parameter in globals, this is a temporary 
  !! work around - this should be fixed when gay berne is up and running

  function getGayBerneCut(atomID) result(cutValue)
    integer, intent(in) :: atomID 
    integer :: gbt1
    real(kind=dp) :: cutValue, l, d

    if (GBMap%currentGBtype == 0) then
       call handleError("GB", "No members in GBMap")
       return
    end if

    gbt1 = GBMap%atidToGBtype(atomID)
    l = GBMap%GBtypes(gbt1)%l
    d = GBMap%GBtypes(gbt1)%d   

    ! sigma is actually sqrt(2)*l  for prolate ellipsoids
    
    cutValue = 2.5_dp*sqrt(2.0_dp)*max(l,d)

  end function getGayBerneCut

  subroutine do_gb_pair(atid1, atid2, d, r, r2, sw, vdwMult, vpair, fpair, &
       pot, A1, A2, f1, t1, t2)
    
    integer, intent(in) :: atid1, atid2
    integer :: gbt1, gbt2, id1, id2
    real (kind=dp), intent(inout) :: r, r2, vdwMult
    real (kind=dp), dimension(3), intent(in) :: d
    real (kind=dp), dimension(3), intent(inout) :: fpair
    real (kind=dp) :: pot, sw, vpair
    real (kind=dp), dimension(9) :: A1, A2
    real (kind=dp), dimension(3) :: f1
    real (kind=dp), dimension(3) :: t1, t2
    real (kind = dp), dimension(3) :: ul1, ul2, rxu1, rxu2, uxu, rhat

    real (kind = dp) :: sigma0, dw, eps0, x2, xa2, xai2, xp2, xpap2, xpapi2
    real (kind = dp) :: e1, e2, eps, sigma, s3, s03, au2, bu2, au, bu, a, b, g, g2
    real (kind = dp) :: U, BigR, R3, R6, R7, R12, R13, H, Hp, fx, fy, fz
    real (kind = dp) :: dUdr, dUda, dUdb, dUdg, pref1, pref2
    logical :: i_is_lj, j_is_lj

    if (.not.haveMixingMap) then
       call createGBMixingMap()
    endif

    gbt1 = GBMap%atidToGBtype(atid1)
    gbt2 = GBMap%atidToGBtype(atid2)    

    i_is_LJ = GBMap%GBTypes(gbt1)%isLJ
    j_is_LJ = GBMap%GBTypes(gbt2)%isLJ

    sigma0 = GBMixingMap(gbt1, gbt2)%sigma0
    dw     = GBMixingMap(gbt1, gbt2)%dw    
    eps0   = GBMixingMap(gbt1, gbt2)%eps0  
    x2     = GBMixingMap(gbt1, gbt2)%x2    
    xa2    = GBMixingMap(gbt1, gbt2)%xa2   
    xai2   = GBMixingMap(gbt1, gbt2)%xai2  
    xp2    = GBMixingMap(gbt1, gbt2)%xp2   
    xpap2  = GBMixingMap(gbt1, gbt2)%xpap2 
    xpapi2 = GBMixingMap(gbt1, gbt2)%xpapi2
    
    ul1(1) = A1(7)
    ul1(2) = A1(8)
    ul1(3) = A1(9)

    ul2(1) = A2(7)
    ul2(2) = A2(8)
    ul2(3) = A2(9)
    
    if (i_is_LJ) then
       a = 0.0_dp
       ul1 = 0.0_dp
    else
       a = d(1)*ul1(1)   + d(2)*ul1(2)   + d(3)*ul1(3)
    endif

    if (j_is_LJ) then
       b = 0.0_dp
       ul2 = 0.0_dp
    else       
       b = d(1)*ul2(1)   + d(2)*ul2(2)   + d(3)*ul2(3)
    endif

    if (i_is_LJ.or.j_is_LJ) then
       g = 0.0_dp
    else
       g = ul1(1)*ul2(1) + ul1(2)*ul2(2) + ul1(3)*ul2(3)
    endif

    au = a / r
    bu = b / r

    au2 = au * au
    bu2 = bu * bu
    g2 = g * g

    H  = (xa2 * au2 + xai2 * bu2 - 2.0_dp*x2*au*bu*g)  / (1.0_dp - x2*g2)
    Hp = (xpap2*au2 + xpapi2*bu2 - 2.0_dp*xp2*au*bu*g) / (1.0_dp - xp2*g2)

    sigma = sigma0 / sqrt(1.0_dp - H)
    e1 = 1.0_dp / sqrt(1.0_dp - x2*g2)
    e2 = 1.0_dp - Hp
    eps = eps0 * (e1**nu) * (e2**mu)
    BigR = dw*sigma0 / (r - sigma + dw*sigma0)
    
    R3 = BigR*BigR*BigR
    R6 = R3*R3
    R7 = R6 * BigR
    R12 = R6*R6
    R13 = R6*R7

    U = vdwMult * 4.0_dp * eps * (R12 - R6)

    s3 = sigma*sigma*sigma
    s03 = sigma0*sigma0*sigma0

    pref1 = - vdwMult * 8.0_dp * eps * mu * (R12 - R6) / (e2 * r)

    pref2 = vdwMult * 8.0_dp * eps * s3 * (6.0_dp*R13 - 3.0_dp*R7) / (dw*r*s03)

    dUdr = - (pref1 * Hp + pref2 * (sigma0*sigma0*r/s3 + H))
    
    dUda = pref1 * (xpap2*au - xp2*bu*g) / (1.0_dp - xp2 * g2) &
         + pref2 * (xa2 * au - x2 *bu*g) / (1.0_dp - x2 * g2)

    dUdb = pref1 * (xpapi2*bu - xp2*au*g) / (1.0_dp - xp2 * g2) &
         + pref2 * (xai2 * bu - x2 *au*g) / (1.0_dp - x2 * g2)

    dUdg = 4.0_dp * eps * nu * (R12 - R6) * x2 * g / (1.0_dp - x2*g2) &
         + 8.0_dp * eps * mu * (R12 - R6) * (xp2*au*bu - Hp*xp2*g) / &
         (1.0_dp - xp2 * g2) / e2 &
         + 8.0_dp * eps * s3 * (3.0_dp * R7 - 6.0_dp * R13) * &  
         (x2 * au * bu - H * x2 * g) / (1.0_dp - x2 * g2) / (dw * s03)

            
    rhat = d / r

    fx = dUdr * rhat(1) + dUda * ul1(1) + dUdb * ul2(1)
    fy = dUdr * rhat(2) + dUda * ul1(2) + dUdb * ul2(2)
    fz = dUdr * rhat(3) + dUda * ul1(3) + dUdb * ul2(3)    

    rxu1 = cross_product(d, ul1)
    rxu2 = cross_product(d, ul2)    
    uxu = cross_product(ul1, ul2)
    
!!$    write(*,*) 'pref = ' , pref1, pref2
!!$    write(*,*) 'rxu1 = ' , rxu1(1), rxu1(2), rxu1(3)
!!$    write(*,*) 'rxu2 = ' , rxu2(1), rxu2(2), rxu2(3)
!!$    write(*,*) 'uxu = ' , uxu(1), uxu(2), uxu(3)
!!$    write(*,*) 'dUda = ', dUda, dudb, dudg, dudr
!!$    write(*,*) 'H = ', H,hp,sigma, e1, e2, BigR
!!$    write(*,*) 'chi = ', xa2, xai2, x2 
!!$    write(*,*) 'chip = ', xpap2, xpapi2, xp2
!!$    write(*,*) 'eps = ', eps0, e1, e2, eps
!!$    write(*,*) 'U =', U, pref1, pref2
!!$    write(*,*) 'f =', fx, fy, fz
!!$    write(*,*) 'au =', au, bu, g
!!$    

    pot = pot + U*sw

    f1(1) = f1(1) + fx*sw
    f1(2) = f1(2) + fy*sw
    f1(3) = f1(3) + fz*sw

    t1(1) = t1(1) + (dUda*rxu1(1) - dUdg*uxu(1))*sw
    t1(2) = t1(2) + (dUda*rxu1(2) - dUdg*uxu(2))*sw
    t1(3) = t1(3) + (dUda*rxu1(3) - dUdg*uxu(3))*sw

    t2(1) = t2(1) + (dUdb*rxu2(1) + dUdg*uxu(1))*sw
    t2(2) = t2(2) + (dUdb*rxu2(2) + dUdg*uxu(2))*sw
    t2(3) = t2(3) + (dUdb*rxu2(3) + dUdg*uxu(3))*sw

    vpair = vpair + U
    
    return
  end subroutine do_gb_pair
  
  subroutine destroyGBTypes()

    GBMap%nGBtypes = 0
    GBMap%currentGBtype = 0
    
    if (associated(GBMap%GBtypes)) then
       deallocate(GBMap%GBtypes)
       GBMap%GBtypes => null()
    end if
    
    if (associated(GBMap%atidToGBtype)) then
       deallocate(GBMap%atidToGBtype)
       GBMap%atidToGBtype => null()
    end if
    
    haveMixingMap = .false.
    
  end subroutine destroyGBTypes

end module gayberne
    
