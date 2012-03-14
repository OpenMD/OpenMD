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

!! gayberne is the Gay-Berne interaction for ellipsoidal particles.  The original 
!! paper (for identical uniaxial particles) is:
!!    J. G. Gay and B. J. Berne, J. Chem. Phys., 74, 3316-3319 (1981).
!! A more-general GB potential for dissimilar uniaxial particles:
!!    D. J. Cleaver, C. M. Care, M. P. Allen and M. P. Neal, Phys. Rev. E, 
!!    54, 559-567 (1996).
!! Further parameterizations can be found in:
!!    A. P. J. Emerson, G. R. Luckhurst and S. G. Whatling, Mol. Phys., 
!!    82, 113-124 (1994).
!! And a nice force expression:
!!    G. R. Luckhurst and R. A. Stephens, Liq. Cryst. 8, 451-464 (1990).
!! Even clearer force and torque expressions:
!!    P. A. Golubkov and P. Y. Ren, J. Chem. Phys., 125, 64103 (2006).
!! New expressions for cross interactions of strength parameters:
!!    J. Wu, X. Zhen, H. Shen, G. Li, and P. Ren, J. Chem. Phys., 
!!    135, 155104 (2011).
!!
!! In this version of the GB interaction, each uniaxial ellipsoidal type 
!! is described using a set of 6 parameters:
!!  d:  range parameter for side-by-side (S) and cross (X) configurations
!!  l:  range parameter for end-to-end (E) configuration
!!  epsilon_X:  well-depth parameter for cross (X) configuration
!!  epsilon_S:  well-depth parameter for side-by-side (S) configuration
!!  epsilon_E:  well depth parameter for end-to-end (E) configuration
!!  dw: "softness" of the potential
!! 
!! Additionally, there are two "universal" paramters to govern the overall 
!! importance of the purely orientational (nu) and the mixed
!! orientational / translational (mu) parts of strength of the interactions. 
!! These parameters have default or "canonical" values, but may be changed
!! as a force field option:
!! nu_: purely orientational part : defaults to 1
!! mu_: mixed orientational / translational part : defaults to 2

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
     real(kind = dp ) :: epsX
     real(kind = dp ) :: epsS
     real(kind = dp ) :: epsE
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
  
  subroutine newGBtype(c_ident, d, l, epsX, epsS, epsE, dw, status)
    
    integer, intent(in) :: c_ident
    real( kind = dp ), intent(in) :: d, l, epsX, epsS, epsE, dw
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
    GBMap%GBtypes(current)%epsX      = epsX
    GBMap%GBtypes(current)%epsS      = epsS
    GBMap%GBtypes(current)%epsE      = epsE
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
          GBMap%GBtypes(current)%epsX      = getEpsilon(myATID)
          GBMap%GBtypes(current)%epsS      = GBMap%GBtypes(current)%epsX
          GBMap%GBtypes(current)%epsE      = GBMap%GBtypes(current)%epsX
          GBMap%GBtypes(current)%dw        = 1.0_dp
          
       endif
       
    end do
    
    haveGBMap = .true.

    
  end subroutine complete_GB_FF

  subroutine createGBMixingMap()
    integer :: nGBtypes, i, j
    real (kind = dp) :: d1, l1, eX1, eS1, eE1, dw1
    real (kind = dp) :: d2, l2, eX2, eS2, eE2, dw2
    real (kind = dp) :: xp, ap2, mi

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
       eX1 = GBMap%GBtypes(i)%epsX
       eS1 = GBMap%GBtypes(i)%epsS
       eE1 = GBMap%GBtypes(i)%epsE
       dw1 = GBMap%GBtypes(i)%dw

       do j = 1, nGBtypes

          d2 = GBMap%GBtypes(j)%d
          l2 = GBMap%GBtypes(j)%l
          eX2 = GBMap%GBtypes(j)%epsX
          eS2 = GBMap%GBtypes(j)%epsS
          eE2 = GBMap%GBtypes(j)%epsE
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
          GBMixingMap(i,j)%eps0 = sqrt(eX1 * eX2)

          mi = 1.0 / mu
      
          GBMixingMap(i,j)%xpap2  = ((eS1**mi) - (eE1**mi)) / &
               ((eS1**mi) + (eE2**mi))
          GBMixingMap(i,j)%xpapi2 = ((eS2**mi) - (eE2**mi)) / &
               ((eS2**mi) + (eE1**mi))
          GBMixingMap(i,j)%xp2    = ((eS1**mi) - (eE1**mi)) * &
               ((eS2**mi) - (eE2**mi)) / &
               ((eS2**mi) + (eE1**mi)) / ((eS1**mi) + (eE2**mi))

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
    
!!$    write(*,*) 'atypes = ',atid1, atid2
!!$    write(*,*) 'sigma0 = ',sigma0 
!!$    write(*,*) 'dw     = ',dw 
!!$    write(*,*) 'eps0   = ',eps0   
!!$    write(*,*) 'x2     = ',x2     
!!$    write(*,*) 'xa2    = ',xa2    
!!$    write(*,*) 'xai2   = ',xai2   
!!$    write(*,*) 'xp2    = ',xp2    
!!$    write(*,*) 'xpap2  = ',xpap2  
!!$    write(*,*) 'xpapi2 = ',xpapi2 


    ul1(1) = A1(7)
    ul1(2) = A1(8)
    ul1(3) = A1(9)

    ul2(1) = A2(7)
    ul2(2) = A2(8)
    ul2(3) = A2(9)
    
!!$    write(*,*) 'ul1 = ', ul1(1), ul1(2), ul1(3)
!!$    write(*,*) 'ul2 = ', ul2(1), ul2(2), ul2(3)

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

!!$    write(*,*) 'au2 = ',au2
!!$    write(*,*) 'bu2 = ',bu2
!!$    write(*,*) 'g2 = ',g2
!!$    write(*,*) 'H = ',H
!!$    write(*,*) 'Hp = ',Hp

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

!!$    write(*,*) 'vdwMult = ', vdwMult
!!$    write(*,*) 'eps = ', eps
!!$    write(*,*) 'mu = ', mu
!!$    write(*,*) 'R12 = ', R12
!!$    write(*,*) 'R6 = ', R6 
!!$    write(*,*) 'R13 = ', R13
!!$    write(*,*) 'R7 = ', R7 
!!$    write(*,*) 'e2 = ', e2 
!!$    write(*,*) 'rij = ', r
!!$    write(*,*) 's3 = ', s3 
!!$    write(*,*) 's03 = ', s03 
!!$    write(*,*) 'dw = ', dw 

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
!!$    write(*,*) 'pref = ',pref1 ,  pref2
!!$    write(*,*) 'dU = ',dUdr ,  dUda, dUdb ,  dUdg

            
    rhat = d / r

    fx = dUdr * rhat(1) + dUda * ul1(1) + dUdb * ul2(1)
    fy = dUdr * rhat(2) + dUda * ul1(2) + dUdb * ul2(2)
    fz = dUdr * rhat(3) + dUda * ul1(3) + dUdb * ul2(3)    

    rxu1 = cross_product(d, ul1)
    rxu2 = cross_product(d, ul2)    
    uxu = cross_product(ul1, ul2)

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

!!$    write(*,*) 'f1 term = ', fx*sw, fy*sw, fz*sw
!!$    write(*,*) 't1 term = ', (dUda*rxu1(1) - dUdg*uxu(1))*sw, &
!!$         (dUda*rxu1(2) - dUdg*uxu(2))*sw, &
!!$         (dUda*rxu1(3) - dUdg*uxu(3))*sw                              
!!$    write(*,*) 't2 term = ', (dUdb*rxu2(1) + dUdg*uxu(1))*sw, &
!!$         (dUdb*rxu2(2) + dUdg*uxu(2))*sw, &
!!$         (dUdb*rxu2(3) + dUdg*uxu(3))*sw
!!$    write(*,*) 'vp term = ', U

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
    
