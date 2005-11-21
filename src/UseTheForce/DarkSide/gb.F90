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


module gayberne
  use force_globals
  use definitions
  use simulation
  use atype_module
  use vector_class
  use status
  use lj
#ifdef IS_MPI
  use mpiSimulation
#endif
  
  implicit none

  private

#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"

  public :: newGBtype
  public :: do_gb_pair
  public :: do_gb_lj_pair
  public :: getGayBerneCut
  public :: destroyGBtypes

  type :: GBtype
     integer          :: atid
     real(kind = dp ) :: sigma
     real(kind = dp ) :: l2b_ratio 
     real(kind = dp ) :: eps
     real(kind = dp ) :: eps_ratio 
     real(kind = dp ) :: mu
     real(kind = dp ) :: nu
     real(kind = dp ) :: sigma_l
     real(kind = dp ) :: eps_l 
  end type GBtype

  type, private :: GBList
     integer               :: nGBtypes = 0
     integer               :: currentGBtype = 0
     type(GBtype), pointer :: GBtypes(:)      => null()
     integer, pointer      :: atidToGBtype(:) => null()
  end type GBList

  type(GBList), save :: GBMap

contains

  subroutine newGBtype(c_ident, sigma, l2b_ratio, eps, eps_ratio, mu, nu, &
       status)
    
    integer, intent(in) :: c_ident
    real( kind = dp ), intent(in) :: sigma, l2b_ratio, eps, eps_ratio
    real( kind = dp ), intent(in) :: mu, nu
    integer, intent(out) :: status

    integer :: nGBTypes, ntypes, myATID
    integer, pointer :: MatchList(:) => null()
    integer :: current, i
    status = 0

    if (.not.associated(GBMap%GBtypes)) then
                               
       call getMatchingElementList(atypes, "is_GayBerne", .true., &
            nGBtypes, MatchList)
       
       GBMap%nGBtypes = nGBtypes

       allocate(GBMap%GBtypes(nGBtypes))

       ntypes = getSize(atypes)
       
       allocate(GBMap%atidToGBtype(ntypes))
       
       !! initialize atidToGBtype to -1 so that we can use this
       !! array to figure out which atom comes first in the GBLJ 
       !! routine

       do i = 1, ntypes
          GBMap%atidToGBtype(i) = -1
       enddo

    endif

    GBMap%currentGBtype = GBMap%currentGBtype + 1
    current = GBMap%currentGBtype

    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)
    GBMap%atidToGBtype(myATID)        = current
    GBMap%GBtypes(current)%atid       = myATID
    GBMap%GBtypes(current)%sigma      = sigma
    GBMap%GBtypes(current)%l2b_ratio  = l2b_ratio
    GBMap%GBtypes(current)%eps        = eps
    GBMap%GBtypes(current)%eps_ratio  = eps_ratio
    GBMap%GBtypes(current)%mu         = mu
    GBMap%GBtypes(current)%nu         = nu
    GBMap%GBtypes(current)%sigma_l    = sigma*l2b_ratio
    GBMap%GBtypes(current)%eps_l      = eps*eps_ratio

    return
  end subroutine newGBtype

  
  !! gay berne cutoff should be a parameter in globals, this is a temporary 
  !! work around - this should be fixed when gay berne is up and running

  function getGayBerneCut(atomID) result(cutValue)
    integer, intent(in) :: atomID 
    integer :: gbt1
    real(kind=dp) :: cutValue, sigma, l2b_ratio

    if (GBMap%currentGBtype == 0) then
       call handleError("GB", "No members in GBMap")
       return
    end if

    gbt1 = GBMap%atidToGBtype(atomID)
    sigma = GBMap%GBtypes(gbt1)%sigma
    l2b_ratio = GBMap%GBtypes(gbt1)%l2b_ratio

    cutValue = l2b_ratio*sigma*2.5_dp
  end function getGayBerneCut

  subroutine do_gb_pair(atom1, atom2, d, r, r2, sw, vpair, fpair, &
       pot, A, f, t, do_pot)
    
    integer, intent(in) :: atom1, atom2
    integer :: atid1, atid2, gbt1, gbt2, id1, id2
    real (kind=dp), intent(inout) :: r, r2
    real (kind=dp), dimension(3), intent(in) :: d
    real (kind=dp), dimension(3), intent(inout) :: fpair
    real (kind=dp) :: pot, sw, vpair
    real (kind=dp), dimension(9,nLocal) :: A
    real (kind=dp), dimension(3,nLocal) :: f
    real (kind=dp), dimension(3,nLocal) :: t
    logical, intent(in) :: do_pot
    real (kind = dp), dimension(3) :: ul1
    real (kind = dp), dimension(3) :: ul2

    real(kind=dp) :: sigma, l2b_ratio, epsilon, eps_ratio, mu, nu, sigma_l, eps_l
    real(kind=dp) :: chi, chiprime, emu, s2
    real(kind=dp) :: r4, rdotu1, rdotu2, u1dotu2, g, gp, gpi, gmu, gmum
    real(kind=dp) :: curlyE, enu, enum, eps, dotsum, dotdiff, ds2, dd2
    real(kind=dp) :: opXdot, omXdot, opXpdot, omXpdot, pref, gfact
    real(kind=dp) :: BigR, Ri, Ri2, Ri6, Ri7, Ri12, Ri13, R126, R137
    real(kind=dp) :: dru1dx, dru1dy, dru1dz
    real(kind=dp) :: dru2dx, dru2dy, dru2dz
    real(kind=dp) :: dBigRdx, dBigRdy, dBigRdz
    real(kind=dp) :: dBigRdu1x, dBigRdu1y, dBigRdu1z
    real(kind=dp) :: dBigRdu2x, dBigRdu2y, dBigRdu2z
    real(kind=dp) :: dUdx, dUdy, dUdz
    real(kind=dp) :: dUdu1x, dUdu1y, dUdu1z, dUdu2x, dUdu2y, dUdu2z
    real(kind=dp) :: dcE, dcEdu1x, dcEdu1y, dcEdu1z, dcEdu2x, dcEdu2y, dcEdu2z
    real(kind=dp) :: depsdu1x, depsdu1y, depsdu1z, depsdu2x, depsdu2y, depsdu2z
    real(kind=dp) :: drdx, drdy, drdz
    real(kind=dp) :: dgdx, dgdy, dgdz
    real(kind=dp) :: dgdu1x, dgdu1y, dgdu1z, dgdu2x, dgdu2y, dgdu2z
    real(kind=dp) :: dgpdx, dgpdy, dgpdz
    real(kind=dp) :: dgpdu1x, dgpdu1y, dgpdu1z, dgpdu2x, dgpdu2y, dgpdu2z
    real(kind=dp) :: line1a, line1bx, line1by, line1bz
    real(kind=dp) :: line2a, line2bx, line2by, line2bz
    real(kind=dp) :: line3a, line3b, line3, line3x, line3y, line3z
    real(kind=dp) :: term1x, term1y, term1z, term1u1x, term1u1y, term1u1z
    real(kind=dp) :: term1u2x, term1u2y, term1u2z
    real(kind=dp) :: term2a, term2b, term2u1x, term2u1y, term2u1z
    real(kind=dp) :: term2u2x, term2u2y, term2u2z
    real(kind=dp) :: yick1, yick2, mess1, mess2
    
#ifdef IS_MPI
    atid1 = atid_Row(atom1)
    atid2 = atid_Col(atom2)
#else
    atid1 = atid(atom1)
    atid2 = atid(atom2)
#endif

    gbt1 = GBMap%atidToGBtype(atid1)
    gbt2 = GBMap%atidToGBtype(atid2)

    if (gbt1 .eq. gbt2) then
       sigma     = GBMap%GBtypes(gbt1)%sigma      
       l2b_ratio = GBMap%GBtypes(gbt1)%l2b_ratio  
       epsilon   = GBMap%GBtypes(gbt1)%eps        
       eps_ratio = GBMap%GBtypes(gbt1)%eps_ratio  
       mu        = GBMap%GBtypes(gbt1)%mu         
       nu        = GBMap%GBtypes(gbt1)%nu         
       sigma_l   = GBMap%GBtypes(gbt1)%sigma_l    
       eps_l     = GBMap%GBtypes(gbt1)%eps_l      
    else
       call handleError("GB", "GB-pair was called with two different GB types! OOPSE can only handle interactions for one GB type at a time.")
    endif

    s2 = (l2b_ratio)**2
    emu = (eps_ratio)**(1.0d0/mu)

    chi = (s2 - 1.0d0)/(s2 + 1.0d0)
    chiprime = (1.0d0 - emu)/(1.0d0 + emu)

    r4 = r2*r2

#ifdef IS_MPI
    ul1(1) = A_Row(7,atom1)
    ul1(2) = A_Row(8,atom1)
    ul1(3) = A_Row(9,atom1)

    ul2(1) = A_Col(7,atom2)
    ul2(2) = A_Col(8,atom2)
    ul2(3) = A_Col(9,atom2)
#else
    ul1(1) = A(7,atom1)
    ul1(2) = A(8,atom1)
    ul1(3) = A(9,atom1)

    ul2(1) = A(7,atom2)
    ul2(2) = A(8,atom2)
    ul2(3) = A(9,atom2)
#endif
    
    dru1dx = ul1(1)
    dru2dx = ul2(1)
    dru1dy = ul1(2)
    dru2dy = ul2(2)
    dru1dz = ul1(3)
    dru2dz = ul2(3)
        
    drdx = d(1) / r
    drdy = d(2) / r
    drdz = d(3) / r
    
    ! do some dot products:
    ! NB the r in these dot products is the actual intermolecular vector,
    ! and is not the unit vector in that direction.
    
    rdotu1 = d(1)*ul1(1) + d(2)*ul1(2) + d(3)*ul1(3)
    rdotu2 = d(1)*ul2(1) + d(2)*ul2(2) + d(3)*ul2(3)
    u1dotu2 = ul1(1)*ul2(1) + ul1(2)*ul2(2) +  ul1(3)*ul2(3)

    ! This stuff is all for the calculation of g(Chi) and dgdx
    ! Line numbers roughly follow the lines in equation A25 of Luckhurst 
    !   et al. Liquid Crystals 8, 451-464 (1990).
    ! We note however, that there are some major typos in that Appendix
    ! of the Luckhurst paper, particularly in equations A23, A29 and A31
    ! We have attempted to correct them below.
    
    dotsum = rdotu1+rdotu2
    dotdiff = rdotu1-rdotu2
    ds2 = dotsum*dotsum
    dd2 = dotdiff*dotdiff
  
    opXdot = 1.0d0 + Chi*u1dotu2
    omXdot = 1.0d0 - Chi*u1dotu2
    opXpdot = 1.0d0 + ChiPrime*u1dotu2
    omXpdot = 1.0d0 - ChiPrime*u1dotu2
  
    line1a = dotsum/opXdot
    line1bx = dru1dx + dru2dx
    line1by = dru1dy + dru2dy
    line1bz = dru1dz + dru2dz
    
    line2a = dotdiff/omXdot
    line2bx = dru1dx - dru2dx
    line2by = dru1dy - dru2dy
    line2bz = dru1dz - dru2dz
    
    term1x = -Chi*(line1a*line1bx + line2a*line2bx)/r2
    term1y = -Chi*(line1a*line1by + line2a*line2by)/r2
    term1z = -Chi*(line1a*line1bz + line2a*line2bz)/r2
    
    line3a = ds2/opXdot
    line3b = dd2/omXdot
    line3 = Chi*(line3a + line3b)/r4
    line3x = d(1)*line3
    line3y = d(2)*line3
    line3z = d(3)*line3
    
    dgdx = term1x + line3x
    dgdy = term1y + line3y
    dgdz = term1z + line3z

    term1u1x = 2.0d0*(line1a+line2a)*d(1)
    term1u1y = 2.0d0*(line1a+line2a)*d(2)
    term1u1z = 2.0d0*(line1a+line2a)*d(3)
    term1u2x = 2.0d0*(line1a-line2a)*d(1)
    term1u2y = 2.0d0*(line1a-line2a)*d(2)
    term1u2z = 2.0d0*(line1a-line2a)*d(3)
    
    term2a = -line3a/opXdot
    term2b =  line3b/omXdot
    
    term2u1x = Chi*ul2(1)*(term2a + term2b)
    term2u1y = Chi*ul2(2)*(term2a + term2b)
    term2u1z = Chi*ul2(3)*(term2a + term2b)
    term2u2x = Chi*ul1(1)*(term2a + term2b)
    term2u2y = Chi*ul1(2)*(term2a + term2b)
    term2u2z = Chi*ul1(3)*(term2a + term2b)
    
    pref = -Chi*0.5d0/r2

    dgdu1x = pref*(term1u1x+term2u1x)
    dgdu1y = pref*(term1u1y+term2u1y)
    dgdu1z = pref*(term1u1z+term2u1z)
    dgdu2x = pref*(term1u2x+term2u2x)
    dgdu2y = pref*(term1u2y+term2u2y)
    dgdu2z = pref*(term1u2z+term2u2z)

    g = 1.0d0 - Chi*(line3a + line3b)/(2.0d0*r2)
  
    BigR = (r - sigma*(g**(-0.5d0)) + sigma)/sigma
    Ri = 1.0d0/BigR
    Ri2 = Ri*Ri
    Ri6 = Ri2*Ri2*Ri2
    Ri7 = Ri6*Ri
    Ri12 = Ri6*Ri6
    Ri13 = Ri6*Ri7

    gfact = (g**(-1.5d0))*0.5d0

    dBigRdx = drdx/sigma + dgdx*gfact
    dBigRdy = drdy/sigma + dgdy*gfact
    dBigRdz = drdz/sigma + dgdz*gfact

    dBigRdu1x = dgdu1x*gfact
    dBigRdu1y = dgdu1y*gfact
    dBigRdu1z = dgdu1z*gfact
    dBigRdu2x = dgdu2x*gfact
    dBigRdu2y = dgdu2y*gfact
    dBigRdu2z = dgdu2z*gfact

    ! Now, we must do it again for g(ChiPrime) and dgpdx

    line1a = dotsum/opXpdot
    line2a = dotdiff/omXpdot
    term1x = -ChiPrime*(line1a*line1bx + line2a*line2bx)/r2
    term1y = -ChiPrime*(line1a*line1by + line2a*line2by)/r2
    term1z = -ChiPrime*(line1a*line1bz + line2a*line2bz)/r2
    line3a = ds2/opXpdot
    line3b = dd2/omXpdot
    line3 = ChiPrime*(line3a + line3b)/r4
    line3x = d(1)*line3
    line3y = d(2)*line3
    line3z = d(3)*line3
    
    dgpdx = term1x + line3x
    dgpdy = term1y + line3y
    dgpdz = term1z + line3z
    
    term1u1x = 2.00d0*(line1a+line2a)*d(1)
    term1u1y = 2.00d0*(line1a+line2a)*d(2) 
    term1u1z = 2.00d0*(line1a+line2a)*d(3) 
    term1u2x = 2.0d0*(line1a-line2a)*d(1)
    term1u2y = 2.0d0*(line1a-line2a)*d(2)
    term1u2z = 2.0d0*(line1a-line2a)*d(3)

    term2a = -line3a/opXpdot
    term2b =  line3b/omXpdot
    
    term2u1x = ChiPrime*ul2(1)*(term2a + term2b)
    term2u1y = ChiPrime*ul2(2)*(term2a + term2b)
    term2u1z = ChiPrime*ul2(3)*(term2a + term2b)
    term2u2x = ChiPrime*ul1(1)*(term2a + term2b)
    term2u2y = ChiPrime*ul1(2)*(term2a + term2b)
    term2u2z = ChiPrime*ul1(3)*(term2a + term2b)
  
    pref = -ChiPrime*0.5d0/r2
    
    dgpdu1x = pref*(term1u1x+term2u1x)
    dgpdu1y = pref*(term1u1y+term2u1y)
    dgpdu1z = pref*(term1u1z+term2u1z)
    dgpdu2x = pref*(term1u2x+term2u2x)
    dgpdu2y = pref*(term1u2y+term2u2y)
    dgpdu2z = pref*(term1u2z+term2u2z)
    
    gp = 1.0d0 - ChiPrime*(line3a + line3b)/(2.0d0*r2)
    gmu = gp**mu
    gpi = 1.0d0 / gp
    gmum = gmu*gpi

    curlyE = 1.0d0/dsqrt(1.0d0 - Chi*Chi*u1dotu2*u1dotu2)
    dcE = (curlyE**3)*Chi*Chi*u1dotu2

    dcEdu1x = dcE*ul2(1)
    dcEdu1y = dcE*ul2(2)
    dcEdu1z = dcE*ul2(3)
    dcEdu2x = dcE*ul1(1)
    dcEdu2y = dcE*ul1(2)
    dcEdu2z = dcE*ul1(3)
    
    enu = curlyE**nu
    enum = enu/curlyE
  
    eps = epsilon*enu*gmu

    yick1 = epsilon*enu*mu*gmum
    yick2 = epsilon*gmu*nu*enum

    depsdu1x = yick1*dgpdu1x + yick2*dcEdu1x
    depsdu1y = yick1*dgpdu1y + yick2*dcEdu1y
    depsdu1z = yick1*dgpdu1z + yick2*dcEdu1z
    depsdu2x = yick1*dgpdu2x + yick2*dcEdu2x
    depsdu2y = yick1*dgpdu2y + yick2*dcEdu2y
    depsdu2z = yick1*dgpdu2z + yick2*dcEdu2z
    
    R126 = Ri12 - Ri6
    R137 = 6.0d0*Ri7 - 12.0d0*Ri13
    
    mess1 = gmu*R137
    mess2 = R126*mu*gmum
    
    dUdx = 4.0d0*epsilon*enu*(mess1*dBigRdx + mess2*dgpdx)*sw
    dUdy = 4.0d0*epsilon*enu*(mess1*dBigRdy + mess2*dgpdy)*sw
    dUdz = 4.0d0*epsilon*enu*(mess1*dBigRdz + mess2*dgpdz)*sw
    
    dUdu1x = 4.0d0*(R126*depsdu1x + eps*R137*dBigRdu1x)*sw
    dUdu1y = 4.0d0*(R126*depsdu1y + eps*R137*dBigRdu1y)*sw
    dUdu1z = 4.0d0*(R126*depsdu1z + eps*R137*dBigRdu1z)*sw
    dUdu2x = 4.0d0*(R126*depsdu2x + eps*R137*dBigRdu2x)*sw
    dUdu2y = 4.0d0*(R126*depsdu2y + eps*R137*dBigRdu2y)*sw
    dUdu2z = 4.0d0*(R126*depsdu2z + eps*R137*dBigRdu2z)*sw
       
#ifdef IS_MPI
    f_Row(1,atom1) = f_Row(1,atom1) + dUdx
    f_Row(2,atom1) = f_Row(2,atom1) + dUdy
    f_Row(3,atom1) = f_Row(3,atom1) + dUdz
    
    f_Col(1,atom2) = f_Col(1,atom2) - dUdx
    f_Col(2,atom2) = f_Col(2,atom2) - dUdy
    f_Col(3,atom2) = f_Col(3,atom2) - dUdz
    
    t_Row(1,atom1) = t_Row(1,atom1) + ul1(3)*dUdu1y - ul1(2)*dUdu1z 
    t_Row(2,atom1) = t_Row(2,atom1) + ul1(1)*dUdu1z - ul1(3)*dUdu1x 
    t_Row(3,atom1) = t_Row(3,atom1) + ul1(2)*dUdu1x - ul1(1)*dUdu1y 
    
    t_Col(1,atom2) = t_Col(1,atom2) + ul2(3)*dUdu2y - ul2(2)*dUdu2z
    t_Col(2,atom2) = t_Col(2,atom2) + ul2(1)*dUdu2z - ul2(3)*dUdu2x 
    t_Col(3,atom2) = t_Col(3,atom2) + ul2(2)*dUdu2x - ul2(1)*dUdu2y
#else
    f(1,atom1) = f(1,atom1) + dUdx
    f(2,atom1) = f(2,atom1) + dUdy
    f(3,atom1) = f(3,atom1) + dUdz
    
    f(1,atom2) = f(1,atom2) - dUdx
    f(2,atom2) = f(2,atom2) - dUdy
    f(3,atom2) = f(3,atom2) - dUdz
    
    t(1,atom1) = t(1,atom1) + ul1(3)*dUdu1y - ul1(2)*dUdu1z 
    t(2,atom1) = t(2,atom1) + ul1(1)*dUdu1z - ul1(3)*dUdu1x 
    t(3,atom1) = t(3,atom1) + ul1(2)*dUdu1x - ul1(1)*dUdu1y 
    
    t(1,atom2) = t(1,atom2) + ul2(3)*dUdu2y - ul2(2)*dUdu2z
    t(2,atom2) = t(2,atom2) + ul2(1)*dUdu2z - ul2(3)*dUdu2x 
    t(3,atom2) = t(3,atom2) + ul2(2)*dUdu2x - ul2(1)*dUdu2y
#endif
   
    if (do_pot) then
#ifdef IS_MPI 
       pot_row(VDW_POT,atom1) = pot_row(VDW_POT,atom1) + 2.0d0*eps*R126*sw
       pot_col(VDW_POT,atom2) = pot_col(VDW_POT,atom2) + 2.0d0*eps*R126*sw
#else
       pot = pot + 4.0*eps*R126*sw
#endif
    endif
    
    vpair = vpair + 4.0*eps*R126
#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif
    
    if (molMembershipList(id1) .ne. molMembershipList(id2)) then
       
       fpair(1) = fpair(1) + dUdx
       fpair(2) = fpair(2) + dUdy
       fpair(3) = fpair(3) + dUdz
       
    endif
    
    return
  end subroutine do_gb_pair

  subroutine do_gb_lj_pair(atom1, atom2, d, r, r2, rcut, sw, vpair, fpair, &
       pot, A, f, t, do_pot)
    
    integer, intent(in) :: atom1, atom2
    integer :: id1, id2
    real (kind=dp), intent(inout) :: r, r2, rcut
    real (kind=dp), dimension(3), intent(in) :: d
    real (kind=dp), dimension(3), intent(inout) :: fpair
    real (kind=dp) :: pot, sw, vpair
    real (kind=dp), dimension(9,nLocal) :: A
    real (kind=dp), dimension(3,nLocal) :: f
    real (kind=dp), dimension(3,nLocal) :: t
    logical, intent(in) :: do_pot
    real (kind = dp), dimension(3) :: ul
    
    real(kind=dp) :: gb_sigma, gb_eps, gb_eps_ratio, gb_mu, gb_l2b_ratio
    real(kind=dp) :: s0, l2, d2, lj2
    real(kind=dp) :: eE, eS, eab, eabf, moom, mum1
    real(kind=dp) :: dx, dy, dz, drdx, drdy, drdz, rdotu
    real(kind=dp) :: mess, sab, dsabdct, depmudct
    real(kind=dp) :: epmu, depmudx, depmudy, depmudz
    real(kind=dp) :: depmudux, depmuduy, depmuduz
    real(kind=dp) :: BigR, dBigRdx, dBigRdy, dBigRdz
    real(kind=dp) :: dBigRdux, dBigRduy, dBigRduz
    real(kind=dp) :: dUdx, dUdy, dUdz, dUdux, dUduy, dUduz, e0
    real(kind=dp) :: Ri, Ri3, Ri6, Ri7, Ri12, Ri13, R126, R137, prefactor
    real(kind=dp) :: chipoalphap2, chioalpha2, ec, epsnot
    real(kind=dp) :: drdotudx, drdotudy, drdotudz    
    real(kind=dp) :: drdotudux, drdotuduy, drdotuduz    
    real(kind=dp) :: ljeps, ljsigma
    integer :: ljt1, ljt2, atid1, atid2, gbt1, gbt2
    logical :: gb_first
    
#ifdef IS_MPI
    atid1 = atid_Row(atom1)
    atid2 = atid_Col(atom2)
#else
    atid1 = atid(atom1)
    atid2 = atid(atom2)
#endif
    
    gbt1 = GBMap%atidToGBtype(atid1)
    gbt2 = GBMap%atidToGBtype(atid2)
    
    if (gbt1 .eq. -1) then
       gb_first = .false.
       if (gbt2 .eq. -1) then
          call handleError("GB", "GBLJ was called without a GB type.")
       endif
    else
       gb_first = .true.
       if (gbt2 .ne. -1) then
          call handleError("GB", "GBLJ was called with two GB types (instead of one).")
       endif
    endif
    
    ri =1/r
    
    dx = d(1)
    dy = d(2)
    dz = d(3)
    
    drdx = dx *ri
    drdy = dy *ri
    drdz = dz *ri
    
    if(gb_first)then
#ifdef IS_MPI
       ul(1) = A_Row(7,atom1)
       ul(2) = A_Row(8,atom1)
       ul(3) = A_Row(9,atom1)
#else
       ul(1) = A(7,atom1)
       ul(2) = A(8,atom1)
       ul(3) = A(9,atom1)       
#endif
       gb_sigma     = GBMap%GBtypes(gbt1)%sigma      
       gb_l2b_ratio = GBMap%GBtypes(gbt1)%l2b_ratio
       gb_eps       = GBMap%GBtypes(gbt1)%eps        
       gb_eps_ratio = GBMap%GBtypes(gbt1)%eps_ratio  
       gb_mu        = GBMap%GBtypes(gbt1)%mu         

       ljsigma = getSigma(atid2)
       ljeps = getEpsilon(atid2)
    else
#ifdef IS_MPI
       ul(1) = A_Col(7,atom2)
       ul(2) = A_Col(8,atom2)
       ul(3) = A_Col(9,atom2)
#else
       ul(1) = A(7,atom2)
       ul(2) = A(8,atom2)
       ul(3) = A(9,atom2)     
#endif
       gb_sigma     = GBMap%GBtypes(gbt2)%sigma      
       gb_l2b_ratio = GBMap%GBtypes(gbt2)%l2b_ratio
       gb_eps       = GBMap%GBtypes(gbt2)%eps        
       gb_eps_ratio = GBMap%GBtypes(gbt2)%eps_ratio  
       gb_mu        = GBMap%GBtypes(gbt2)%mu         

       ljsigma = getSigma(atid1)
       ljeps = getEpsilon(atid1)
    endif  
 
    rdotu = (dx*ul(1)+dy*ul(2)+dz*ul(3))*ri
   
    drdotudx = ul(1)*ri-rdotu*dx*ri*ri
    drdotudy = ul(2)*ri-rdotu*dy*ri*ri
    drdotudz = ul(3)*ri-rdotu*dz*ri*ri
    drdotudux = drdx
    drdotuduy = drdy
    drdotuduz = drdz

    l2 = (gb_sigma*gb_l2b_ratio)**2
    d2 = gb_sigma**2
    lj2 = ljsigma**2
    s0 = dsqrt(d2 + lj2)

    chioalpha2 = (l2 - d2)/(l2 + lj2)

    eE = dsqrt(gb_eps*gb_eps_ratio*ljeps)
    eS = dsqrt(gb_eps*ljeps)
    moom =  1.0d0 / gb_mu
    mum1 = gb_mu-1
    chipoalphap2 = 1 - (eE/eS)**moom

    !! mess matches cleaver (eq 20)
    
    mess = 1-rdotu*rdotu*chioalpha2
    sab = 1.0d0/dsqrt(mess)

    dsabdct = s0*sab*sab*sab*rdotu*chioalpha2
       
    eab = 1-chipoalphap2*rdotu*rdotu
    eabf = eS*(eab**gb_mu)

    depmudct = -2*eS*chipoalphap2*gb_mu*rdotu*(eab**mum1)
        
    BigR = (r - sab*s0 + s0)/s0
    dBigRdx = (drdx -dsabdct*drdotudx)/s0
    dBigRdy = (drdy -dsabdct*drdotudy)/s0
    dBigRdz = (drdz -dsabdct*drdotudz)/s0
    dBigRdux = (-dsabdct*drdotudux)/s0
    dBigRduy = (-dsabdct*drdotuduy)/s0
    dBigRduz = (-dsabdct*drdotuduz)/s0
    
    depmudx = depmudct*drdotudx
    depmudy = depmudct*drdotudy
    depmudz = depmudct*drdotudz
    depmudux = depmudct*drdotudux
    depmuduy = depmudct*drdotuduy
    depmuduz = depmudct*drdotuduz
    
    Ri = 1.0d0/BigR
    Ri3 = Ri*Ri*Ri
    Ri6 = Ri3*Ri3
    Ri7 = Ri6*Ri
    Ri12 = Ri6*Ri6
    Ri13 = Ri6*Ri7
    R126 = Ri12 - Ri6
    R137 = 6.0d0*Ri7 - 12.0d0*Ri13
    
    prefactor = 4.0d0
    
    dUdx = prefactor*(eabf*R137*dBigRdx + R126*depmudx)*sw
    dUdy = prefactor*(eabf*R137*dBigRdy + R126*depmudy)*sw
    dUdz = prefactor*(eabf*R137*dBigRdz + R126*depmudz)*sw

    dUdux = prefactor*(eabf*R137*dBigRdux + R126*depmudux)*sw
    dUduy = prefactor*(eabf*R137*dBigRduy + R126*depmuduy)*sw
    dUduz = prefactor*(eabf*R137*dBigRduz + R126*depmuduz)*sw
    
#ifdef IS_MPI
    f_Row(1,atom1) = f_Row(1,atom1) + dUdx
    f_Row(2,atom1) = f_Row(2,atom1) + dUdy
    f_Row(3,atom1) = f_Row(3,atom1) + dUdz
    
    f_Col(1,atom2) = f_Col(1,atom2) - dUdx
    f_Col(2,atom2) = f_Col(2,atom2) - dUdy
    f_Col(3,atom2) = f_Col(3,atom2) - dUdz
    
    if (gb_first) then
       t_Row(1,atom1) = t_Row(1,atom1) - ul(2)*dUduz + ul(3)*dUduy
       t_Row(2,atom1) = t_Row(2,atom1) - ul(3)*dUdux + ul(1)*dUduz
       t_Row(3,atom1) = t_Row(3,atom1) - ul(1)*dUduy + ul(2)*dUdux
    else
       t_Col(1,atom2) = t_Col(1,atom2) - ul(2)*dUduz + ul(3)*dUduy
       t_Col(2,atom2) = t_Col(2,atom2) - ul(3)*dUdux + ul(1)*dUduz
       t_Col(3,atom2) = t_Col(3,atom2) - ul(1)*dUduy + ul(2)*dUdux
    endif
#else    
    f(1,atom1) = f(1,atom1) + dUdx
    f(2,atom1) = f(2,atom1) + dUdy
    f(3,atom1) = f(3,atom1) + dUdz
    
    f(1,atom2) = f(1,atom2) - dUdx
    f(2,atom2) = f(2,atom2) - dUdy
    f(3,atom2) = f(3,atom2) - dUdz
    
    ! torques are cross products:    

    if (gb_first) then
       t(1,atom1) = t(1,atom1) - ul(2)*dUduz + ul(3)*dUduy
       t(2,atom1) = t(2,atom1) - ul(3)*dUdux + ul(1)*dUduz
       t(3,atom1) = t(3,atom1) - ul(1)*dUduy + ul(2)*dUdux
    else
       t(1,atom2) = t(1,atom2) - ul(2)*dUduz + ul(3)*dUduy
       t(2,atom2) = t(2,atom2) - ul(3)*dUdux + ul(1)*dUduz
       t(3,atom2) = t(3,atom2) - ul(1)*dUduy + ul(2)*dUdux
    endif

#endif
       
    if (do_pot) then
#ifdef IS_MPI 
       pot_row(VDW_POT,atom1) = pot_row(VDW_POT,atom1) + 2.0d0*eabf*R126*sw
       pot_col(VDW_POT,atom2) = pot_col(VDW_POT,atom2) + 2.0d0*eabf*R126*sw
#else
       pot = pot + prefactor*eabf*R126*sw
#endif
    endif
    
    vpair = vpair + 4.0*eabf*R126
#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif
    
    If (Molmembershiplist(Id1) .Ne. Molmembershiplist(Id2)) Then
       
       Fpair(1) = Fpair(1) + Dudx
       Fpair(2) = Fpair(2) + Dudy
       Fpair(3) = Fpair(3) + Dudz
       
    Endif
    
    return
    
  end subroutine do_gb_lj_pair

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
    
  end subroutine destroyGBTypes

end module gayberne
    
