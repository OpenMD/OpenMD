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


module gb_pair
  use force_globals
  use definitions
  use simulation
#ifdef IS_MPI
  use mpiSimulation
#endif
  
  implicit none

  PRIVATE

  logical, save :: gb_pair_initialized = .false.
  real(kind=dp), save :: gb_sigma
  real(kind=dp), save :: gb_l2b_ratio
  real(kind=dp), save :: gb_eps
  real(kind=dp), save :: gb_eps_ratio
  real(kind=dp), save :: gb_mu
  real(kind=dp), save :: gb_nu

  public :: check_gb_pair_FF
  public :: set_gb_pair_params
  public :: do_gb_pair
  public :: getGayBerneCut

contains

  subroutine check_gb_pair_FF(status)
    integer :: status 
    status = -1
    if (gb_pair_initialized) status = 0
    return
  end subroutine check_gb_pair_FF

  subroutine set_gb_pair_params(sigma, l2b_ratio, eps, eps_ratio, mu, nu)
    real( kind = dp ), intent(in) :: sigma, l2b_ratio, eps, eps_ratio
    real( kind = dp ), intent(in) :: mu, nu
   
    gb_sigma = sigma
    gb_l2b_ratio = l2b_ratio
    gb_eps = eps
    gb_eps_ratio = eps_ratio
    gb_mu = mu
    gb_nu = nu

    gb_pair_initialized = .true.
    return
  end subroutine set_gb_pair_params

  function getGayBerneCut(atomID) result(cutValue)
    integer, intent(in) :: atomID !! nah... we don't need to use this...
    real(kind=dp) :: cutValue

    cutValue = gb_sigma*2.5_dp
  end function getGayBerneCut

  subroutine do_gb_pair(atom1, atom2, d, r, r2, sw, vpair, fpair, &
       pot, A, f, t, do_pot)
    
    integer, intent(in) :: atom1, atom2
    integer :: id1, id2
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
    
    s2 = (gb_l2b_ratio)**2
    emu = (gb_eps_ratio)**(1.0d0/gb_mu)

    chi = (s2 - 1.0d0)/(s2 + 1.0d0)
    chiprime = (1.0d0 - emu)/(1.0d0 + emu)

    r4 = r2*r2

#ifdef IS_MPI
    ul1(1) = A_Row(3,atom1)
    ul1(2) = A_Row(6,atom1)
    ul1(3) = A_Row(9,atom1)

    ul2(1) = A_Col(3,atom2)
    ul2(2) = A_Col(6,atom2)
    ul2(3) = A_Col(9,atom2)
#else
    ul1(1) = A(3,atom1)
    ul1(2) = A(6,atom1)
    ul1(3) = A(9,atom1)

    ul2(1) = A(3,atom2)
    ul2(2) = A(6,atom2)
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

    term1u1x = (line1a+line2a)*dru1dx
    term1u1y = (line1a+line2a)*dru1dy
    term1u1z = (line1a+line2a)*dru1dz
    term1u2x = (line1a-line2a)*dru2dx
    term1u2y = (line1a-line2a)*dru2dy
    term1u2z = (line1a-line2a)*dru2dz
    
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
  
    BigR = (r - gb_sigma*(g**(-0.5d0)) + gb_sigma)/gb_sigma
    Ri = 1.0d0/BigR
    Ri2 = Ri*Ri
    Ri6 = Ri2*Ri2*Ri2
    Ri7 = Ri6*Ri
    Ri12 = Ri6*Ri6
    Ri13 = Ri6*Ri7

    gfact = (g**(-1.5d0))*0.5d0

    dBigRdx = drdx/gb_sigma + dgdx*gfact
    dBigRdy = drdy/gb_sigma + dgdy*gfact
    dBigRdz = drdz/gb_sigma + dgdz*gfact

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
    
    term1u1x = (line1a+line2a)*dru1dx
    term1u1y = (line1a+line2a)*dru1dy
    term1u1z = (line1a+line2a)*dru1dz
    term1u2x = (line1a-line2a)*dru2dx
    term1u2y = (line1a-line2a)*dru2dy
    term1u2z = (line1a-line2a)*dru2dz

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
    gmu = gp**gb_mu
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
    
    enu = curlyE**gb_nu
    enum = enu/curlyE
  
    eps = gb_eps*enu*gmu

    yick1 = gb_eps*enu*gb_mu*gmum
    yick2 = gb_eps*gmu*gb_nu*enum

    depsdu1x = yick1*dgpdu1x + yick2*dcEdu1x
    depsdu1y = yick1*dgpdu1y + yick2*dcEdu1y
    depsdu1z = yick1*dgpdu1z + yick2*dcEdu1z
    depsdu2x = yick1*dgpdu2x + yick2*dcEdu2x
    depsdu2y = yick1*dgpdu2y + yick2*dcEdu2y
    depsdu2z = yick1*dgpdu2z + yick2*dcEdu2z
    
    R126 = Ri12 - Ri6
    R137 = 6.0d0*Ri7 - 12.0d0*Ri13
    
    mess1 = gmu*R137
    mess2 = R126*gb_mu*gmum
    
    dUdx = 4.0d0*gb_eps*enu*(mess1*dBigRdx + mess2*dgpdx)*sw
    dUdy = 4.0d0*gb_eps*enu*(mess1*dBigRdy + mess2*dgpdy)*sw
    dUdz = 4.0d0*gb_eps*enu*(mess1*dBigRdz + mess2*dgpdz)*sw
    
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
    
    t_Row(1,atom1) = t_Row(1,atom1) - ul1(2)*dUdu1z + ul1(3)*dUdu1y
    t_Row(2,atom1) = t_Row(2,atom1) - ul1(3)*dUdu1x + ul1(1)*dUdu1z
    t_Row(3,atom1) = t_Row(3,atom1) - ul1(1)*dUdu1y + ul1(2)*dUdu1x
    
    t_Col(1,atom2) = t_Col(1,atom2) - ul2(2)*dUdu2z + ul2(3)*dUdu2y
    t_Col(2,atom2) = t_Col(2,atom2) - ul2(3)*dUdu2x + ul2(1)*dUdu2z
    t_Col(3,atom2) = t_Col(3,atom2) - ul2(1)*dUdu2y + ul2(2)*dUdu2x
#else
    f(1,atom1) = f(1,atom1) + dUdx
    f(2,atom1) = f(2,atom1) + dUdy
    f(3,atom1) = f(3,atom1) + dUdz
    
    f(1,atom2) = f(1,atom2) - dUdx
    f(2,atom2) = f(2,atom2) - dUdy
    f(3,atom2) = f(3,atom2) - dUdz
    
    t(1,atom1) = t(1,atom1) - ul1(2)*dUdu1z + ul1(3)*dUdu1y
    t(2,atom1) = t(2,atom1) - ul1(3)*dUdu1x + ul1(1)*dUdu1z
    t(3,atom1) = t(3,atom1) - ul1(1)*dUdu1y + ul1(2)*dUdu1x
    
    t(1,atom2) = t(1,atom2) - ul2(2)*dUdu2z + ul2(3)*dUdu2y
    t(2,atom2) = t(2,atom2) - ul2(3)*dUdu2x + ul2(1)*dUdu2z
    t(3,atom2) = t(3,atom2) - ul2(1)*dUdu2y + ul2(2)*dUdu2x
#endif
            
    if (do_pot) then
#ifdef IS_MPI 
       pot_row(atom1) = pot_row(atom1) + 2.0d0*eps*R126*sw
       pot_col(atom2) = pot_col(atom2) + 2.0d0*eps*R126*sw
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

end module gb_pair
