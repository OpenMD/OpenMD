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

module electrostatic_module

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
#include "UseTheForce/DarkSide/fElectrostaticSummationMethod.h"
#include "UseTheForce/DarkSide/fElectrostaticScreeningMethod.h"


  !! these prefactors convert the multipole interactions into kcal / mol
  !! all were computed assuming distances are measured in angstroms
  !! Charge-Charge, assuming charges are measured in electrons
  real(kind=dp), parameter :: pre11 = 332.0637778_dp
  !! Charge-Dipole, assuming charges are measured in electrons, and
  !! dipoles are measured in debyes
  real(kind=dp), parameter :: pre12 = 69.13373_dp
  !! Dipole-Dipole, assuming dipoles are measured in debyes
  real(kind=dp), parameter :: pre22 = 14.39325_dp
  !! Charge-Quadrupole, assuming charges are measured in electrons, and
  !! quadrupoles are measured in 10^-26 esu cm^2
  !! This unit is also known affectionately as an esu centi-barn.
  real(kind=dp), parameter :: pre14 = 69.13373_dp

  real(kind=dp), parameter :: zero = 0.0_dp
  
  !! conversions for the simulation box dipole moment
  real(kind=dp), parameter :: chargeToC = 1.60217733e-19_dp
  real(kind=dp), parameter :: angstromToM = 1.0e-10_dp
  real(kind=dp), parameter :: debyeToCm = 3.33564095198e-30_dp

  !! number of points for electrostatic splines 
  integer, parameter :: np = 100

  !! variables to handle different summation methods for long-range 
  !! electrostatics:
  integer, save :: summationMethod = NONE
  integer, save :: screeningMethod = UNDAMPED
  logical, save :: summationMethodChecked = .false.
  real(kind=DP), save :: defaultCutoff = 0.0_DP
  real(kind=DP), save :: defaultCutoff2 = 0.0_DP
  logical, save :: haveDefaultCutoff = .false.
  real(kind=DP), save :: dampingAlpha = 0.0_DP
  real(kind=DP), save :: alpha2 = 0.0_DP 
  real(kind=DP), save :: alpha4 = 0.0_DP 
  real(kind=DP), save :: alpha6 = 0.0_DP 
  real(kind=DP), save :: alpha8 = 0.0_DP 
  logical, save :: haveDampingAlpha = .false.
  real(kind=DP), save :: dielectric = 1.0_DP
  logical, save :: haveDielectric = .false.
  real(kind=DP), save :: constEXP = 0.0_DP
  real(kind=dp), save :: rcuti = 0.0_DP
  real(kind=dp), save :: rcuti2 = 0.0_DP
  real(kind=dp), save :: rcuti3 = 0.0_DP
  real(kind=dp), save :: rcuti4 = 0.0_DP
  real(kind=dp), save :: alphaPi = 0.0_DP
  real(kind=dp), save :: invRootPi = 0.0_DP
  real(kind=dp), save :: rrf = 1.0_DP
  real(kind=dp), save :: rt = 1.0_DP
  real(kind=dp), save :: rrfsq = 1.0_DP
  real(kind=dp), save :: preRF = 0.0_DP
  real(kind=dp), save :: preRF2 = 0.0_DP
  real(kind=dp), save :: erfcVal = 1.0_DP
  real(kind=dp), save :: derfcVal = 0.0_DP
  type(cubicSpline), save :: erfcSpline
  logical, save :: haveElectroSpline = .false.
  real(kind=dp), save :: c1 = 1.0_DP
  real(kind=dp), save :: c2 = 1.0_DP
  real(kind=dp), save :: c3 = 0.0_DP
  real(kind=dp), save :: c4 = 0.0_DP
  real(kind=dp), save :: c5 = 0.0_DP
  real(kind=dp), save :: c6 = 0.0_DP
  real(kind=dp), save :: c1c = 1.0_DP
  real(kind=dp), save :: c2c = 1.0_DP
  real(kind=dp), save :: c3c = 0.0_DP
  real(kind=dp), save :: c4c = 0.0_DP
  real(kind=dp), save :: c5c = 0.0_DP
  real(kind=dp), save :: c6c = 0.0_DP
  real(kind=dp), save :: one_third = 1.0_DP / 3.0_DP

#if defined(__IFC) || defined(__PGI)
! error function for ifc version > 7.
  real(kind=dp), external :: erfc
#endif
  
  public :: setElectrostaticSummationMethod
  public :: setScreeningMethod
  public :: setElectrostaticCutoffRadius
  public :: setDampingAlpha
  public :: setReactionFieldDielectric
  public :: buildElectroSpline
  public :: newElectrostaticType
  public :: setCharge
  public :: setDipoleMoment
  public :: setSplitDipoleDistance
  public :: setQuadrupoleMoments
  public :: doElectrostaticPair
  public :: getCharge
  public :: getDipoleMoment
  public :: destroyElectrostaticTypes
  public :: self_self
  public :: rf_self_excludes
  public :: accumulate_box_dipole

  type :: Electrostatic
     integer :: c_ident
     logical :: is_Charge = .false.
     logical :: is_Dipole = .false.
     logical :: is_SplitDipole = .false.
     logical :: is_Quadrupole = .false.
     logical :: is_Tap = .false.
     real(kind=DP) :: charge = 0.0_DP
     real(kind=DP) :: dipole_moment = 0.0_DP
     real(kind=DP) :: split_dipole_distance = 0.0_DP
     real(kind=DP), dimension(3) :: quadrupole_moments = 0.0_DP
  end type Electrostatic

  type(Electrostatic), dimension(:), allocatable :: ElectrostaticMap

  logical, save :: hasElectrostaticMap

contains

  subroutine setElectrostaticSummationMethod(the_ESM)
    integer, intent(in) :: the_ESM    

    if ((the_ESM .le. 0) .or. (the_ESM .gt. REACTION_FIELD)) then
       call handleError("setElectrostaticSummationMethod", "Unsupported Summation Method")
    endif

    summationMethod = the_ESM

  end subroutine setElectrostaticSummationMethod

  subroutine setScreeningMethod(the_SM)
    integer, intent(in) :: the_SM    
    screeningMethod = the_SM
  end subroutine setScreeningMethod

  subroutine setElectrostaticCutoffRadius(thisRcut, thisRsw) 
    real(kind=dp), intent(in) :: thisRcut
    real(kind=dp), intent(in) :: thisRsw
    defaultCutoff = thisRcut
    defaultCutoff2 = defaultCutoff*defaultCutoff
    rrf = defaultCutoff
    rt = thisRsw
    haveDefaultCutoff = .true.
  end subroutine setElectrostaticCutoffRadius

  subroutine setDampingAlpha(thisAlpha)
    real(kind=dp), intent(in) :: thisAlpha
    dampingAlpha = thisAlpha
    alpha2 = dampingAlpha*dampingAlpha
    alpha4 = alpha2*alpha2
    alpha6 = alpha4*alpha2
    alpha8 = alpha4*alpha4
    haveDampingAlpha = .true.
  end subroutine setDampingAlpha
  
  subroutine setReactionFieldDielectric(thisDielectric)
    real(kind=dp), intent(in) :: thisDielectric
    dielectric = thisDielectric
    haveDielectric = .true.
  end subroutine setReactionFieldDielectric

  subroutine buildElectroSpline()
    real( kind = dp ), dimension(np) :: xvals, yvals
    real( kind = dp ) :: dx, rmin, rval
    integer :: i

    rmin = 0.0_dp

    dx = (defaultCutoff-rmin) / dble(np-1)
    
    do i = 1, np
       rval = rmin + dble(i-1)*dx
       xvals(i) = rval
       yvals(i) = erfc(dampingAlpha*rval)
    enddo

    call newSpline(erfcSpline, xvals, yvals, .true.)

    haveElectroSpline = .true.
  end subroutine buildElectroSpline

  subroutine newElectrostaticType(c_ident, is_Charge, is_Dipole, &
       is_SplitDipole, is_Quadrupole, is_Tap, status)

    integer, intent(in) :: c_ident
    logical, intent(in) :: is_Charge
    logical, intent(in) :: is_Dipole
    logical, intent(in) :: is_SplitDipole
    logical, intent(in) :: is_Quadrupole
    logical, intent(in) :: is_Tap
    integer, intent(out) :: status
    integer :: nAtypes, myATID, i, j

    status = 0
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    !! Be simple-minded and assume that we need an ElectrostaticMap that
    !! is the same size as the total number of atom types

    if (.not.allocated(ElectrostaticMap)) then

       nAtypes = getSize(atypes)

       if (nAtypes == 0) then
          status = -1
          return
       end if

       allocate(ElectrostaticMap(nAtypes))

    end if

    if (myATID .gt. size(ElectrostaticMap)) then
       status = -1
       return
    endif

    ! set the values for ElectrostaticMap for this atom type:

    ElectrostaticMap(myATID)%c_ident = c_ident
    ElectrostaticMap(myATID)%is_Charge = is_Charge
    ElectrostaticMap(myATID)%is_Dipole = is_Dipole
    ElectrostaticMap(myATID)%is_SplitDipole = is_SplitDipole
    ElectrostaticMap(myATID)%is_Quadrupole = is_Quadrupole
    ElectrostaticMap(myATID)%is_Tap = is_Tap

    hasElectrostaticMap = .true.

  end subroutine newElectrostaticType

  subroutine setCharge(c_ident, charge, status)
    integer, intent(in) :: c_ident
    real(kind=dp), intent(in) :: charge
    integer, intent(out) :: status
    integer :: myATID

    status = 0
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    if (.not.hasElectrostaticMap) then
       call handleError("electrostatic", "no ElectrostaticMap was present before first call of setCharge!")
       status = -1
       return
    end if

    if (myATID .gt. size(ElectrostaticMap)) then
       call handleError("electrostatic", "ElectrostaticMap was found to be too small during setCharge!")
       status = -1
       return
    endif

    if (.not.ElectrostaticMap(myATID)%is_Charge) then
       call handleError("electrostatic", "Attempt to setCharge of an atom type that is not a charge!")
       status = -1
       return
    endif

    ElectrostaticMap(myATID)%charge = charge
  end subroutine setCharge

  subroutine setDipoleMoment(c_ident, dipole_moment, status)
    integer, intent(in) :: c_ident
    real(kind=dp), intent(in) :: dipole_moment
    integer, intent(out) :: status
    integer :: myATID

    status = 0
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    if (.not.hasElectrostaticMap) then
       call handleError("electrostatic", "no ElectrostaticMap was present before first call of setDipoleMoment!")
       status = -1
       return
    end if

    if (myATID .gt. size(ElectrostaticMap)) then
       call handleError("electrostatic", "ElectrostaticMap was found to be too small during setDipoleMoment!")
       status = -1
       return
    endif

    if (.not.ElectrostaticMap(myATID)%is_Dipole) then
       call handleError("electrostatic", "Attempt to setDipoleMoment of an atom type that is not a dipole!")
       status = -1
       return
    endif

    ElectrostaticMap(myATID)%dipole_moment = dipole_moment
  end subroutine setDipoleMoment

  subroutine setSplitDipoleDistance(c_ident, split_dipole_distance, status)
    integer, intent(in) :: c_ident
    real(kind=dp), intent(in) :: split_dipole_distance
    integer, intent(out) :: status
    integer :: myATID

    status = 0
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    if (.not.hasElectrostaticMap) then
       call handleError("electrostatic", "no ElectrostaticMap was present before first call of setSplitDipoleDistance!")
       status = -1
       return
    end if

    if (myATID .gt. size(ElectrostaticMap)) then
       call handleError("electrostatic", "ElectrostaticMap was found to be too small during setSplitDipoleDistance!")
       status = -1
       return
    endif

    if (.not.ElectrostaticMap(myATID)%is_SplitDipole) then
       call handleError("electrostatic", "Attempt to setSplitDipoleDistance of an atom type that is not a splitDipole!")
       status = -1
       return
    endif

    ElectrostaticMap(myATID)%split_dipole_distance = split_dipole_distance
  end subroutine setSplitDipoleDistance

  subroutine setQuadrupoleMoments(c_ident, quadrupole_moments, status)
    integer, intent(in) :: c_ident
    real(kind=dp), intent(in), dimension(3) :: quadrupole_moments
    integer, intent(out) :: status
    integer :: myATID, i, j

    status = 0
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    if (.not.hasElectrostaticMap) then
       call handleError("electrostatic", "no ElectrostaticMap was present before first call of setQuadrupoleMoments!")
       status = -1
       return
    end if

    if (myATID .gt. size(ElectrostaticMap)) then
       call handleError("electrostatic", "ElectrostaticMap was found to be too small during setQuadrupoleMoments!")
       status = -1
       return
    endif

    if (.not.ElectrostaticMap(myATID)%is_Quadrupole) then
       call handleError("electrostatic", "Attempt to setQuadrupoleMoments of an atom type that is not a quadrupole!")
       status = -1
       return
    endif

    do i = 1, 3
       ElectrostaticMap(myATID)%quadrupole_moments(i) = &
            quadrupole_moments(i)
    enddo

  end subroutine setQuadrupoleMoments


  function getCharge(atid) result (c)
    integer, intent(in) :: atid
    integer :: localError
    real(kind=dp) :: c

    if (.not.hasElectrostaticMap) then
       call handleError("electrostatic", "no ElectrostaticMap was present before first call of getCharge!")
       return
    end if

    if (.not.ElectrostaticMap(atid)%is_Charge) then
       call handleError("electrostatic", "getCharge was called for an atom type that isn't a charge!")
       return
    endif

    c = ElectrostaticMap(atid)%charge
  end function getCharge

  function getDipoleMoment(atid) result (dm)
    integer, intent(in) :: atid
    integer :: localError
    real(kind=dp) :: dm

    if (.not.hasElectrostaticMap) then
       call handleError("electrostatic", "no ElectrostaticMap was present before first call of getDipoleMoment!")
       return
    end if

    if (.not.ElectrostaticMap(atid)%is_Dipole) then
       call handleError("electrostatic", "getDipoleMoment was called for an atom type that isn't a dipole!")
       return
    endif

    dm = ElectrostaticMap(atid)%dipole_moment
  end function getDipoleMoment

  subroutine checkSummationMethod()

    if (.not.haveDefaultCutoff) then
       call handleError("checkSummationMethod", "no Default Cutoff set!")
    endif

    rcuti = 1.0_dp / defaultCutoff
    rcuti2 = rcuti*rcuti
    rcuti3 = rcuti2*rcuti
    rcuti4 = rcuti2*rcuti2

    if (screeningMethod .eq. DAMPED) then
       if (.not.haveDampingAlpha) then
          call handleError("checkSummationMethod", "no Damping Alpha set!")
       endif
       
       if (.not.haveDefaultCutoff) then
          call handleError("checkSummationMethod", "no Default Cutoff set!")
       endif

       constEXP = exp(-alpha2*defaultCutoff2)
       invRootPi = 0.56418958354775628695_dp
       alphaPi = 2.0_dp*dampingAlpha*invRootPi

       c1c = erfc(dampingAlpha*defaultCutoff) * rcuti
       c2c = alphaPi*constEXP*rcuti + c1c*rcuti
       c3c = 2.0_dp*alphaPi*alpha2 + 3.0_dp*c2c*rcuti
       c4c = 4.0_dp*alphaPi*alpha4 + 5.0_dp*c3c*rcuti2
       c5c = 8.0_dp*alphaPi*alpha6 + 7.0_dp*c4c*rcuti2
       c6c = 16.0_dp*alphaPi*alpha8 + 9.0_dp*c5c*rcuti2
    else
       c1c = rcuti
       c2c = c1c*rcuti
       c3c = 3.0_dp*c2c*rcuti
       c4c = 5.0_dp*c3c*rcuti2
       c5c = 7.0_dp*c4c*rcuti2
       c6c = 9.0_dp*c5c*rcuti2
    endif

    if (summationMethod .eq. REACTION_FIELD) then
       if (haveDielectric) then
          defaultCutoff2 = defaultCutoff*defaultCutoff
          preRF = (dielectric-1.0_dp) / &
               ((2.0_dp*dielectric+1.0_dp)*defaultCutoff2*defaultCutoff)
          preRF2 = 2.0_dp*preRF
       else 
          call handleError("checkSummationMethod", "Dielectric not set")
       endif
       
    endif

    if (.not.haveElectroSpline) then
       call buildElectroSpline()
    end if

    summationMethodChecked = .true.
  end subroutine checkSummationMethod


  subroutine doElectrostaticPair(atom1, atom2, me1, me2, d, rij, r2, rcut, sw, &
       electroMult, vpair, fpair, pot, eF1, eF2, f1, t1, t2, do_pot)

    logical, intent(in) :: do_pot

    integer, intent(in) :: atom1, atom2, me1, me2
    integer :: localError

    real(kind=dp), intent(in) :: rij, r2, sw, rcut, electroMult
    real(kind=dp), intent(in), dimension(3) :: d
    real(kind=dp), intent(inout) :: vpair
    real(kind=dp), intent(inout), dimension(3) :: fpair    

    real( kind = dp ) :: pot
    real( kind = dp ), dimension(9) :: eF1, eF2  ! eFrame = electroFrame
    real( kind = dp ), dimension(3) :: f1
    real( kind = dp ), dimension(3,nLocal) :: felec
    real( kind = dp ), dimension(3) :: t1, t2

    real (kind = dp), dimension(3) :: ux_i, uy_i, uz_i
    real (kind = dp), dimension(3) :: ux_j, uy_j, uz_j
    real (kind = dp), dimension(3) :: dudux_i, duduy_i, duduz_i
    real (kind = dp), dimension(3) :: dudux_j, duduy_j, duduz_j

    logical :: i_is_Charge, i_is_Dipole, i_is_SplitDipole, i_is_Quadrupole
    logical :: j_is_Charge, j_is_Dipole, j_is_SplitDipole, j_is_Quadrupole
    logical :: i_is_Tap, j_is_Tap
    integer :: id1, id2
    real (kind=dp) :: q_i, q_j, mu_i, mu_j, d_i, d_j
    real (kind=dp) :: qxx_i, qyy_i, qzz_i
    real (kind=dp) :: qxx_j, qyy_j, qzz_j
    real (kind=dp) :: cx_i, cy_i, cz_i
    real (kind=dp) :: cx_j, cy_j, cz_j
    real (kind=dp) :: cx2, cy2, cz2
    real (kind=dp) :: ct_i, ct_j, ct_ij, a0, a1
    real (kind=dp) :: riji, ri, ri2, ri3, ri4
    real (kind=dp) :: pref, vterm, epot, dudr, vterm1, vterm2 
    real (kind=dp) :: xhat, yhat, zhat
    real (kind=dp) :: dudx, dudy, dudz
    real (kind=dp) :: scale, sc2, bigR
    real (kind=dp) :: varEXP
    real (kind=dp) :: pot_term
    real (kind=dp) :: preVal, rfVal
    real (kind=dp) :: c2ri, c3ri, c4rij
    real (kind=dp) :: cti3, ctj3, ctidotj
    real (kind=dp) :: preSw, preSwSc
    real (kind=dp) :: xhatdot2, yhatdot2, zhatdot2
    real (kind=dp) :: xhatc4, yhatc4, zhatc4

    if (.not.summationMethodChecked) then
       call checkSummationMethod()
    endif

    !! some variables we'll need independent of electrostatic type:

    riji = 1.0_dp / rij
   
    xhat = d(1) * riji
    yhat = d(2) * riji
    zhat = d(3) * riji

    !! logicals
    i_is_Charge = ElectrostaticMap(me1)%is_Charge
    i_is_Dipole = ElectrostaticMap(me1)%is_Dipole
    i_is_SplitDipole = ElectrostaticMap(me1)%is_SplitDipole
    i_is_Quadrupole = ElectrostaticMap(me1)%is_Quadrupole
    i_is_Tap = ElectrostaticMap(me1)%is_Tap

    j_is_Charge = ElectrostaticMap(me2)%is_Charge
    j_is_Dipole = ElectrostaticMap(me2)%is_Dipole
    j_is_SplitDipole = ElectrostaticMap(me2)%is_SplitDipole
    j_is_Quadrupole = ElectrostaticMap(me2)%is_Quadrupole
    j_is_Tap = ElectrostaticMap(me2)%is_Tap

    if (i_is_Charge) then
       q_i = ElectrostaticMap(me1)%charge      
    endif

    if (i_is_Dipole) then
       mu_i = ElectrostaticMap(me1)%dipole_moment

       uz_i(1) = eF1(3)
       uz_i(2) = eF1(6)
       uz_i(3) = eF1(9)

       ct_i = uz_i(1)*xhat + uz_i(2)*yhat + uz_i(3)*zhat

       if (i_is_SplitDipole) then
          d_i = ElectrostaticMap(me1)%split_dipole_distance
       endif
       duduz_i = zero
    endif

    if (i_is_Quadrupole) then
       qxx_i = ElectrostaticMap(me1)%quadrupole_moments(1)
       qyy_i = ElectrostaticMap(me1)%quadrupole_moments(2)
       qzz_i = ElectrostaticMap(me1)%quadrupole_moments(3)

       ux_i(1) = eF1(1)
       ux_i(2) = eF1(4)
       ux_i(3) = eF1(7)
       uy_i(1) = eF1(2)
       uy_i(2) = eF1(5)
       uy_i(3) = eF1(8)
       uz_i(1) = eF1(3)
       uz_i(2) = eF1(6)
       uz_i(3) = eF1(9)

       cx_i = ux_i(1)*xhat + ux_i(2)*yhat + ux_i(3)*zhat
       cy_i = uy_i(1)*xhat + uy_i(2)*yhat + uy_i(3)*zhat
       cz_i = uz_i(1)*xhat + uz_i(2)*yhat + uz_i(3)*zhat
       dudux_i = zero
       duduy_i = zero
       duduz_i = zero
    endif

    if (j_is_Charge) then
       q_j = ElectrostaticMap(me2)%charge      
    endif

    if (j_is_Dipole) then
       mu_j = ElectrostaticMap(me2)%dipole_moment

       uz_j(1) = eF2(3)
       uz_j(2) = eF2(6)
       uz_j(3) = eF2(9)

       ct_j = uz_j(1)*xhat + uz_j(2)*yhat + uz_j(3)*zhat

       if (j_is_SplitDipole) then
          d_j = ElectrostaticMap(me2)%split_dipole_distance
       endif
       duduz_j = zero
    endif

    if (j_is_Quadrupole) then
       qxx_j = ElectrostaticMap(me2)%quadrupole_moments(1)
       qyy_j = ElectrostaticMap(me2)%quadrupole_moments(2)
       qzz_j = ElectrostaticMap(me2)%quadrupole_moments(3)

       ux_j(1) = eF2(1)
       ux_j(2) = eF2(4)
       ux_j(3) = eF2(7)
       uy_j(1) = eF2(2)
       uy_j(2) = eF2(5)
       uy_j(3) = eF2(8)
       uz_j(1) = eF2(3)
       uz_j(2) = eF2(6)
       uz_j(3) = eF2(9)

       cx_j = ux_j(1)*xhat + ux_j(2)*yhat + ux_j(3)*zhat
       cy_j = uy_j(1)*xhat + uy_j(2)*yhat + uy_j(3)*zhat
       cz_j = uz_j(1)*xhat + uz_j(2)*yhat + uz_j(3)*zhat
       dudux_j = zero
       duduy_j = zero
       duduz_j = zero
    endif
   
    epot = zero
    dudx = zero
    dudy = zero
    dudz = zero  

    if (i_is_Charge) then

       if (j_is_Charge) then
          if (screeningMethod .eq. DAMPED) then
             ! assemble the damping variables
             call lookupUniformSpline1d(erfcSpline, rij, erfcVal, derfcVal)
             c1 = erfcVal*riji
             c2 = (-derfcVal + c1)*riji
          else
             c1 = riji
             c2 = c1*riji
          endif

          preVal = electroMult * pre11 * q_i * q_j

          if (summationMethod .eq. SHIFTED_POTENTIAL) then
             vterm = preVal * (c1 - c1c)
             
             dudr  = -sw * preVal * c2
  
          elseif (summationMethod .eq. SHIFTED_FORCE) then
             vterm = preVal * ( c1 - c1c + c2c*(rij - defaultCutoff) )
             
             dudr  = sw * preVal * (c2c - c2)
  
          elseif (summationMethod .eq. REACTION_FIELD) then
             rfVal = electroMult * preRF*rij*rij
             vterm = preVal * ( riji + rfVal )
             
             dudr  = sw * preVal * ( 2.0_dp*rfVal - riji )*riji
  
          else
             vterm = preVal * riji*erfcVal
             
             dudr  = - sw * preVal * c2
  
          endif

          vpair = vpair + vterm
          epot = epot + sw*vterm

          dudx = dudx + dudr * xhat
          dudy = dudy + dudr * yhat
          dudz = dudz + dudr * zhat

       endif

       if (j_is_Dipole) then
          ! pref is used by all the possible methods
          pref = electroMult * pre12 * q_i * mu_j
          preSw = sw*pref

          if (summationMethod .eq. REACTION_FIELD) then
             ri2 = riji * riji
             ri3 = ri2 * riji
    
             vterm = - pref * ct_j * ( ri2 - preRF2*rij )
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             dudx = dudx - preSw*( ri3*(uz_j(1) - 3.0_dp*ct_j*xhat) - &
                  preRF2*uz_j(1) )
             dudy = dudy - preSw*( ri3*(uz_j(2) - 3.0_dp*ct_j*yhat) - &
                  preRF2*uz_j(2) )
             dudz = dudz - preSw*( ri3*(uz_j(3) - 3.0_dp*ct_j*zhat) - &
                  preRF2*uz_j(3) )         
             duduz_j(1) = duduz_j(1) - preSw * xhat * ( ri2 - preRF2*rij )
             duduz_j(2) = duduz_j(2) - preSw * yhat * ( ri2 - preRF2*rij )
             duduz_j(3) = duduz_j(3) - preSw * zhat * ( ri2 - preRF2*rij )

          else
             ! determine the inverse r used if we have split dipoles
             if (j_is_SplitDipole) then
                BigR = sqrt(r2 + 0.25_dp * d_j * d_j)
                ri = 1.0_dp / BigR
                scale = rij * ri
             else
                ri = riji
                scale = 1.0_dp
             endif

             sc2 = scale * scale

             if (screeningMethod .eq. DAMPED) then
                ! assemble the damping variables
                call lookupUniformSpline1d(erfcSpline, rij, erfcVal, derfcVal)
                c1 = erfcVal*ri
                c2 = (-derfcVal + c1)*ri
                c3 = -2.0_dp*derfcVal*alpha2 + 3.0_dp*c2*ri
             else
                c1 = ri
                c2 = c1*ri
                c3 = 3.0_dp*c2*ri
             endif
             
             c2ri = c2*ri

             ! calculate the potential
             pot_term =  scale * c2
             vterm = -pref * ct_j * pot_term
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             ! calculate derivatives for forces and torques
             dudx = dudx - preSw*( uz_j(1)*c2ri - ct_j*xhat*sc2*c3 )
             dudy = dudy - preSw*( uz_j(2)*c2ri - ct_j*yhat*sc2*c3 )
             dudz = dudz - preSw*( uz_j(3)*c2ri - ct_j*zhat*sc2*c3 )
                          
             duduz_j(1) = duduz_j(1) - preSw * pot_term * xhat
             duduz_j(2) = duduz_j(2) - preSw * pot_term * yhat
             duduz_j(3) = duduz_j(3) - preSw * pot_term * zhat

          endif
       endif

       if (j_is_Quadrupole) then
          ! first precalculate some necessary variables
          cx2 = cx_j * cx_j
          cy2 = cy_j * cy_j
          cz2 = cz_j * cz_j
          pref =  electroMult * pre14 * q_i * one_third
          
          if (screeningMethod .eq. DAMPED) then
             ! assemble the damping variables
             call lookupUniformSpline1d(erfcSpline, rij, erfcVal, derfcVal)
             c1 = erfcVal*riji
             c2 = (-derfcVal + c1)*riji
             c3 = -2.0_dp*derfcVal*alpha2 + 3.0_dp*c2*riji
             c4 = -4.0_dp*derfcVal*alpha4 + 5.0_dp*c3*riji*riji
          else
             c1 = riji
             c2 = c1*riji
             c3 = 3.0_dp*c2*riji
             c4 = 5.0_dp*c3*riji*riji
          endif

          ! precompute variables for convenience
          preSw = sw*pref
          c2ri = c2*riji
          c3ri = c3*riji
          c4rij = c4*rij
          xhatdot2 = 2.0_dp*xhat*c3
          yhatdot2 = 2.0_dp*yhat*c3
          zhatdot2 = 2.0_dp*zhat*c3
          xhatc4 = xhat*c4rij
          yhatc4 = yhat*c4rij
          zhatc4 = zhat*c4rij

          ! calculate the potential
          pot_term = ( qxx_j*(cx2*c3 - c2ri) + qyy_j*(cy2*c3 - c2ri) + &
               qzz_j*(cz2*c3 - c2ri) )
          vterm = pref * pot_term
          vpair = vpair + vterm
          epot = epot + sw*vterm

          ! calculate derivatives for the forces and torques
          dudx = dudx - preSw * ( &
               qxx_j*(cx2*xhatc4 - (2.0_dp*cx_j*ux_j(1) + xhat)*c3ri) + &
               qyy_j*(cy2*xhatc4 - (2.0_dp*cy_j*uy_j(1) + xhat)*c3ri) + &
               qzz_j*(cz2*xhatc4 - (2.0_dp*cz_j*uz_j(1) + xhat)*c3ri) ) 
          dudy = dudy - preSw * ( &
               qxx_j*(cx2*yhatc4 - (2.0_dp*cx_j*ux_j(2) + yhat)*c3ri) + &
               qyy_j*(cy2*yhatc4 - (2.0_dp*cy_j*uy_j(2) + yhat)*c3ri) + &
               qzz_j*(cz2*yhatc4 - (2.0_dp*cz_j*uz_j(2) + yhat)*c3ri) ) 
          dudz = dudz - preSw * ( &
               qxx_j*(cx2*zhatc4 - (2.0_dp*cx_j*ux_j(3) + zhat)*c3ri) + &
               qyy_j*(cy2*zhatc4 - (2.0_dp*cy_j*uy_j(3) + zhat)*c3ri) + &
               qzz_j*(cz2*zhatc4 - (2.0_dp*cz_j*uz_j(3) + zhat)*c3ri) ) 
          
          dudux_j(1) = dudux_j(1) + preSw*(qxx_j*cx_j*xhatdot2)
          dudux_j(2) = dudux_j(2) + preSw*(qxx_j*cx_j*yhatdot2)
          dudux_j(3) = dudux_j(3) + preSw*(qxx_j*cx_j*zhatdot2)
          
          duduy_j(1) = duduy_j(1) + preSw*(qyy_j*cy_j*xhatdot2)
          duduy_j(2) = duduy_j(2) + preSw*(qyy_j*cy_j*yhatdot2)
          duduy_j(3) = duduy_j(3) + preSw*(qyy_j*cy_j*zhatdot2)
          
          duduz_j(1) = duduz_j(1) + preSw*(qzz_j*cz_j*xhatdot2)
          duduz_j(2) = duduz_j(2) + preSw*(qzz_j*cz_j*yhatdot2)
          duduz_j(3) = duduz_j(3) + preSw*(qzz_j*cz_j*zhatdot2)

           
       endif
    endif
    
    if (i_is_Dipole) then 

       if (j_is_Charge) then
          ! variables used by all the methods
          pref = electroMult * pre12 * q_j * mu_i                       
          preSw = sw*pref

          if (summationMethod .eq. REACTION_FIELD) then

             ri2 = riji * riji
             ri3 = ri2 * riji

             vterm = pref * ct_i * ( ri2 - preRF2*rij )
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             dudx = dudx + preSw * ( ri3*(uz_i(1) - 3.0_dp*ct_i*xhat) - &
                  preRF2*uz_i(1) )
             dudy = dudy + preSw * ( ri3*(uz_i(2) - 3.0_dp*ct_i*yhat) - &
                  preRF2*uz_i(2) )
             dudz = dudz + preSw * ( ri3*(uz_i(3) - 3.0_dp*ct_i*zhat) - &
                  preRF2*uz_i(3) )
             
             duduz_i(1) = duduz_i(1) + preSw * xhat * ( ri2 - preRF2*rij )
             duduz_i(2) = duduz_i(2) + preSw * yhat * ( ri2 - preRF2*rij )
             duduz_i(3) = duduz_i(3) + preSw * zhat * ( ri2 - preRF2*rij )

          else
             ! determine inverse r if we are using split dipoles
             if (i_is_SplitDipole) then
                BigR = sqrt(r2 + 0.25_dp * d_i * d_i)
                ri = 1.0_dp / BigR
                scale = rij * ri
             else
                ri = riji
                scale = 1.0_dp
             endif
 
             sc2 = scale * scale
              
             if (screeningMethod .eq. DAMPED) then
                ! assemble the damping variables
                call lookupUniformSpline1d(erfcSpline, rij, erfcVal, derfcVal)
                c1 = erfcVal*ri
                c2 = (-derfcVal + c1)*ri
                c3 = -2.0_dp*derfcVal*alpha2 + 3.0_dp*c2*ri
             else
                c1 = ri
                c2 = c1*ri
                c3 = 3.0_dp*c2*ri
             endif
            
             c2ri = c2*ri

             ! calculate the potential
             pot_term = c2 * scale
             vterm = pref * ct_i * pot_term
             vpair = vpair + vterm
             epot = epot + sw*vterm

             ! calculate derivatives for the forces and torques
             dudx = dudx + preSw * ( uz_i(1)*c2ri - ct_i*xhat*sc2*c3 )
             dudy = dudy + preSw * ( uz_i(2)*c2ri - ct_i*yhat*sc2*c3 )
             dudz = dudz + preSw * ( uz_i(3)*c2ri - ct_i*zhat*sc2*c3 )

             duduz_i(1) = duduz_i(1) + preSw * pot_term * xhat
             duduz_i(2) = duduz_i(2) + preSw * pot_term * yhat
             duduz_i(3) = duduz_i(3) + preSw * pot_term * zhat
             
          endif
       endif
       
       if (j_is_Dipole) then
          ! variables used by all methods
          ct_ij = uz_i(1)*uz_j(1) + uz_i(2)*uz_j(2) + uz_i(3)*uz_j(3)
          pref = electroMult * pre22 * mu_i * mu_j
          preSw = sw*pref

          if (summationMethod .eq. REACTION_FIELD) then
             ri2 = riji * riji
             ri3 = ri2 * riji
             ri4 = ri2 * ri2

             vterm = pref*( ri3*(ct_ij - 3.0_dp * ct_i * ct_j) - &
                  preRF2*ct_ij )
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             a1 = 5.0_dp * ct_i * ct_j - ct_ij
             
             dudx = dudx + preSw*3.0_dp*ri4*(a1*xhat-ct_i*uz_j(1)-ct_j*uz_i(1))
             dudy = dudy + preSw*3.0_dp*ri4*(a1*yhat-ct_i*uz_j(2)-ct_j*uz_i(2))
             dudz = dudz + preSw*3.0_dp*ri4*(a1*zhat-ct_i*uz_j(3)-ct_j*uz_i(3))
             
             duduz_i(1) = duduz_i(1) + preSw*(ri3*(uz_j(1)-3.0_dp*ct_j*xhat) &
                  - preRF2*uz_j(1))
             duduz_i(2) = duduz_i(2) + preSw*(ri3*(uz_j(2)-3.0_dp*ct_j*yhat) &
                  - preRF2*uz_j(2))
             duduz_i(3) = duduz_i(3) + preSw*(ri3*(uz_j(3)-3.0_dp*ct_j*zhat) &
                  - preRF2*uz_j(3))
             duduz_j(1) = duduz_j(1) + preSw*(ri3*(uz_i(1)-3.0_dp*ct_i*xhat) &
                  - preRF2*uz_i(1))
             duduz_j(2) = duduz_j(2) + preSw*(ri3*(uz_i(2)-3.0_dp*ct_i*yhat) &
                  - preRF2*uz_i(2))
             duduz_j(3) = duduz_j(3) + preSw*(ri3*(uz_i(3)-3.0_dp*ct_i*zhat) &
                  - preRF2*uz_i(3))

          else
             if (i_is_SplitDipole) then
                if (j_is_SplitDipole) then
                   BigR = sqrt(r2 + 0.25_dp * d_i * d_i + 0.25_dp * d_j * d_j)
                else
                   BigR = sqrt(r2 + 0.25_dp * d_i * d_i)
                endif
                ri = 1.0_dp / BigR
                scale = rij * ri                
             else
                if (j_is_SplitDipole) then
                   BigR = sqrt(r2 + 0.25_dp * d_j * d_j)
                   ri = 1.0_dp / BigR
                   scale = rij * ri                             
                else                
                   ri = riji
                   scale = 1.0_dp
                endif
             endif

             if (screeningMethod .eq. DAMPED) then
                ! assemble the damping variables
                call lookupUniformSpline1d(erfcSpline, rij, erfcVal, derfcVal)
                c1 = erfcVal*ri
                c2 = (-derfcVal + c1)*ri
                c3 = -2.0_dp*derfcVal*alpha2 + 3.0_dp*c2*ri
                c4 = -4.0_dp*derfcVal*alpha4 + 5.0_dp*c3*ri*ri
             else
                c1 = ri
                c2 = c1*ri
                c3 = 3.0_dp*c2*ri
                c4 = 5.0_dp*c3*ri*ri
             endif

             ! precompute variables for convenience
             sc2 = scale * scale
             cti3 = ct_i*sc2*c3
             ctj3 = ct_j*sc2*c3
             ctidotj = ct_i * ct_j * sc2        
             preSwSc = preSw*scale
             c2ri = c2*ri
             c3ri = c3*ri
             c4rij = c4*rij


             ! calculate the potential 
             pot_term = (ct_ij*c2ri - ctidotj*c3)
             vterm = pref * pot_term
             vpair = vpair + vterm
             epot = epot + sw*vterm

             ! calculate derivatives for the forces and torques
             dudx = dudx + preSwSc * ( ctidotj*xhat*c4rij - &
                  (ct_i*uz_j(1) + ct_j*uz_i(1) + ct_ij*xhat)*c3ri )
             dudy = dudy + preSwSc * ( ctidotj*yhat*c4rij - &
                  (ct_i*uz_j(2) + ct_j*uz_i(2) + ct_ij*yhat)*c3ri )
             dudz = dudz + preSwSc * ( ctidotj*zhat*c4rij - &
                  (ct_i*uz_j(3) + ct_j*uz_i(3) + ct_ij*zhat)*c3ri )

             duduz_i(1) = duduz_i(1) + preSw * ( uz_j(1)*c2ri - ctj3*xhat )
             duduz_i(2) = duduz_i(2) + preSw * ( uz_j(2)*c2ri - ctj3*yhat )
             duduz_i(3) = duduz_i(3) + preSw * ( uz_j(3)*c2ri - ctj3*zhat )
             
             duduz_j(1) = duduz_j(1) + preSw * ( uz_i(1)*c2ri - cti3*xhat )
             duduz_j(2) = duduz_j(2) + preSw * ( uz_i(2)*c2ri - cti3*yhat )
             duduz_j(3) = duduz_j(3) + preSw * ( uz_i(3)*c2ri - cti3*zhat )

          endif
       endif
    endif

    if (i_is_Quadrupole) then
       if (j_is_Charge) then
          ! precompute some necessary variables
          cx2 = cx_i * cx_i
          cy2 = cy_i * cy_i
          cz2 = cz_i * cz_i
          pref = electroMult * pre14 * q_j * one_third

          if (screeningMethod .eq. DAMPED) then
             ! assemble the damping variables
             call lookupUniformSpline1d(erfcSpline, rij, erfcVal, derfcVal)
             c1 = erfcVal*riji
             c2 = (-derfcVal + c1)*riji
             c3 = -2.0_dp*derfcVal*alpha2 + 3.0_dp*c2*riji
             c4 = -4.0_dp*derfcVal*alpha4 + 5.0_dp*c3*riji*riji
          else
             c1 = riji
             c2 = c1*riji
             c3 = 3.0_dp*c2*riji
             c4 = 5.0_dp*c3*riji*riji
          endif
          
          ! precompute some variables for convenience
          preSw = sw*pref
          c2ri = c2*riji
          c3ri = c3*riji
          c4rij = c4*rij
          xhatdot2 = 2.0_dp*xhat*c3
          yhatdot2 = 2.0_dp*yhat*c3
          zhatdot2 = 2.0_dp*zhat*c3
          xhatc4 = xhat*c4rij
          yhatc4 = yhat*c4rij
          zhatc4 = zhat*c4rij

          ! calculate the potential
          pot_term = ( qxx_i * (cx2*c3 - c2ri) + qyy_i * (cy2*c3 - c2ri) + &
               qzz_i * (cz2*c3 - c2ri) )

          vterm = pref * pot_term
          vpair = vpair + vterm
          epot = epot + sw*vterm

          ! calculate the derivatives for the forces and torques
          dudx = dudx - preSw * ( &
               qxx_i*(cx2*xhatc4 - (2.0_dp*cx_i*ux_i(1) + xhat)*c3ri) + &
               qyy_i*(cy2*xhatc4 - (2.0_dp*cy_i*uy_i(1) + xhat)*c3ri) + &
               qzz_i*(cz2*xhatc4 - (2.0_dp*cz_i*uz_i(1) + xhat)*c3ri) ) 
          dudy = dudy - preSw * ( &
               qxx_i*(cx2*yhatc4 - (2.0_dp*cx_i*ux_i(2) + yhat)*c3ri) + &
               qyy_i*(cy2*yhatc4 - (2.0_dp*cy_i*uy_i(2) + yhat)*c3ri) + &
               qzz_i*(cz2*yhatc4 - (2.0_dp*cz_i*uz_i(2) + yhat)*c3ri) ) 
          dudz = dudz - preSw * ( &
               qxx_i*(cx2*zhatc4 - (2.0_dp*cx_i*ux_i(3) + zhat)*c3ri) + &
               qyy_i*(cy2*zhatc4 - (2.0_dp*cy_i*uy_i(3) + zhat)*c3ri) + &
               qzz_i*(cz2*zhatc4 - (2.0_dp*cz_i*uz_i(3) + zhat)*c3ri) ) 
          
          dudux_i(1) = dudux_i(1) + preSw*(qxx_i*cx_i*xhatdot2)
          dudux_i(2) = dudux_i(2) + preSw*(qxx_i*cx_i*yhatdot2)
          dudux_i(3) = dudux_i(3) + preSw*(qxx_i*cx_i*zhatdot2)
          
          duduy_i(1) = duduy_i(1) + preSw*(qyy_i*cy_i*xhatdot2)
          duduy_i(2) = duduy_i(2) + preSw*(qyy_i*cy_i*yhatdot2)
          duduy_i(3) = duduy_i(3) + preSw*(qyy_i*cy_i*zhatdot2)
          
          duduz_i(1) = duduz_i(1) + preSw*(qzz_i*cz_i*xhatdot2)
          duduz_i(2) = duduz_i(2) + preSw*(qzz_i*cz_i*yhatdot2)
          duduz_i(3) = duduz_i(3) + preSw*(qzz_i*cz_i*zhatdot2)
       endif
    endif

    pot = pot + epot

    f1(1) = f1(1) + dudx
    f1(2) = f1(2) + dudy
    f1(3) = f1(3) + dudz

    if (i_is_Dipole .or. i_is_Quadrupole) then
       t1(1) = t1(1) - uz_i(2)*duduz_i(3) + uz_i(3)*duduz_i(2)
       t1(2) = t1(2) - uz_i(3)*duduz_i(1) + uz_i(1)*duduz_i(3)
       t1(3) = t1(3) - uz_i(1)*duduz_i(2) + uz_i(2)*duduz_i(1)
    endif
    if (i_is_Quadrupole) then
       t1(1) = t1(1) - ux_i(2)*dudux_i(3) + ux_i(3)*dudux_i(2)
       t1(2) = t1(2) - ux_i(3)*dudux_i(1) + ux_i(1)*dudux_i(3)
       t1(3) = t1(3) - ux_i(1)*dudux_i(2) + ux_i(2)*dudux_i(1)

       t1(1) = t1(1) - uy_i(2)*duduy_i(3) + uy_i(3)*duduy_i(2)
       t1(2) = t1(2) - uy_i(3)*duduy_i(1) + uy_i(1)*duduy_i(3)
       t1(3) = t1(3) - uy_i(1)*duduy_i(2) + uy_i(2)*duduy_i(1)
    endif

    if (j_is_Dipole .or. j_is_Quadrupole) then
       t2(1) = t2(1) - uz_j(2)*duduz_j(3) + uz_j(3)*duduz_j(2)
       t2(2) = t2(2) - uz_j(3)*duduz_j(1) + uz_j(1)*duduz_j(3)
       t2(3) = t2(3) - uz_j(1)*duduz_j(2) + uz_j(2)*duduz_j(1)
    endif
    if (j_is_Quadrupole) then
       t2(1) = t2(1) - ux_j(2)*dudux_j(3) + ux_j(3)*dudux_j(2)
       t2(2) = t2(2) - ux_j(3)*dudux_j(1) + ux_j(1)*dudux_j(3)
       t2(3) = t2(3) - ux_j(1)*dudux_j(2) + ux_j(2)*dudux_j(1)

       t2(1) = t2(1) - uy_j(2)*duduy_j(3) + uy_j(3)*duduy_j(2)
       t2(2) = t2(2) - uy_j(3)*duduy_j(1) + uy_j(1)*duduy_j(3)
       t2(3) = t2(3) - uy_j(1)*duduy_j(2) + uy_j(2)*duduy_j(1)
    endif

    return
  end subroutine doElectrostaticPair

  subroutine destroyElectrostaticTypes()

    if(allocated(ElectrostaticMap)) deallocate(ElectrostaticMap)

  end subroutine destroyElectrostaticTypes
       
  subroutine self_self(atom1, eFrame, skch, mypot, t, do_pot)
    logical, intent(in) :: do_pot
    integer, intent(in) :: atom1
    integer :: atid1
    real(kind=dp), dimension(9,nLocal) :: eFrame
    real(kind=dp), dimension(3,nLocal) :: t
    real(kind=dp) :: mu1, chg1, c1e
    real(kind=dp) :: preVal, epot, mypot, skch
    real(kind=dp) :: eix, eiy, eiz, self

    ! this is a local only array, so we use the local atom type id's:
    atid1 = atid(atom1)

    if (.not.summationMethodChecked) then
       call checkSummationMethod()
    endif
    
    if (summationMethod .eq. REACTION_FIELD) then
       if (ElectrostaticMap(atid1)%is_Dipole) then
          mu1 = getDipoleMoment(atid1)
          
          preVal = pre22 * preRF2 * mu1*mu1
          mypot = mypot - 0.5_dp*preVal
          
          ! The self-correction term adds into the reaction field vector
          
          eix = preVal * eFrame(3,atom1)
          eiy = preVal * eFrame(6,atom1)
          eiz = preVal * eFrame(9,atom1)
          
          ! once again, this is self-self, so only the local arrays are needed
          ! even for MPI jobs:
          
          t(1,atom1)=t(1,atom1) - eFrame(6,atom1)*eiz + &
               eFrame(9,atom1)*eiy
          t(2,atom1)=t(2,atom1) - eFrame(9,atom1)*eix + &
               eFrame(3,atom1)*eiz
          t(3,atom1)=t(3,atom1) - eFrame(3,atom1)*eiy + &
               eFrame(6,atom1)*eix
          
       endif

    elseif ( (summationMethod .eq. SHIFTED_FORCE) .or. &
         (summationMethod .eq. SHIFTED_POTENTIAL) ) then
       if (ElectrostaticMap(atid1)%is_Charge) then
          chg1 = getCharge(atid1)         
          if (screeningMethod .eq. DAMPED) then
             self = - 0.5_dp * (c1c+alphaPi) * chg1 * (chg1+skch) * pre11
          else             
             self = - 0.5_dp * rcuti * chg1 * (chg1+skch) * pre11
          endif

          mypot = mypot + self
       endif
    endif
    
    

    return
  end subroutine self_self

  subroutine rf_self_excludes(atom1, atom2, sw, electroMult, eFrame, d, &
       rij, vpair, myPot, f, t, do_pot)
    logical, intent(in) :: do_pot
    integer, intent(in) :: atom1
    integer, intent(in) :: atom2
    logical :: i_is_Charge, j_is_Charge
    logical :: i_is_Dipole, j_is_Dipole
    integer :: atid1
    integer :: atid2
    real(kind=dp), intent(in) :: rij
    real(kind=dp), intent(in) :: sw, electroMult
    real(kind=dp), intent(in), dimension(3) :: d
    real(kind=dp), intent(inout) :: vpair
    real(kind=dp), dimension(9,nLocal) :: eFrame
    real(kind=dp), dimension(3,nLocal) :: f
    real(kind=dp), dimension(3,nLocal) :: t
    real (kind = dp), dimension(3) :: duduz_i
    real (kind = dp), dimension(3) :: duduz_j
    real (kind = dp), dimension(3) :: uz_i
    real (kind = dp), dimension(3) :: uz_j
    real(kind=dp) :: q_i, q_j, mu_i, mu_j
    real(kind=dp) :: xhat, yhat, zhat
    real(kind=dp) :: ct_i, ct_j
    real(kind=dp) :: ri2, ri3, riji, vterm
    real(kind=dp) :: pref, preVal, rfVal, myPot
    real(kind=dp) :: dudx, dudy, dudz, dudr

    if (.not.summationMethodChecked) then
       call checkSummationMethod()
    endif

    dudx = zero
    dudy = zero
    dudz = zero

    riji = 1.0_dp/rij

    xhat = d(1) * riji
    yhat = d(2) * riji
    zhat = d(3) * riji

    ! this is a local only array, so we use the local atom type id's:
    atid1 = atid(atom1)
    atid2 = atid(atom2)
    i_is_Charge = ElectrostaticMap(atid1)%is_Charge
    j_is_Charge = ElectrostaticMap(atid2)%is_Charge
    i_is_Dipole = ElectrostaticMap(atid1)%is_Dipole
    j_is_Dipole = ElectrostaticMap(atid2)%is_Dipole

    if (i_is_Charge.and.j_is_Charge) then
       q_i = ElectrostaticMap(atid1)%charge
       q_j = ElectrostaticMap(atid2)%charge
       
       preVal = electroMult * pre11 * q_i * q_j
       rfVal = preRF*rij*rij
       vterm = preVal * rfVal
       
       myPot = myPot + sw*vterm
       
       dudr  = sw*preVal * 2.0_dp*rfVal*riji
       
       dudx = dudx + dudr * xhat
       dudy = dudy + dudr * yhat
       dudz = dudz + dudr * zhat
       
    elseif (i_is_Charge.and.j_is_Dipole) then
       q_i = ElectrostaticMap(atid1)%charge
       mu_j = ElectrostaticMap(atid2)%dipole_moment
       uz_j(1) = eFrame(3,atom2)
       uz_j(2) = eFrame(6,atom2)
       uz_j(3) = eFrame(9,atom2)
       ct_j = uz_j(1)*xhat + uz_j(2)*yhat + uz_j(3)*zhat
       
       ri2 = riji * riji
       ri3 = ri2 * riji
       
       pref = electroMult * pre12 * q_i * mu_j
       vterm = - pref * ct_j * ( ri2 - preRF2*rij )
       myPot = myPot + sw*vterm
       
       dudx = dudx - sw*pref*( ri3*(uz_j(1)-3.0_dp*ct_j*xhat) &
            - preRF2*uz_j(1) )
       dudy = dudy - sw*pref*( ri3*(uz_j(2)-3.0_dp*ct_j*yhat) &
            - preRF2*uz_j(2) )
       dudz = dudz - sw*pref*( ri3*(uz_j(3)-3.0_dp*ct_j*zhat) &
            - preRF2*uz_j(3) )
       
       duduz_j(1) = duduz_j(1) - sw * pref * xhat * ( ri2 - preRF2*rij )
       duduz_j(2) = duduz_j(2) - sw * pref * yhat * ( ri2 - preRF2*rij )
       duduz_j(3) = duduz_j(3) - sw * pref * zhat * ( ri2 - preRF2*rij )
       
    elseif (i_is_Dipole.and.j_is_Charge) then
       mu_i = ElectrostaticMap(atid1)%dipole_moment
       q_j = ElectrostaticMap(atid2)%charge
       uz_i(1) = eFrame(3,atom1)
       uz_i(2) = eFrame(6,atom1)
       uz_i(3) = eFrame(9,atom1)
       ct_i = uz_i(1)*xhat + uz_i(2)*yhat + uz_i(3)*zhat
       
       ri2 = riji * riji
       ri3 = ri2 * riji
       
       pref = electroMult * pre12 * q_j * mu_i
       vterm = pref * ct_i * ( ri2 - preRF2*rij )
       myPot = myPot + sw*vterm
       
       dudx = dudx + sw*pref*( ri3*(uz_i(1)-3.0_dp*ct_i*xhat) &
            - preRF2*uz_i(1) )
       dudy = dudy + sw*pref*( ri3*(uz_i(2)-3.0_dp*ct_i*yhat) &
            - preRF2*uz_i(2) )
       dudz = dudz + sw*pref*( ri3*(uz_i(3)-3.0_dp*ct_i*zhat) &
            - preRF2*uz_i(3) )
       
       duduz_i(1) = duduz_i(1) + sw * pref * xhat * ( ri2 - preRF2*rij )
       duduz_i(2) = duduz_i(2) + sw * pref * yhat * ( ri2 - preRF2*rij )
       duduz_i(3) = duduz_i(3) + sw * pref * zhat * ( ri2 - preRF2*rij )
       
    endif
       

    ! accumulate the forces and torques resulting from the self term
    f(1,atom1) = f(1,atom1) + dudx
    f(2,atom1) = f(2,atom1) + dudy
    f(3,atom1) = f(3,atom1) + dudz
    
    f(1,atom2) = f(1,atom2) - dudx
    f(2,atom2) = f(2,atom2) - dudy
    f(3,atom2) = f(3,atom2) - dudz
    
    if (i_is_Dipole) then
       t(1,atom1)=t(1,atom1) - uz_i(2)*duduz_i(3) + uz_i(3)*duduz_i(2)
       t(2,atom1)=t(2,atom1) - uz_i(3)*duduz_i(1) + uz_i(1)*duduz_i(3)
       t(3,atom1)=t(3,atom1) - uz_i(1)*duduz_i(2) + uz_i(2)*duduz_i(1)
    elseif (j_is_Dipole) then
       t(1,atom2)=t(1,atom2) - uz_j(2)*duduz_j(3) + uz_j(3)*duduz_j(2)
       t(2,atom2)=t(2,atom2) - uz_j(3)*duduz_j(1) + uz_j(1)*duduz_j(3)
       t(3,atom2)=t(3,atom2) - uz_j(1)*duduz_j(2) + uz_j(2)*duduz_j(1)
    endif

    return
  end subroutine rf_self_excludes

  subroutine accumulate_box_dipole(atom1, eFrame, d, pChg, nChg, pChgPos, &
       nChgPos, dipVec, pChgCount, nChgCount)
    integer, intent(in) :: atom1
    logical :: i_is_Charge
    logical :: i_is_Dipole
    integer :: atid1
    integer :: pChgCount
    integer :: nChgCount
    real(kind=dp), intent(in), dimension(3) :: d
    real(kind=dp), dimension(9,nLocal) :: eFrame
    real(kind=dp) :: pChg
    real(kind=dp) :: nChg
    real(kind=dp), dimension(3) :: pChgPos
    real(kind=dp), dimension(3) :: nChgPos
    real(kind=dp), dimension(3) :: dipVec
    real(kind=dp), dimension(3) :: uz_i
    real(kind=dp), dimension(3) :: pos
    real(kind=dp) :: q_i, mu_i
    real(kind=dp) :: pref, preVal

    if (.not.summationMethodChecked) then
       call checkSummationMethod()
    endif

    ! this is a local only array, so we use the local atom type id's:
    atid1 = atid(atom1)
    i_is_Charge = ElectrostaticMap(atid1)%is_Charge
    i_is_Dipole = ElectrostaticMap(atid1)%is_Dipole
    
    if (i_is_Charge) then
       q_i = ElectrostaticMap(atid1)%charge
       ! convert to the proper units
       q_i = q_i * chargeToC
       pos = d * angstromToM

       if (q_i.le.0.0_dp) then
          nChg = nChg - q_i
          nChgPos(1) = nChgPos(1) + pos(1)
          nChgPos(2) = nChgPos(2) + pos(2)
          nChgPos(3) = nChgPos(3) + pos(3)
          nChgCount = nChgCount + 1

       else
          pChg = pChg + q_i
          pChgPos(1) = pChgPos(1) + pos(1)
          pChgPos(2) = pChgPos(2) + pos(2)
          pChgPos(3) = pChgPos(3) + pos(3)
          pChgCount = pChgCount + 1

       endif

    endif
    
    if (i_is_Dipole) then
       mu_i = ElectrostaticMap(atid1)%dipole_moment
       uz_i(1) = eFrame(3,atom1)
       uz_i(2) = eFrame(6,atom1)
       uz_i(3) = eFrame(9,atom1)
       ! convert to the proper units
       mu_i = mu_i * debyeToCm

       dipVec(1) = dipVec(1) + uz_i(1)*mu_i
       dipVec(2) = dipVec(2) + uz_i(2)*mu_i
       dipVec(3) = dipVec(3) + uz_i(3)*mu_i

    endif
   
    return
  end subroutine accumulate_box_dipole

end module electrostatic_module
