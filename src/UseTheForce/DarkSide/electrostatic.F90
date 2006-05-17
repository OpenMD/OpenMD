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
#ifdef IS_MPI
  use mpiSimulation
#endif
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
  real(kind=dp), save :: f0 = 1.0_DP
  real(kind=dp), save :: f1 = 1.0_DP
  real(kind=dp), save :: f2 = 0.0_DP
  real(kind=dp), save :: f3 = 0.0_DP
  real(kind=dp), save :: f4 = 0.0_DP
  real(kind=dp), save :: f0c = 1.0_DP
  real(kind=dp), save :: f1c = 1.0_DP
  real(kind=dp), save :: f2c = 0.0_DP
  real(kind=dp), save :: f3c = 0.0_DP
  real(kind=dp), save :: f4c = 0.0_DP
  real(kind=dp), save :: df0 = 0.0_DP
  type(cubicSpline), save :: f0spline
  logical, save :: haveElectroSpline = .false.
  real(kind=dp), save :: one_third = 1.0_DP / 3.0_DP

#if defined(__IFC) || defined(__PGI)
! error function for ifc version > 7.
  double precision, external :: derfc
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
       yvals(i) = derfc(dampingAlpha*rval)
    enddo

    call newSpline(f0spline, xvals, yvals, .true.)

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
       f0c = derfc(dampingAlpha*defaultCutoff)
       f1c = alphaPi*defaultCutoff*constEXP + f0c
       f2c = alphaPi*2.0_dp*alpha2*constEXP
       f3c = alphaPi*2.0_dp*alpha2*constEXP*defaultCutoff2*defaultCutoff
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


  subroutine doElectrostaticPair(atom1, atom2, d, rij, r2, rcut, sw, &
       vpair, fpair, pot, eFrame, f, t, do_pot)

    logical, intent(in) :: do_pot

    integer, intent(in) :: atom1, atom2
    integer :: localError

    real(kind=dp), intent(in) :: rij, r2, sw, rcut
    real(kind=dp), intent(in), dimension(3) :: d
    real(kind=dp), intent(inout) :: vpair
    real(kind=dp), intent(inout), dimension(3) :: fpair    

    real( kind = dp ) :: pot
    real( kind = dp ), dimension(9,nLocal) :: eFrame
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), dimension(3,nLocal) :: felec
    real( kind = dp ), dimension(3,nLocal) :: t

    real (kind = dp), dimension(3) :: ux_i, uy_i, uz_i
    real (kind = dp), dimension(3) :: ux_j, uy_j, uz_j
    real (kind = dp), dimension(3) :: dudux_i, duduy_i, duduz_i
    real (kind = dp), dimension(3) :: dudux_j, duduy_j, duduz_j

    logical :: i_is_Charge, i_is_Dipole, i_is_SplitDipole, i_is_Quadrupole
    logical :: j_is_Charge, j_is_Dipole, j_is_SplitDipole, j_is_Quadrupole
    logical :: i_is_Tap, j_is_Tap
    integer :: me1, me2, id1, id2
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
    real (kind=dp) :: cti3, ctj3, ctidotj
    real (kind=dp) :: ri7damp, ri5damp, prei3, prei4
    real (kind=dp) :: xhatdot2, yhatdot2, zhatdot2
    real (kind=dp) :: xhatdot5, yhatdot5, zhatdot5

    if (.not.summationMethodChecked) then
       call checkSummationMethod()
    endif

#ifdef IS_MPI
    me1 = atid_Row(atom1)
    me2 = atid_Col(atom2)
#else
    me1 = atid(atom1)
    me2 = atid(atom2)
#endif

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
#ifdef IS_MPI
       uz_i(1) = eFrame_Row(3,atom1)
       uz_i(2) = eFrame_Row(6,atom1)
       uz_i(3) = eFrame_Row(9,atom1)
#else
       uz_i(1) = eFrame(3,atom1)
       uz_i(2) = eFrame(6,atom1)
       uz_i(3) = eFrame(9,atom1)
#endif
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
#ifdef IS_MPI
       ux_i(1) = eFrame_Row(1,atom1)
       ux_i(2) = eFrame_Row(4,atom1)
       ux_i(3) = eFrame_Row(7,atom1)
       uy_i(1) = eFrame_Row(2,atom1)
       uy_i(2) = eFrame_Row(5,atom1)
       uy_i(3) = eFrame_Row(8,atom1)
       uz_i(1) = eFrame_Row(3,atom1)
       uz_i(2) = eFrame_Row(6,atom1)
       uz_i(3) = eFrame_Row(9,atom1)
#else
       ux_i(1) = eFrame(1,atom1)
       ux_i(2) = eFrame(4,atom1)
       ux_i(3) = eFrame(7,atom1)
       uy_i(1) = eFrame(2,atom1)
       uy_i(2) = eFrame(5,atom1)
       uy_i(3) = eFrame(8,atom1)
       uz_i(1) = eFrame(3,atom1)
       uz_i(2) = eFrame(6,atom1)
       uz_i(3) = eFrame(9,atom1)
#endif
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
#ifdef IS_MPI
       uz_j(1) = eFrame_Col(3,atom2)
       uz_j(2) = eFrame_Col(6,atom2)
       uz_j(3) = eFrame_Col(9,atom2)
#else
       uz_j(1) = eFrame(3,atom2)
       uz_j(2) = eFrame(6,atom2)
       uz_j(3) = eFrame(9,atom2)
#endif
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
#ifdef IS_MPI
       ux_j(1) = eFrame_Col(1,atom2)
       ux_j(2) = eFrame_Col(4,atom2)
       ux_j(3) = eFrame_Col(7,atom2)
       uy_j(1) = eFrame_Col(2,atom2)
       uy_j(2) = eFrame_Col(5,atom2)
       uy_j(3) = eFrame_Col(8,atom2)
       uz_j(1) = eFrame_Col(3,atom2)
       uz_j(2) = eFrame_Col(6,atom2)
       uz_j(3) = eFrame_Col(9,atom2)
#else
       ux_j(1) = eFrame(1,atom2)
       ux_j(2) = eFrame(4,atom2)
       ux_j(3) = eFrame(7,atom2)
       uy_j(1) = eFrame(2,atom2)
       uy_j(2) = eFrame(5,atom2)
       uy_j(3) = eFrame(8,atom2)
       uz_j(1) = eFrame(3,atom2)
       uz_j(2) = eFrame(6,atom2)
       uz_j(3) = eFrame(9,atom2)
#endif
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
             call lookupUniformSpline1d(f0spline, rij, f0, df0)
             f1 = -rij * df0 + f0
          endif

          preVal = pre11 * q_i * q_j

          if (summationMethod .eq. SHIFTED_POTENTIAL) then
             vterm = preVal * (riji*f0 - rcuti*f0c)
             
             dudr  = -sw * preVal * riji * riji * f1
  
          elseif (summationMethod .eq. SHIFTED_FORCE) then
             vterm = preVal * ( riji*f0 - rcuti*f0c + &
                  f1c*rcuti2*(rij-defaultCutoff) )
             
             dudr  = -sw*preVal * (riji*riji*f1 - rcuti2*f1c)
  
          elseif (summationMethod .eq. REACTION_FIELD) then
             rfVal = preRF*rij*rij
             vterm = preVal * ( riji + rfVal )
             
             dudr  = sw * preVal * ( 2.0_dp*rfVal - riji )*riji
  
          else
             vterm = preVal * riji*f0
             
             dudr  = - sw * preVal * riji*riji*f1
  
          endif

          vpair = vpair + vterm
          epot = epot + sw*vterm

          dudx = dudx + dudr * xhat
          dudy = dudy + dudr * yhat
          dudz = dudz + dudr * zhat

       endif

       if (j_is_Dipole) then
          if (screeningMethod .eq. DAMPED) then
             ! assemble the damping variables
             call lookupUniformSpline1d(f0spline, rij, f0, df0)
             f1 = -rij * df0 + f0
             f3 = -2.0_dp*alpha2*df0*rij*rij*rij
          endif

          pref = pre12 * q_i * mu_j

          if (summationMethod .eq. REACTION_FIELD) then
             ri2 = riji * riji
             ri3 = ri2 * riji
    
             vterm = - pref * ct_j * ( ri2 - preRF2*rij )
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             dudx = dudx - sw*pref*( ri3*(uz_j(1) - 3.0_dp*ct_j*xhat) - &
                                     preRF2*uz_j(1) )
             dudy = dudy - sw*pref*( ri3*(uz_j(2) - 3.0_dp*ct_j*yhat) - &
                                     preRF2*uz_j(2) )
             dudz = dudz - sw*pref*( ri3*(uz_j(3) - 3.0_dp*ct_j*zhat) - &
                                     preRF2*uz_j(3) )         
             duduz_j(1) = duduz_j(1) - sw*pref * xhat * ( ri2 - preRF2*rij )
             duduz_j(2) = duduz_j(2) - sw*pref * yhat * ( ri2 - preRF2*rij )
             duduz_j(3) = duduz_j(3) - sw*pref * zhat * ( ri2 - preRF2*rij )

          else
             if (j_is_SplitDipole) then
                BigR = sqrt(r2 + 0.25_dp * d_j * d_j)
                ri = 1.0_dp / BigR
                scale = rij * ri
             else
                ri = riji
                scale = 1.0_dp
             endif
             
             ri2 = ri * ri
             ri3 = ri2 * ri
             sc2 = scale * scale

             pot_term =  ri2 * scale * f1
             vterm = -pref * ct_j * pot_term
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             prei3 = sw*pref*ri3
             ri5damp = 3.0_dp*f1 + f3
             
             dudx = dudx - prei3 * ( uz_j(1)*f1 - ct_j*xhat*sc2*ri5damp )
             dudy = dudy - prei3 * ( uz_j(2)*f1 - ct_j*yhat*sc2*ri5damp )
             dudz = dudz - prei3 * ( uz_j(3)*f1 - ct_j*zhat*sc2*ri5damp )
                          
             duduz_j(1) = duduz_j(1) - sw*pref * pot_term * xhat
             duduz_j(2) = duduz_j(2) - sw*pref * pot_term * yhat
             duduz_j(3) = duduz_j(3) - sw*pref * pot_term * zhat

          endif
       endif

       if (j_is_Quadrupole) then
          if (screeningMethod .eq. DAMPED) then
             ! assemble the damping variables
             call lookupUniformSpline1d(f0spline, rij, f0, df0)
             f1 = -rij * df0 + f0
             f2 = -2.0_dp*alpha2*df0
             f3 = f2*r2*rij
             f4 = 0.4_dp*alpha2*f3*r2
          endif
          ri5damp = f1 + f3*one_third
          ri7damp = ri5damp + f4

          ri2 = riji * riji
          ri3 = ri2 * riji
          cx2 = cx_j * cx_j
          cy2 = cy_j * cy_j
          cz2 = cz_j * cz_j

          pref =  pre14 * q_i * one_third

          pot_term = ri3*( qxx_j*(3.0_dp*cx2*ri5damp - f1) + &
               qyy_j*(3.0_dp*cy2*ri5damp - f1) + &
               qzz_j*(3.0_dp*cz2*ri5damp - f1) )
          vterm = pref * pot_term
          vpair = vpair + vterm
          epot = epot + sw*vterm

          ! precompute variables for convenience (and obfuscation unfortunatly)
          prei3 = 3.0_dp*sw*pref*ri3
          prei4 = prei3*riji
          xhatdot2 = xhat*2.0_dp * ri5damp
          yhatdot2 = yhat*2.0_dp * ri5damp
          zhatdot2 = zhat*2.0_dp * ri5damp
          xhatdot5 = xhat*5.0_dp * ri7damp
          yhatdot5 = yhat*5.0_dp * ri7damp
          zhatdot5 = zhat*5.0_dp * ri7damp

          dudx = dudx - prei4 * ( &
               qxx_j*(cx2*xhatdot5 - (2.0_dp*cx_j*ux_j(1) + xhat)*ri5damp) + &
               qyy_j*(cy2*xhatdot5 - (2.0_dp*cy_j*uy_j(1) + xhat)*ri5damp) + &
               qzz_j*(cz2*xhatdot5 - (2.0_dp*cz_j*uz_j(1) + xhat)*ri5damp) ) 
          dudy = dudy - prei4 * ( &
               qxx_j*(cx2*yhatdot5 - (2.0_dp*cx_j*ux_j(2) + yhat)*ri5damp) + &
               qyy_j*(cy2*yhatdot5 - (2.0_dp*cy_j*uy_j(2) + yhat)*ri5damp) + &
               qzz_j*(cz2*yhatdot5 - (2.0_dp*cz_j*uz_j(2) + yhat)*ri5damp) ) 
          dudz = dudz - prei4 * ( &
               qxx_j*(cx2*zhatdot5 - (2.0_dp*cx_j*ux_j(3) + zhat)*ri5damp) + &
               qyy_j*(cy2*zhatdot5 - (2.0_dp*cy_j*uy_j(3) + zhat)*ri5damp) + &
               qzz_j*(cz2*zhatdot5 - (2.0_dp*cz_j*uz_j(3) + zhat)*ri5damp) ) 
          
          dudux_j(1) = dudux_j(1) + prei3*(qxx_j*cx_j*xhatdot2)
          dudux_j(2) = dudux_j(2) + prei3*(qxx_j*cx_j*yhatdot2)
          dudux_j(3) = dudux_j(3) + prei3*(qxx_j*cx_j*zhatdot2)
          
          duduy_j(1) = duduy_j(1) + prei3*(qyy_j*cy_j*xhatdot2)
          duduy_j(2) = duduy_j(2) + prei3*(qyy_j*cy_j*yhatdot2)
          duduy_j(3) = duduy_j(3) + prei3*(qyy_j*cy_j*zhatdot2)
          
          duduz_j(1) = duduz_j(1) + prei3*(qzz_j*cz_j*xhatdot2)
          duduz_j(2) = duduz_j(2) + prei3*(qzz_j*cz_j*yhatdot2)
          duduz_j(3) = duduz_j(3) + prei3*(qzz_j*cz_j*zhatdot2)

           
       endif
    endif
    
    if (i_is_Dipole) then 

       if (j_is_Charge) then
          if (screeningMethod .eq. DAMPED) then
             ! assemble the damping variables
             call lookupUniformSpline1d(f0spline, rij, f0, df0)
             f1 = -rij * df0 + f0
             f3 = -2.0_dp*alpha2*df0*r2*rij
          endif
          
          pref = pre12 * q_j * mu_i
          
          if (summationMethod .eq. REACTION_FIELD) then

             ri2 = riji * riji
             ri3 = ri2 * riji

             vterm = pref * ct_i * ( ri2 - preRF2*rij )
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             dudx = dudx + sw*pref * ( ri3*(uz_i(1) - 3.0_dp*ct_i*xhat) - &
                  preRF2*uz_i(1) )
             dudy = dudy + sw*pref * ( ri3*(uz_i(2) - 3.0_dp*ct_i*yhat) - &
                  preRF2*uz_i(2) )
             dudz = dudz + sw*pref * ( ri3*(uz_i(3) - 3.0_dp*ct_i*zhat) - &
                  preRF2*uz_i(3) )
             
             duduz_i(1) = duduz_i(1) + sw*pref * xhat * ( ri2 - preRF2*rij )
             duduz_i(2) = duduz_i(2) + sw*pref * yhat * ( ri2 - preRF2*rij )
             duduz_i(3) = duduz_i(3) + sw*pref * zhat * ( ri2 - preRF2*rij )

          else
             if (i_is_SplitDipole) then
                BigR = sqrt(r2 + 0.25_dp * d_i * d_i)
                ri = 1.0_dp / BigR
                scale = rij * ri
             else
                ri = riji
                scale = 1.0_dp
             endif
             
             ri2 = ri * ri
             ri3 = ri2 * ri
             sc2 = scale * scale

             pot_term = ri2 * f1 * scale
             vterm = pref * ct_i * pot_term
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             prei3 = sw*pref*ri3
             ri5damp = 3.0_dp*f1 + f3
             
             dudx = dudx + prei3 * ( uz_i(1)*f1 - ct_i*xhat*sc2*ri5damp )
             dudy = dudy + prei3 * ( uz_i(2)*f1 - ct_i*yhat*sc2*ri5damp )
             dudz = dudz + prei3 * ( uz_i(3)*f1 - ct_i*zhat*sc2*ri5damp )

             duduz_i(1) = duduz_i(1) + sw*pref * pot_term * xhat
             duduz_i(2) = duduz_i(2) + sw*pref * pot_term * yhat
             duduz_i(3) = duduz_i(3) + sw*pref * pot_term * zhat
             
          endif
       endif
       
       if (j_is_Dipole) then
!!$          if (screeningMethod .eq. DAMPED) then
!!$             ! assemble the damping variables
!!$             call lookupUniformSpline1d(f0spline, rij, f0, df0)
!!$             f1 = -rij * df0 + f0
!!$             f2 = -2.0_dp*alpha2*df0
!!$             f3 = f2*r2*rij
!!$             f4 = 0.4_dp*alpha2*f3*r2
!!$          endif

          ct_ij = uz_i(1)*uz_j(1) + uz_i(2)*uz_j(2) + uz_i(3)*uz_j(3)
          
          ri2 = riji * riji
          ri3 = ri2 * riji
          ri4 = ri2 * ri2
          
          pref = pre22 * mu_i * mu_j

          if (summationMethod .eq. REACTION_FIELD) then
             vterm = pref*( ri3*(ct_ij - 3.0_dp * ct_i * ct_j) - &
                  preRF2*ct_ij )
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             a1 = 5.0_dp * ct_i * ct_j - ct_ij
             
             dudx = dudx + sw*pref*3.0_dp*ri4 &
                             * (a1*xhat-ct_i*uz_j(1)-ct_j*uz_i(1))
             dudy = dudy + sw*pref*3.0_dp*ri4 &
                             * (a1*yhat-ct_i*uz_j(2)-ct_j*uz_i(2))
             dudz = dudz + sw*pref*3.0_dp*ri4 &
                             * (a1*zhat-ct_i*uz_j(3)-ct_j*uz_i(3))
             
             duduz_i(1) = duduz_i(1) + sw*pref*(ri3*(uz_j(1)-3.0_dp*ct_j*xhat) &
                  - preRF2*uz_j(1))
             duduz_i(2) = duduz_i(2) + sw*pref*(ri3*(uz_j(2)-3.0_dp*ct_j*yhat) &
                  - preRF2*uz_j(2))
             duduz_i(3) = duduz_i(3) + sw*pref*(ri3*(uz_j(3)-3.0_dp*ct_j*zhat) &
                  - preRF2*uz_j(3))
             duduz_j(1) = duduz_j(1) + sw*pref*(ri3*(uz_i(1)-3.0_dp*ct_i*xhat) &
                  - preRF2*uz_i(1))
             duduz_j(2) = duduz_j(2) + sw*pref*(ri3*(uz_i(2)-3.0_dp*ct_i*yhat) &
                  - preRF2*uz_i(2))
             duduz_j(3) = duduz_j(3) + sw*pref*(ri3*(uz_i(3)-3.0_dp*ct_i*zhat) &
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
             
             sc2 = scale * scale

             pot_term = (ct_ij - 3.0_dp * ct_i * ct_j * sc2)
!!$             vterm = pref * ( ri3*pot_term*f1 + (ct_i * ct_j * sc2)*f2 )
             vterm = pref * ri3 * pot_term
             vpair = vpair + vterm
             epot = epot + sw*vterm
             
             ! precompute variables for convenience (and obfuscation 
             ! unfortunatly)
!!$             ri5damp = f1 + f3*one_third
!!$             ri7damp = 5.0_dp*(ri5damp + f4)
             prei3 = sw*pref*ri3
             prei4 = 3.0_dp*sw*pref*ri4*scale
!!$             cti3 = 3.0_dp*ct_i*sc2*ri5damp
!!$             ctj3 = 3.0_dp*ct_j*sc2*ri5damp
             cti3 = 3.0_dp*ct_i*sc2
             ctj3 = 3.0_dp*ct_j*sc2
             ctidotj = ct_i * ct_j * sc2

             dudx = dudx + prei4 * ( 5.0_dp*ctidotj*xhat - &
                  (ct_i*uz_j(1) + ct_j*uz_i(1) + ct_ij*xhat) )
             dudy = dudy + prei4 * ( 5.0_dp*ctidotj*yhat - &
                  (ct_i*uz_j(2) + ct_j*uz_i(2) + ct_ij*yhat) )
             dudz = dudz + prei4 * ( 5.0_dp*ctidotj*zhat - &
                  (ct_i*uz_j(3) + ct_j*uz_i(3) + ct_ij*zhat) )

             duduz_i(1) = duduz_i(1) + prei3 * ( uz_j(1) - ctj3*xhat )
             duduz_i(2) = duduz_i(2) + prei3 * ( uz_j(2) - ctj3*yhat )
             duduz_i(3) = duduz_i(3) + prei3 * ( uz_j(3) - ctj3*zhat )
             
             duduz_j(1) = duduz_j(1) + prei3 * ( uz_i(1) - cti3*xhat )
             duduz_j(2) = duduz_j(2) + prei3 * ( uz_i(2) - cti3*yhat )
             duduz_j(3) = duduz_j(3) + prei3 * ( uz_i(3) - cti3*zhat )

!!$             dudx = dudx + prei4 * ( ctidotj*xhat*ri7damp - &
!!$                  (ct_i*uz_j(1) + ct_j*uz_i(1) + ct_ij*xhat)*ri5damp )
!!$             dudy = dudy + prei4 * ( ctidotj*yhat*ri7damp - &
!!$                  (ct_i*uz_j(2) + ct_j*uz_i(2) + ct_ij*yhat)*ri5damp )
!!$             dudz = dudz + prei4 * ( ctidotj*zhat*ri7damp - &
!!$                  (ct_i*uz_j(3) + ct_j*uz_i(3) + ct_ij*zhat)*ri5damp )
!!$
!!$             duduz_i(1) = duduz_i(1) + prei3 * ( uz_j(1)*f1 - ctj3*xhat )
!!$             duduz_i(2) = duduz_i(2) + prei3 * ( uz_j(2)*f1 - ctj3*yhat )
!!$             duduz_i(3) = duduz_i(3) + prei3 * ( uz_j(3)*f1 - ctj3*zhat )
!!$             
!!$             duduz_j(1) = duduz_j(1) + prei3 * ( uz_i(1)*f1 - cti3*xhat )
!!$             duduz_j(2) = duduz_j(2) + prei3 * ( uz_i(2)*f1 - cti3*yhat )
!!$             duduz_j(3) = duduz_j(3) + prei3 * ( uz_i(3)*f1 - cti3*zhat )

          endif
       endif
    endif

    if (i_is_Quadrupole) then
       if (j_is_Charge) then
          if (screeningMethod .eq. DAMPED) then
             ! assemble the damping variables
             call lookupUniformSpline1d(f0spline, rij, f0, df0)
             f1 = -rij * df0 + f0
             f2 = -2.0_dp*alpha2*df0
             f3 = f2*r2*rij
             f4 = 0.4_dp*alpha2*f3*r2
          endif
          ri5damp = f1 + f3*one_third
          ri7damp = ri5damp + f4

          ri2 = riji * riji
          ri3 = ri2 * riji
          ri4 = ri2 * ri2
          cx2 = cx_i * cx_i
          cy2 = cy_i * cy_i
          cz2 = cz_i * cz_i

          pref = pre14 * q_j * one_third

          pot_term = ri3 * ( qxx_i * (3.0_dp*cx2*ri5damp - f1) + &
                             qyy_i * (3.0_dp*cy2*ri5damp - f1) + &
                             qzz_i * (3.0_dp*cz2*ri5damp - f1) )

          vterm = pref * pot_term
          vpair = vpair + vterm
          epot = epot + sw*vterm
 
          ! precompute variables for convenience (and obfuscation unfortunatly)
          prei3 = 3.0_dp*sw*pref*ri3
          prei4 = prei3*riji
          xhatdot2 = xhat*2.0_dp * ri5damp
          yhatdot2 = yhat*2.0_dp * ri5damp
          zhatdot2 = zhat*2.0_dp * ri5damp
          xhatdot5 = xhat*5.0_dp * ri7damp
          yhatdot5 = yhat*5.0_dp * ri7damp
          zhatdot5 = zhat*5.0_dp * ri7damp

          dudx = dudx - prei4 * ( &
               qxx_i*(cx2*xhatdot5 - (2.0_dp*cx_i*ux_i(1) + xhat)*ri5damp) + &
               qyy_i*(cy2*xhatdot5 - (2.0_dp*cy_i*uy_i(1) + xhat)*ri5damp) + &
               qzz_i*(cz2*xhatdot5 - (2.0_dp*cz_i*uz_i(1) + xhat)*ri5damp) ) 
          dudy = dudy - prei4 * ( &
               qxx_i*(cx2*yhatdot5 - (2.0_dp*cx_i*ux_i(2) + yhat)*ri5damp) + &
               qyy_i*(cy2*yhatdot5 - (2.0_dp*cy_i*uy_i(2) + yhat)*ri5damp) + &
               qzz_i*(cz2*yhatdot5 - (2.0_dp*cz_i*uz_i(2) + yhat)*ri5damp) ) 
          dudz = dudz - prei4 * ( &
               qxx_i*(cx2*zhatdot5 - (2.0_dp*cx_i*ux_i(3) + zhat)*ri5damp) + &
               qyy_i*(cy2*zhatdot5 - (2.0_dp*cy_i*uy_i(3) + zhat)*ri5damp) + &
               qzz_i*(cz2*zhatdot5 - (2.0_dp*cz_i*uz_i(3) + zhat)*ri5damp) ) 
          
          dudux_i(1) = dudux_i(1) + prei3*(qxx_i*cx_i*xhatdot2)
          dudux_i(2) = dudux_i(2) + prei3*(qxx_i*cx_i*yhatdot2)
          dudux_i(3) = dudux_i(3) + prei3*(qxx_i*cx_i*zhatdot2)
          
          duduy_i(1) = duduy_i(1) + prei3*(qyy_i*cy_i*xhatdot2)
          duduy_i(2) = duduy_i(2) + prei3*(qyy_i*cy_i*yhatdot2)
          duduy_i(3) = duduy_i(3) + prei3*(qyy_i*cy_i*zhatdot2)
          
          duduz_i(1) = duduz_i(1) + prei3*(qzz_i*cz_i*xhatdot2)
          duduz_i(2) = duduz_i(2) + prei3*(qzz_i*cz_i*yhatdot2)
          duduz_i(3) = duduz_i(3) + prei3*(qzz_i*cz_i*zhatdot2)
       endif
    endif


    if (do_pot) then
#ifdef IS_MPI 
       pot_row(ELECTROSTATIC_POT,atom1) = pot_row(ELECTROSTATIC_POT,atom1) + 0.5_dp*epot
       pot_col(ELECTROSTATIC_POT,atom2) = pot_col(ELECTROSTATIC_POT,atom2) + 0.5_dp*epot
#else
       pot = pot + epot
#endif
    endif

#ifdef IS_MPI
    f_Row(1,atom1) = f_Row(1,atom1) + dudx
    f_Row(2,atom1) = f_Row(2,atom1) + dudy
    f_Row(3,atom1) = f_Row(3,atom1) + dudz

    f_Col(1,atom2) = f_Col(1,atom2) - dudx
    f_Col(2,atom2) = f_Col(2,atom2) - dudy
    f_Col(3,atom2) = f_Col(3,atom2) - dudz

    if (i_is_Dipole .or. i_is_Quadrupole) then
       t_Row(1,atom1)=t_Row(1,atom1) - uz_i(2)*duduz_i(3) + uz_i(3)*duduz_i(2)
       t_Row(2,atom1)=t_Row(2,atom1) - uz_i(3)*duduz_i(1) + uz_i(1)*duduz_i(3)
       t_Row(3,atom1)=t_Row(3,atom1) - uz_i(1)*duduz_i(2) + uz_i(2)*duduz_i(1)
    endif
    if (i_is_Quadrupole) then
       t_Row(1,atom1)=t_Row(1,atom1) - ux_i(2)*dudux_i(3) + ux_i(3)*dudux_i(2)
       t_Row(2,atom1)=t_Row(2,atom1) - ux_i(3)*dudux_i(1) + ux_i(1)*dudux_i(3)
       t_Row(3,atom1)=t_Row(3,atom1) - ux_i(1)*dudux_i(2) + ux_i(2)*dudux_i(1)

       t_Row(1,atom1)=t_Row(1,atom1) - uy_i(2)*duduy_i(3) + uy_i(3)*duduy_i(2)
       t_Row(2,atom1)=t_Row(2,atom1) - uy_i(3)*duduy_i(1) + uy_i(1)*duduy_i(3)
       t_Row(3,atom1)=t_Row(3,atom1) - uy_i(1)*duduy_i(2) + uy_i(2)*duduy_i(1)
    endif

    if (j_is_Dipole .or. j_is_Quadrupole) then
       t_Col(1,atom2)=t_Col(1,atom2) - uz_j(2)*duduz_j(3) + uz_j(3)*duduz_j(2)
       t_Col(2,atom2)=t_Col(2,atom2) - uz_j(3)*duduz_j(1) + uz_j(1)*duduz_j(3)
       t_Col(3,atom2)=t_Col(3,atom2) - uz_j(1)*duduz_j(2) + uz_j(2)*duduz_j(1)
    endif
    if (j_is_Quadrupole) then
       t_Col(1,atom2)=t_Col(1,atom2) - ux_j(2)*dudux_j(3) + ux_j(3)*dudux_j(2)
       t_Col(2,atom2)=t_Col(2,atom2) - ux_j(3)*dudux_j(1) + ux_j(1)*dudux_j(3)
       t_Col(3,atom2)=t_Col(3,atom2) - ux_j(1)*dudux_j(2) + ux_j(2)*dudux_j(1)

       t_Col(1,atom2)=t_Col(1,atom2) - uy_j(2)*duduy_j(3) + uy_j(3)*duduy_j(2)
       t_Col(2,atom2)=t_Col(2,atom2) - uy_j(3)*duduy_j(1) + uy_j(1)*duduy_j(3)
       t_Col(3,atom2)=t_Col(3,atom2) - uy_j(1)*duduy_j(2) + uy_j(2)*duduy_j(1)
    endif

#else
    f(1,atom1) = f(1,atom1) + dudx
    f(2,atom1) = f(2,atom1) + dudy
    f(3,atom1) = f(3,atom1) + dudz

    f(1,atom2) = f(1,atom2) - dudx
    f(2,atom2) = f(2,atom2) - dudy
    f(3,atom2) = f(3,atom2) - dudz

    if (i_is_Dipole .or. i_is_Quadrupole) then
       t(1,atom1)=t(1,atom1) - uz_i(2)*duduz_i(3) + uz_i(3)*duduz_i(2)
       t(2,atom1)=t(2,atom1) - uz_i(3)*duduz_i(1) + uz_i(1)*duduz_i(3)
       t(3,atom1)=t(3,atom1) - uz_i(1)*duduz_i(2) + uz_i(2)*duduz_i(1)
    endif
    if (i_is_Quadrupole) then
       t(1,atom1)=t(1,atom1) - ux_i(2)*dudux_i(3) + ux_i(3)*dudux_i(2)
       t(2,atom1)=t(2,atom1) - ux_i(3)*dudux_i(1) + ux_i(1)*dudux_i(3)
       t(3,atom1)=t(3,atom1) - ux_i(1)*dudux_i(2) + ux_i(2)*dudux_i(1)

       t(1,atom1)=t(1,atom1) - uy_i(2)*duduy_i(3) + uy_i(3)*duduy_i(2)
       t(2,atom1)=t(2,atom1) - uy_i(3)*duduy_i(1) + uy_i(1)*duduy_i(3)
       t(3,atom1)=t(3,atom1) - uy_i(1)*duduy_i(2) + uy_i(2)*duduy_i(1)
    endif

    if (j_is_Dipole .or. j_is_Quadrupole) then
       t(1,atom2)=t(1,atom2) - uz_j(2)*duduz_j(3) + uz_j(3)*duduz_j(2)
       t(2,atom2)=t(2,atom2) - uz_j(3)*duduz_j(1) + uz_j(1)*duduz_j(3)
       t(3,atom2)=t(3,atom2) - uz_j(1)*duduz_j(2) + uz_j(2)*duduz_j(1)
    endif
    if (j_is_Quadrupole) then
       t(1,atom2)=t(1,atom2) - ux_j(2)*dudux_j(3) + ux_j(3)*dudux_j(2)
       t(2,atom2)=t(2,atom2) - ux_j(3)*dudux_j(1) + ux_j(1)*dudux_j(3)
       t(3,atom2)=t(3,atom2) - ux_j(1)*dudux_j(2) + ux_j(2)*dudux_j(1)

       t(1,atom2)=t(1,atom2) - uy_j(2)*duduy_j(3) + uy_j(3)*duduy_j(2)
       t(2,atom2)=t(2,atom2) - uy_j(3)*duduy_j(1) + uy_j(1)*duduy_j(3)
       t(3,atom2)=t(3,atom2) - uy_j(1)*duduy_j(2) + uy_j(2)*duduy_j(1)
    endif

#endif

#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif

    if (molMembershipList(id1) .ne. molMembershipList(id2)) then

       fpair(1) = fpair(1) + dudx
       fpair(2) = fpair(2) + dudy
       fpair(3) = fpair(3) + dudz

    endif

    return
  end subroutine doElectrostaticPair

  subroutine destroyElectrostaticTypes()

    if(allocated(ElectrostaticMap)) deallocate(ElectrostaticMap)

  end subroutine destroyElectrostaticTypes

  subroutine self_self(atom1, eFrame, mypot, t, do_pot)
    logical, intent(in) :: do_pot
    integer, intent(in) :: atom1
    integer :: atid1
    real(kind=dp), dimension(9,nLocal) :: eFrame
    real(kind=dp), dimension(3,nLocal) :: t
    real(kind=dp) :: mu1, c1
    real(kind=dp) :: preVal, epot, mypot
    real(kind=dp) :: eix, eiy, eiz

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
          c1 = getCharge(atid1)
          
          if (screeningMethod .eq. DAMPED) then
             mypot = mypot - (f0c * rcuti * 0.5_dp + &
                  dampingAlpha*invRootPi) * c1 * c1    
             
          else             
             mypot = mypot - (rcuti * 0.5_dp * c1 * c1)
             
          endif
       endif
    endif
    
    return
  end subroutine self_self

  subroutine rf_self_excludes(atom1, atom2, sw, eFrame, d, rij, vpair, myPot, &
       f, t, do_pot)
    logical, intent(in) :: do_pot
    integer, intent(in) :: atom1
    integer, intent(in) :: atom2
    logical :: i_is_Charge, j_is_Charge
    logical :: i_is_Dipole, j_is_Dipole
    integer :: atid1
    integer :: atid2
    real(kind=dp), intent(in) :: rij
    real(kind=dp), intent(in) :: sw
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
       
       preVal = pre11 * q_i * q_j
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
       
       pref = pre12 * q_i * mu_j
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
       
       pref = pre12 * q_j * mu_i
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

end module electrostatic_module
