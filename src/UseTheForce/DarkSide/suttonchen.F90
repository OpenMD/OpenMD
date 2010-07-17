!!
!! Copyright (c) 2005, 2009 The University of Notre Dame. All Rights Reserved.
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

!! Implements Sutton-Chen Metallic Potential
!! See A.P.SUTTON and J.CHEN,PHIL MAG LETT 61,139-146,1990


module suttonchen
  use definitions
  use simulation
  use force_globals
  use status
  use atype_module
  use vector_class
  use fForceOptions
  use interpolation
  implicit none
  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"
  
  !! number of points for the spline approximations
  INTEGER, PARAMETER :: np = 3000
  
  logical, save :: SC_FF_initialized = .false.
  integer, save :: SC_Mixing_Policy
  real(kind = dp), save :: SC_rcut
  logical, save :: haveRcut = .false.
  logical, save :: haveMixingMap = .false.
  logical, save :: useGeometricDistanceMixing = .false.
  logical, save :: cleanArrays = .true.
  logical, save :: arraysAllocated = .false.
  

  character(len = statusMsgSize) :: errMesg
  integer :: sc_err
  
  character(len = 200) :: errMsg
  character(len=*), parameter :: RoutineName =  "Sutton-Chen MODULE"
  
  type, private :: SCtype
     integer       :: atid       
     real(kind=dp) :: c 
     real(kind=dp) :: m
     real(kind=dp) :: n
     real(kind=dp) :: alpha
     real(kind=dp) :: epsilon
     real(kind=dp) :: sc_rcut
  end type SCtype
  
  
  type, private :: SCTypeList
     integer           :: nSCTypes = 0
     integer           :: currentSCtype = 0     
     type (SCtype), pointer :: SCtypes(:) => null()
     integer, pointer       :: atidToSCtype(:) => null()
  end type SCTypeList
  
  type (SCTypeList), save :: SCList
  
  type:: MixParameters
     real(kind=DP) :: alpha
     real(kind=DP) :: epsilon
     real(kind=DP) :: m 
     real(Kind=DP) :: n
     real(kind=dp) :: rCut
     real(kind=dp) :: vCut
     logical       :: rCutWasSet = .false.
     type(cubicSpline) :: V
     type(cubicSpline) :: phi
  end type MixParameters
  
  type(MixParameters), dimension(:,:), allocatable :: MixingMap
  
  public :: setCutoffSC
  public :: do_SC_pair
  public :: newSCtype
  public :: calc_SC_prepair_rho
  public :: calc_SC_preforce_Frho
  public :: destroySCtypes
  public :: getSCCut
  ! public :: setSCDefaultCutoff
  ! public :: setSCUniformCutoff
  
  
contains
  
  
  subroutine newSCtype(c_ident,c,m,n,alpha,epsilon,status)
    real (kind = dp ) :: c ! Density Scaling
    real (kind = dp ) :: m ! Density Exponent
    real (kind = dp ) :: n ! Pair Potential Exponent
    real (kind = dp ) :: alpha ! Length Scaling
    real (kind = dp ) :: epsilon ! Energy Scaling
    integer           :: c_ident
    integer           :: status
    integer           :: nAtypes,nSCTypes,myATID
    integer           :: maxVals
    integer           :: alloc_stat
    integer           :: current
    integer,pointer   :: Matchlist(:) => null()
    
    status = 0
    
    
    !! Assume that atypes has already been set and get the total number of types in atypes
    
    ! check to see if this is the first time into 
    if (.not.associated(SCList%SCTypes)) then
       call getMatchingElementList(atypes, "is_SC", .true., nSCtypes, MatchList)
       SCList%nSCtypes = nSCtypes
       allocate(SCList%SCTypes(nSCTypes))
       nAtypes = getSize(atypes)
       allocate(SCList%atidToSCType(nAtypes))
       SCList%atidToSCType = -1
    end if
    
    SCList%currentSCType = SCList%currentSCType + 1
    current = SCList%currentSCType
    
    myATID =  getFirstMatchingElement(atypes, "c_ident", c_ident)
    SCList%atidToSCType(myATID) = current
    
    
    SCList%SCTypes(current)%atid         = c_ident
    SCList%SCTypes(current)%alpha        = alpha
    SCList%SCTypes(current)%c            = c
    SCList%SCTypes(current)%m            = m
    SCList%SCTypes(current)%n            = n
    SCList%SCTypes(current)%epsilon      = epsilon
  end subroutine newSCtype
  
  
  subroutine destroySCTypes()
    if (associated(SCList%SCtypes)) then
       deallocate(SCList%SCTypes)
       SCList%SCTypes=>null()
    end if
    if (associated(SCList%atidToSCtype)) then
       deallocate(SCList%atidToSCtype)
       SCList%atidToSCtype=>null()
    end if
    ! Reset Capacity
    SCList%nSCTypes = 0
    SCList%currentSCtype=0
    
  end subroutine destroySCTypes
  
  function getSCCut(atomID) result(cutValue)
    integer, intent(in) :: atomID
    integer :: scID
    real(kind=dp) :: cutValue
    
    scID = SCList%atidToSCType(atomID)
    cutValue = 2.0_dp * SCList%SCTypes(scID)%alpha
  end function getSCCut
  
  subroutine createMixingMap()
    integer :: nSCtypes, i, j, k
    real ( kind = dp ) :: e1, e2, m1, m2, alpha1, alpha2, n1, n2
    real ( kind = dp ) :: epsilon, m, n, alpha, rCut, vCut, dr, r
    real ( kind = dp ), dimension(np) :: rvals, vvals, phivals
    
    if (SCList%currentSCtype == 0) then
       call handleError("SuttonChen", "No members in SCMap")
       return
    end if
    
    nSCtypes = SCList%nSCtypes
    
    if (.not. allocated(MixingMap)) then
       allocate(MixingMap(nSCtypes, nSCtypes))
    endif
    useGeometricDistanceMixing = usesGeometricDistanceMixing()
    do i = 1, nSCtypes
       
       e1 = SCList%SCtypes(i)%epsilon
       m1 = SCList%SCtypes(i)%m
       n1 = SCList%SCtypes(i)%n
       alpha1 = SCList%SCtypes(i)%alpha
       
       do j = 1, nSCtypes
          
          e2 = SCList%SCtypes(j)%epsilon
          m2 = SCList%SCtypes(j)%m
          n2 = SCList%SCtypes(j)%n
          alpha2 = SCList%SCtypes(j)%alpha
          
          if (useGeometricDistanceMixing) then
             alpha = sqrt(alpha1 * alpha2) ! SC formulation
          else
             alpha = 0.5_dp * (alpha1 + alpha2) ! Goddard formulation
          endif
          rcut = 2.0_dp * alpha
          epsilon = sqrt(e1 * e2)
          m = 0.5_dp*(m1+m2)
          n = 0.5_dp*(n1+n2)
          
          dr = (rCut) / dble(np-1)
          rvals(1) = 0.0_dp
          vvals(1) = 0.0_dp
          phivals(1) = 0.0_dp
          
          do k = 2, np
             r = dble(k-1)*dr
             rvals(k) = r
             vvals(k) = epsilon*((alpha/r)**n)
             phivals(k) = (alpha/r)**m
          enddo
          
          vCut = epsilon*((alpha/rCut)**n)
          
          call newSpline(MixingMap(i,j)%V, rvals, vvals, .true.)
          call newSpline(MixingMap(i,j)%phi, rvals, phivals, .true.)
          
          MixingMap(i,j)%epsilon = epsilon
          MixingMap(i,j)%m = m
          MixingMap(i,j)%n = n
          MixingMap(i,j)%alpha = alpha
          MixingMap(i,j)%rCut = rcut
          MixingMap(i,j)%vCut = vCut
       enddo
    enddo
    
    haveMixingMap = .true.
    
  end subroutine createMixingMap
  

  
  subroutine setCutoffSC(rcut)
    real(kind=dp) :: rcut
    SC_rcut = rcut
  end subroutine setCutoffSC
  
  
  !! Calculates rho_r
  subroutine calc_sc_prepair_rho(atid1, atid2, d, r, rijsq, rho_i_at_j, rho_j_at_i)
    integer :: atid1, atid2
    real(kind = dp), dimension(3) :: d
    real(kind = dp), intent(inout)               :: r
    real(kind = dp), intent(inout)               :: rijsq
    ! value of electron density rho do to atom i at atom j
    real(kind = dp), intent(inout) :: rho_i_at_j
    ! value of electron density rho do to atom j at atom i
    real(kind = dp), intent(inout) :: rho_j_at_i
    ! we don't use the derivatives, dummy variables
    real( kind = dp) :: drho
    integer :: sc_err
    
    integer :: myid_atom1 ! SC atid
    integer :: myid_atom2 
    
    ! check to see if we need to be cleaned at the start of a force loop
    
    if (.not.haveMixingMap) call createMixingMap()
    haveMixingMap = .true.
    
    Myid_atom1 = SCList%atidtoSCtype(Atid1)
    Myid_atom2 = SCList%atidtoSCtype(Atid2)
    
    call lookupUniformSpline(MixingMap(myid_atom1,myid_atom2)%phi, r, &
         rho_i_at_j)
    rho_j_at_i = rho_i_at_j
       
  end subroutine calc_sc_prepair_rho
  
  
  !! Calculate the rho_a for all local atoms
  subroutine calc_sc_preforce_Frho(nlocal, pot, particle_pot)
    integer :: nlocal
    real(kind=dp) :: pot
    real(kind=dp), dimension(nlocal) :: particle_pot
    integer :: i,j
    integer :: atom
    integer :: atype1
    integer :: atid1
    integer :: myid
    
    !! Calculate F(rho) and derivative
    do atom = 1, nlocal
       Myid = SCList%atidtoSctype(Atid(atom))
       ! Myid is set to -1 for non SC atoms.
       ! Punt if we are a non-SC atom type.
       if (Myid == -1) then
          frho(atom) = 0.0_dp
          dfrhodrho(atom) = 0.0_dp
       else
          frho(atom) = - SCList%SCTypes(Myid)%c * &
               SCList%SCTypes(Myid)%epsilon * sqrt(rho(atom))
          
          dfrhodrho(atom) = 0.5_dp*frho(atom)/rho(atom)
       end if
       pot = pot + frho(atom)
       particle_pot(atom) = particle_pot(atom) + frho(atom)
    enddo
    
  end subroutine calc_sc_preforce_Frho
  
  !! Does Sutton-Chen  pairwise Force calculation.  
  subroutine do_sc_pair(atid1, atid2, d, rij, r2, sw, vpair, &
       pot, f1, rho_i, rho_j, dfrhodrho_i, dfrhodrho_j, &
       fshift_i, fshift_j)
    !Arguments    
    integer, intent(in) ::  atid1, atid2
    real( kind = dp ), intent(in) :: rij, r2
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), dimension(3) :: f1
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout) :: dfrhodrho_i, dfrhodrho_j 
    real( kind = dp ), intent(inout) :: rho_i, rho_j 
    real( kind = dp ), intent(inout):: fshift_i, fshift_j   
    
    real( kind = dp ) :: drdx, drdy, drdz
    real( kind = dp ) :: dvpdr
    real( kind = dp ) :: rhtmp, drhodr
    real( kind = dp ) :: dudr
    real( kind = dp ) :: Fx,Fy,Fz
    real( kind = dp ) :: pot_temp, vptmp
    real( kind = dp ) :: rcij, vcij
    
    integer :: id1, id2
    integer  :: mytype_atom1
    integer  :: mytype_atom2
    
    mytype_atom1 = SCList%atidToSCType(atid1)
    mytype_atom2 = SCList%atidTOSCType(atid2)
    
    drdx = d(1)/rij
    drdy = d(2)/rij
    drdz = d(3)/rij
    
    rcij = MixingMap(mytype_atom1,mytype_atom2)%rCut
    vcij = MixingMap(mytype_atom1,mytype_atom2)%vCut
    
    call lookupUniformSpline1d(MixingMap(mytype_atom1, mytype_atom2)%phi, &
         rij, rhtmp, drhodr)
    
    call lookupUniformSpline1d(MixingMap(mytype_atom1, mytype_atom2)%V, &
         rij, vptmp, dvpdr)
    
    dudr = drhodr*(dfrhodrho_i + dfrhodrho_j) + dvpdr

    pot_temp = vptmp - vcij
    
    vpair = vpair + pot_temp
    
    fx = dudr * drdx
    fy = dudr * drdy
    fz = dudr * drdz
        
    ! particle_pot is the difference between the full potential 
    ! and the full potential without the presence of a particular
    ! particle (atom1).
    !
    ! This reduces the density at other particle locations, so
    ! we need to recompute the density at atom2 assuming atom1
    ! didn't contribute.  This then requires recomputing the
    ! density functional for atom2 as well.
    !
    ! Most of the particle_pot heavy lifting comes from the
    ! pair interaction, and will be handled by vpair.
    
    fshift_i = -SCList%SCTypes(mytype_atom1)%c * &
         SCList%SCTypes(mytype_atom1)%epsilon * sqrt(rho_i-rhtmp)
    fshift_j = -SCList%SCTypes(mytype_atom2)%c * &
         SCList%SCTypes(mytype_atom2)%epsilon * sqrt(rho_j-rhtmp)         
    
    pot = pot + pot_temp
    
    f1(1) = f1(1) + fx
    f1(2) = f1(2) + fy
    f1(3) = f1(3) + fz
    
  end subroutine do_sc_pair
end module suttonchen
