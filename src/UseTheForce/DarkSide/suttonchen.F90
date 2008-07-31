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

!! Impliments Sutton-Chen Metallic Potential
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
#ifdef IS_MPI
  use mpiSimulation
#endif
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
  

  !! Arrays for derivatives used in force calculation
  real( kind = dp), dimension(:), allocatable :: rho
  real( kind = dp), dimension(:), allocatable :: frho
  real( kind = dp), dimension(:), allocatable :: dfrhodrho
  
  !! Arrays for MPI storage
#ifdef IS_MPI
  real( kind = dp),save, dimension(:), allocatable :: dfrhodrho_col
  real( kind = dp),save, dimension(:), allocatable :: dfrhodrho_row
  real( kind = dp),save, dimension(:), allocatable :: frho_row
  real( kind = dp),save, dimension(:), allocatable :: frho_col
  real( kind = dp),save, dimension(:), allocatable :: rho_row
  real( kind = dp),save, dimension(:), allocatable :: rho_col
  real( kind = dp),save, dimension(:), allocatable :: rho_tmp
#endif

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
  public :: clean_SC
  public :: destroySCtypes
  public :: getSCCut
 ! public :: setSCDefaultCutoff
 ! public :: setSCUniformCutoff
 

contains


  subroutine newSCtype(c_ident,c,m,n,alpha,epsilon,status)
    real (kind = dp )                      :: c ! Density Scaling
    real (kind = dp )                      :: m ! Density Exponent
    real (kind = dp )                      :: n ! Pair Potential Exponent
    real (kind = dp )                      :: alpha ! Length Scaling
    real (kind = dp )                      :: epsilon ! Energy Scaling
    integer                                :: c_ident
    integer                                :: status
    integer                                :: nAtypes,nSCTypes,myATID
    integer                                :: maxVals
    integer                                :: alloc_stat
    integer                                :: current
    integer,pointer                        :: Matchlist(:) => null()

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
             alpha = sqrt(alpha1 * alpha2) !SC formulation
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
  

  !! routine checks to see if array is allocated, deallocates array if allocated
  !! and then creates the array to the required size
  subroutine allocateSC()
    integer :: status

#ifdef IS_MPI
    integer :: nAtomsInRow
    integer :: nAtomsInCol
#endif
    integer :: alloc_stat

    
    status = 0
#ifdef IS_MPI
    nAtomsInRow = getNatomsInRow(plan_atom_row)
    nAtomsInCol = getNatomsInCol(plan_atom_col)
#endif

    if (allocated(frho)) deallocate(frho)
    allocate(frho(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0) then
       status = -1
    end if

    if (allocated(rho)) deallocate(rho)
    allocate(rho(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
    end if

    if (allocated(dfrhodrho)) deallocate(dfrhodrho)
    allocate(dfrhodrho(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
    end if

#ifdef IS_MPI

    if (allocated(rho_tmp)) deallocate(rho_tmp)
    allocate(rho_tmp(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
    end if


    if (allocated(frho_row)) deallocate(frho_row)
    allocate(frho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
    end if
    if (allocated(rho_row)) deallocate(rho_row)
    allocate(rho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
    end if
    if (allocated(dfrhodrho_row)) deallocate(dfrhodrho_row)
    allocate(dfrhodrho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
    end if


    ! Now do column arrays

    if (allocated(frho_col)) deallocate(frho_col)
    allocate(frho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
    end if
    if (allocated(rho_col)) deallocate(rho_col)
    allocate(rho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
    end if
    if (allocated(dfrhodrho_col)) deallocate(dfrhodrho_col)
    allocate(dfrhodrho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
    end if

#endif
    if (status == -1) then
       call handleError("SuttonChen:allocateSC","Error in allocating SC arrays")
    end if
    arraysAllocated = .true.
  end subroutine allocateSC

  subroutine setCutoffSC(rcut)
    real(kind=dp) :: rcut
    SC_rcut = rcut
  end subroutine setCutoffSC
  
  !! This array allocates module arrays if needed and builds mixing map.
  subroutine clean_SC()
    if (.not.arraysAllocated) call allocateSC()
    if (.not.haveMixingMap) call createMixingMap()
    ! clean non-IS_MPI first
    frho = 0.0_dp
    rho  = 0.0_dp
    dfrhodrho = 0.0_dp
    ! clean MPI if needed
#ifdef IS_MPI
    frho_row = 0.0_dp
    frho_col = 0.0_dp
    rho_row  = 0.0_dp
    rho_col  = 0.0_dp
    rho_tmp  = 0.0_dp
    dfrhodrho_row = 0.0_dp
    dfrhodrho_col = 0.0_dp
#endif
  end subroutine clean_SC

  !! Calculates rho_r
  subroutine calc_sc_prepair_rho(atom1, atom2, d, r, rijsq, rcut)
    integer :: atom1,atom2
    real(kind = dp), dimension(3) :: d
    real(kind = dp), intent(inout)               :: r
    real(kind = dp), intent(inout)               :: rijsq, rcut
    ! value of electron density rho do to atom i at atom j
    real(kind = dp) :: rho_i_at_j
    ! value of electron density rho do to atom j at atom i
    real(kind = dp) :: rho_j_at_i

    ! we don't use the derivatives, dummy variables
    real( kind = dp) :: drho
    integer :: sc_err
    
    integer :: atid1,atid2 ! Global atid    
    integer :: myid_atom1 ! SC atid
    integer :: myid_atom2 


    ! check to see if we need to be cleaned at the start of a force loop

    if (cleanArrays) call clean_SC()
    cleanArrays = .false.
    
    

#ifdef IS_MPI
    Atid1 = Atid_row(Atom1)
    Atid2 = Atid_col(Atom2)
#else
    Atid1 = Atid(Atom1)
    Atid2 = Atid(Atom2)
#endif

    Myid_atom1 = SCList%atidtoSCtype(Atid1)
    Myid_atom2 = SCList%atidtoSCtype(Atid2)
    
    call lookupUniformSpline(MixingMap(myid_atom1,myid_atom2)%phi, r, &
         rho_i_at_j)
    rho_j_at_i = rho_i_at_j

#ifdef  IS_MPI
    rho_col(atom2) = rho_col(atom2) + rho_i_at_j
    rho_row(atom1) = rho_row(atom1) + rho_j_at_i
#else
    rho(atom2) = rho(atom2) + rho_i_at_j
    rho(atom1) = rho(atom1) + rho_j_at_i
#endif
    
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

    !! Scatter the electron density from  pre-pair calculation back to 
    !! local atoms

    
    if (cleanArrays) call clean_SC()
    cleanArrays = .false.
    
      
#ifdef IS_MPI
    call scatter(rho_row,rho,plan_atom_row,sc_err)
    if (sc_err /= 0 ) then
       write(errMsg,*) " Error scattering rho_row into rho"
       call handleError(RoutineName,errMesg)
    endif
    call scatter(rho_col,rho_tmp,plan_atom_col,sc_err)
    if (sc_err /= 0 ) then
       write(errMsg,*) " Error scattering rho_col into rho"
       call handleError(RoutineName,errMesg)

    endif

    rho(1:nlocal) = rho(1:nlocal) + rho_tmp(1:nlocal)
#endif
    
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

#ifdef IS_MPI
    !! communicate f(rho) and derivatives back into row and column arrays
    call gather(frho,frho_row,plan_atom_row, sc_err)
    if (sc_err /=  0) then
       call handleError("calc_sc_forces()","MPI gather frho_row failure")
    endif
    call gather(dfrhodrho,dfrhodrho_row,plan_atom_row, sc_err)
    if (sc_err /=  0) then
       call handleError("calc_sc_forces()","MPI gather dfrhodrho_row failure")
    endif
    call gather(frho,frho_col,plan_atom_col, sc_err)
    if (sc_err /=  0) then
       call handleError("calc_sc_forces()","MPI gather frho_col failure")
    endif
    call gather(dfrhodrho,dfrhodrho_col,plan_atom_col, sc_err)
    if (sc_err /=  0) then
       call handleError("calc_sc_forces()","MPI gather dfrhodrho_col failure")
    endif
#endif
    
    
  end subroutine calc_sc_preforce_Frho  
  
  !! Does Sutton-Chen  pairwise Force calculation.  
  subroutine do_sc_pair(atom1, atom2, d, rij, r2, rcut, sw, vpair, fpair, &
       pot, f, do_pot)
    !Arguments    
    integer, intent(in) ::  atom1, atom2
    real( kind = dp ), intent(in) :: rij, r2, rcut
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair

    logical, intent(in) :: do_pot

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
    integer  :: atid1, atid2
    !Local Variables
    
    cleanArrays = .true.

#ifdef IS_MPI
       atid1 = atid_row(atom1)
       atid2 = atid_col(atom2)
#else
       atid1 = atid(atom1)
       atid2 = atid(atom2)
#endif

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
       
#ifdef IS_MPI
       dudr = drhodr*(dfrhodrho_row(atom1) + dfrhodrho_col(atom2)) + dvpdr
#else
       dudr = drhodr*(dfrhodrho(atom1)+ dfrhodrho(atom2)) + dvpdr
#endif
              
       pot_temp = vptmp - vcij

       vpair = vpair + pot_temp
 
       fx = dudr * drdx
       fy = dudr * drdy
       fz = dudr * drdz

#ifdef IS_MPI
       if (do_pot) then
          pot_Row(METALLIC_POT,atom1) = pot_Row(METALLIC_POT,atom1) + (pot_temp)*0.5
          pot_Col(METALLIC_POT,atom2) = pot_Col(METALLIC_POT,atom2) + (pot_temp)*0.5
       end if

       f_Row(1,atom1) = f_Row(1,atom1) + fx
       f_Row(2,atom1) = f_Row(2,atom1) + fy
       f_Row(3,atom1) = f_Row(3,atom1) + fz

       f_Col(1,atom2) = f_Col(1,atom2) - fx
       f_Col(2,atom2) = f_Col(2,atom2) - fy
       f_Col(3,atom2) = f_Col(3,atom2) - fz
#else

       if(do_pot) then
          pot = pot + pot_temp
       end if

       f(1,atom1) = f(1,atom1) + fx
       f(2,atom1) = f(2,atom1) + fy
       f(3,atom1) = f(3,atom1) + fz

       f(1,atom2) = f(1,atom2) - fx
       f(2,atom2) = f(2,atom2) - fy
       f(3,atom2) = f(3,atom2) - fz
#endif

       
#ifdef IS_MPI
       id1 = AtomRowToGlobal(atom1)
       id2 = AtomColToGlobal(atom2)
#else
       id1 = atom1
       id2 = atom2
#endif

       if (molMembershipList(id1) .ne. molMembershipList(id2)) then

          fpair(1) = fpair(1) + fx
          fpair(2) = fpair(2) + fy
          fpair(3) = fpair(3) + fz

       endif
  end subroutine do_sc_pair
end module suttonchen
