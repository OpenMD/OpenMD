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

module eam
  use simulation
  use force_globals
  use status
  use atype_module
  use vector_class
#ifdef IS_MPI
  use mpiSimulation
#endif
  implicit none
  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"

  INTEGER, PARAMETER :: DP = selected_real_kind(15)

  logical, save :: EAM_FF_initialized = .false.
  integer, save :: EAM_Mixing_Policy
  real(kind = dp), save :: EAM_rcut
  logical, save :: haveRcut = .false.

  character(len = statusMsgSize) :: errMesg
  integer :: eam_err

  character(len = 200) :: errMsg
  character(len=*), parameter :: RoutineName =  "EAM MODULE"
  !! Logical that determines if eam arrays should be zeroed
  logical :: cleanme = .true.
  logical :: nmflag  = .false.


  type, private :: EAMtype
     integer           :: eam_atype       
     real( kind = DP ) :: eam_dr          
     integer           :: eam_nr           
     integer           :: eam_nrho          
     real( kind = DP ) :: eam_lattice        
     real( kind = DP ) :: eam_drho      
     real( kind = DP ) :: eam_rcut     
     integer           :: eam_atype_map

     real( kind = DP ), pointer, dimension(:) :: eam_rvals        => null()
     real( kind = DP ), pointer, dimension(:) :: eam_rhovals      => null()
     real( kind = DP ), pointer, dimension(:) :: eam_F_rho        => null()
     real( kind = DP ), pointer, dimension(:) :: eam_Z_r          => null()
     real( kind = DP ), pointer, dimension(:) :: eam_rho_r        => null()
     real( kind = DP ), pointer, dimension(:) :: eam_phi_r        => null()
     real( kind = DP ), pointer, dimension(:) :: eam_F_rho_pp     => null()
     real( kind = DP ), pointer, dimension(:) :: eam_Z_r_pp       => null()
     real( kind = DP ), pointer, dimension(:) :: eam_rho_r_pp     => null()
     real( kind = DP ), pointer, dimension(:) :: eam_phi_r_pp     => null()
  end type EAMtype


  !! Arrays for derivatives used in force calculation
  real( kind = dp), dimension(:), allocatable :: frho
  real( kind = dp), dimension(:), allocatable :: rho

  real( kind = dp), dimension(:), allocatable :: dfrhodrho
  real( kind = dp), dimension(:), allocatable :: d2frhodrhodrho


  !! Arrays for MPI storage
#ifdef IS_MPI
  real( kind = dp),save, dimension(:), allocatable :: dfrhodrho_col
  real( kind = dp),save, dimension(:), allocatable :: dfrhodrho_row
  real( kind = dp),save, dimension(:), allocatable :: frho_row
  real( kind = dp),save, dimension(:), allocatable :: frho_col
  real( kind = dp),save, dimension(:), allocatable :: rho_row
  real( kind = dp),save, dimension(:), allocatable :: rho_col
  real( kind = dp),save, dimension(:), allocatable :: rho_tmp
  real( kind = dp),save, dimension(:), allocatable :: d2frhodrhodrho_col
  real( kind = dp),save, dimension(:), allocatable :: d2frhodrhodrho_row
#endif

  type, private :: EAMTypeList
     integer           :: n_eam_types = 0
     integer           :: currentAddition = 0

     type (EAMtype), pointer  :: EAMParams(:) => null()
     integer, pointer         :: atidToEAMType(:) => null()
  end type EAMTypeList


  type (eamTypeList), save :: EAMList

  !! standard eam stuff  


  public :: init_EAM_FF
  public :: setCutoffEAM
  public :: do_eam_pair
  public :: newEAMtype
  public :: calc_eam_prepair_rho
  public :: calc_eam_preforce_Frho
  public :: clean_EAM
  public :: destroyEAMTypes
  public :: getEAMCut

contains


  subroutine newEAMtype(lattice_constant,eam_nrho,eam_drho,eam_nr,&
       eam_dr,rcut,eam_Z_r,eam_rho_r,eam_F_rho,&
       c_ident,status)
    real (kind = dp )                      :: lattice_constant
    integer                                :: eam_nrho
    real (kind = dp )                      :: eam_drho
    integer                                :: eam_nr
    real (kind = dp )                      :: eam_dr
    real (kind = dp )                      :: rcut
    real (kind = dp ), dimension(eam_nr)   :: eam_Z_r
    real (kind = dp ), dimension(eam_nr)   :: eam_rho_r
    real (kind = dp ), dimension(eam_nrho) :: eam_F_rho
    integer                                :: c_ident
    integer                                :: status

    integer                                :: nAtypes,nEAMTypes,myATID
    integer                                :: maxVals
    integer                                :: alloc_stat
    integer                                :: current
    integer,pointer                        :: Matchlist(:) => null()

    status = 0


    !! Assume that atypes has already been set and get the total number of types in atypes
    !! Also assume that every member of atypes is a EAM model.


    ! check to see if this is the first time into 
    if (.not.associated(EAMList%EAMParams)) then
       call getMatchingElementList(atypes, "is_EAM", .true., nEAMtypes, MatchList)
       EAMList%n_eam_types = nEAMtypes
       allocate(EAMList%EAMParams(nEAMTypes))
       nAtypes = getSize(atypes)
       allocate(EAMList%atidToEAMType(nAtypes))
    end if

    EAMList%currentAddition = EAMList%currentAddition + 1
    current = EAMList%currentAddition

    myATID =  getFirstMatchingElement(atypes, "c_ident", c_ident)
    EAMList%atidToEAMType(myATID) = current

    call allocate_EAMType(eam_nrho,eam_nr,EAMList%EAMParams(current),stat=alloc_stat)
    if (alloc_stat /= 0) then
       status = -1
       return
    end if

  
    EAMList%EAMParams(current)%eam_atype    = c_ident
    EAMList%EAMParams(current)%eam_lattice  = lattice_constant
    EAMList%EAMParams(current)%eam_nrho     = eam_nrho
    EAMList%EAMParams(current)%eam_drho     = eam_drho
    EAMList%EAMParams(current)%eam_nr       = eam_nr
    EAMList%EAMParams(current)%eam_dr       = eam_dr
    EAMList%EAMParams(current)%eam_rcut     = rcut
    EAMList%EAMParams(current)%eam_Z_r      = eam_Z_r
    EAMList%EAMParams(current)%eam_rho_r    = eam_rho_r
    EAMList%EAMParams(current)%eam_F_rho    = eam_F_rho

  end subroutine newEAMtype


  ! kills all eam types entered and sets simulation to uninitalized
  subroutine destroyEAMtypes()
    integer :: i
    type(EAMType), pointer :: tempEAMType=>null()

    do i = 1, EAMList%n_eam_types
       tempEAMType => eamList%EAMParams(i)
       call deallocate_EAMType(tempEAMType)
    end do
    if(associated( eamList%EAMParams)) deallocate( eamList%EAMParams)
    eamList%EAMParams => null()

    eamList%n_eam_types = 0
    eamList%currentAddition = 0

  end subroutine destroyEAMtypes

  function getEAMCut(atomID) result(cutValue)
    integer, intent(in) :: atomID
    integer :: eamID
    real(kind=dp) :: cutValue
    
    eamID = EAMList%atidToEAMType(atomID)
    cutValue = EAMList%EAMParams(eamID)%eam_rcut
  end function getEAMCut

  subroutine init_EAM_FF(status)
    integer :: status
    integer :: i,j
    real(kind=dp) :: current_rcut_max
    integer :: alloc_stat
    integer :: number_r, number_rho


    status = 0
    if (EAMList%currentAddition == 0) then
       call handleError("init_EAM_FF","No members in EAMList")
       status = -1
       return
    end if


    do i = 1, EAMList%currentAddition

       ! Build array of r values

       do j = 1,EAMList%EAMParams(i)%eam_nr
          EAMList%EAMParams(i)%eam_rvals(j) = &
               real(j-1,kind=dp)* &
               EAMList%EAMParams(i)%eam_dr
       end do
       ! Build array of rho values
       do j = 1,EAMList%EAMParams(i)%eam_nrho
          EAMList%EAMParams(i)%eam_rhovals(j) = &
               real(j-1,kind=dp)* &
               EAMList%EAMParams(i)%eam_drho
       end do
       ! convert from eV to kcal / mol:
       EAMList%EAMParams(i)%eam_F_rho = EAMList%EAMParams(i)%eam_F_rho * 23.06054E0_DP

       ! precompute the pair potential and get it into kcal / mol:
       EAMList%EAMParams(i)%eam_phi_r(1) = 0.0E0_DP
       do j = 2, EAMList%EAMParams(i)%eam_nr
          EAMList%EAMParams(i)%eam_phi_r(j) = (EAMList%EAMParams(i)%eam_Z_r(j)**2)/EAMList%EAMParams(i)%eam_rvals(j)
          EAMList%EAMParams(i)%eam_phi_r(j) =  EAMList%EAMParams(i)%eam_phi_r(j)*331.999296E0_DP
       enddo
    end do


    do i = 1,  EAMList%currentAddition
       number_r   = EAMList%EAMParams(i)%eam_nr
       number_rho = EAMList%EAMParams(i)%eam_nrho

       call eam_spline(number_r, EAMList%EAMParams(i)%eam_rvals, &
            EAMList%EAMParams(i)%eam_rho_r, &
            EAMList%EAMParams(i)%eam_rho_r_pp, &
            0.0E0_DP, 0.0E0_DP, 'N')
       call eam_spline(number_r, EAMList%EAMParams(i)%eam_rvals, &
            EAMList%EAMParams(i)%eam_Z_r, &
            EAMList%EAMParams(i)%eam_Z_r_pp, &
            0.0E0_DP, 0.0E0_DP, 'N')
       call eam_spline(number_rho, EAMList%EAMParams(i)%eam_rhovals, &
            EAMList%EAMParams(i)%eam_F_rho, &
            EAMList%EAMParams(i)%eam_F_rho_pp, &
            0.0E0_DP, 0.0E0_DP, 'N')
       call eam_spline(number_r, EAMList%EAMParams(i)%eam_rvals, &
            EAMList%EAMParams(i)%eam_phi_r, &
            EAMList%EAMParams(i)%eam_phi_r_pp, &
            0.0E0_DP, 0.0E0_DP, 'N')
    enddo

    !       current_rcut_max = EAMList%EAMParams(1)%eam_rcut
    !! find the smallest rcut for any eam atype
    !       do i = 2, EAMList%currentAddition 
    !          current_rcut_max =max(current_rcut_max,EAMList%EAMParams(i)%eam_rcut)
    !       end do

    !       EAM_rcut = current_rcut_max
    !       EAM_rcut_orig = current_rcut_max
    !       do i = 1, EAMList%currentAddition
    !          EAMList%EAMParam(i)s%eam_atype_map(eam_atype(i)) = i
    !       end do
    !! Allocate arrays for force calculation

    call allocateEAM(alloc_stat)
    if (alloc_stat /= 0 ) then
       write(*,*) "allocateEAM failed"
       status = -1
       return
    endif

  end subroutine init_EAM_FF

  !! routine checks to see if array is allocated, deallocates array if allocated
  !! and then creates the array to the required size
  subroutine allocateEAM(status)
    integer, intent(out) :: status

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
       return
    end if
    if (allocated(rho)) deallocate(rho)
    allocate(rho(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if

    if (allocated(dfrhodrho)) deallocate(dfrhodrho)
    allocate(dfrhodrho(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if

    if (allocated(d2frhodrhodrho)) deallocate(d2frhodrhodrho)
    allocate(d2frhodrhodrho(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if

#ifdef IS_MPI

    if (allocated(rho_tmp)) deallocate(rho_tmp)
    allocate(rho_tmp(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if


    if (allocated(frho_row)) deallocate(frho_row)
    allocate(frho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if
    if (allocated(rho_row)) deallocate(rho_row)
    allocate(rho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if
    if (allocated(dfrhodrho_row)) deallocate(dfrhodrho_row)
    allocate(dfrhodrho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if
    if (allocated(d2frhodrhodrho_row)) deallocate(d2frhodrhodrho_row)
    allocate(d2frhodrhodrho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if


    ! Now do column arrays

    if (allocated(frho_col)) deallocate(frho_col)
    allocate(frho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if
    if (allocated(rho_col)) deallocate(rho_col)
    allocate(rho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if
    if (allocated(dfrhodrho_col)) deallocate(dfrhodrho_col)
    allocate(dfrhodrho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if
    if (allocated(d2frhodrhodrho_col)) deallocate(d2frhodrhodrho_col)
    allocate(d2frhodrhodrho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       status = -1
       return
    end if

#endif

  end subroutine allocateEAM

  !! C sets rcut to be the largest cutoff of any atype 
  !! present in this simulation. Doesn't include all atypes
  !! sim knows about, just those in the simulation.
  subroutine setCutoffEAM(rcut, status)
    real(kind=dp) :: rcut
    integer :: status
    status = 0

    EAM_rcut = rcut

  end subroutine setCutoffEAM



  subroutine clean_EAM()

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
  end subroutine clean_EAM



  subroutine allocate_EAMType(eam_n_rho,eam_n_r,thisEAMType,stat)
    integer, intent(in)          :: eam_n_rho
    integer, intent(in)          :: eam_n_r
    type (EAMType)               :: thisEAMType
    integer, optional   :: stat
    integer             :: alloc_stat



    if (present(stat)) stat = 0

    allocate(thisEAMType%eam_rvals(eam_n_r),stat=alloc_stat)   
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if
    allocate(thisEAMType%eam_rhovals(eam_n_rho),stat=alloc_stat)   
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if
    allocate(thisEAMType%eam_F_rho(eam_n_rho),stat=alloc_stat)   
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if
    allocate(thisEAMType%eam_Z_r(eam_n_r),stat=alloc_stat)        
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if
    allocate(thisEAMType%eam_rho_r(eam_n_r),stat=alloc_stat)      
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if
    allocate(thisEAMType%eam_phi_r(eam_n_r),stat=alloc_stat)      
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if
    allocate(thisEAMType%eam_F_rho_pp(eam_n_rho),stat=alloc_stat)   
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if
    allocate(thisEAMType%eam_Z_r_pp(eam_n_r),stat=alloc_stat)   
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if
    allocate(thisEAMType%eam_rho_r_pp(eam_n_r),stat=alloc_stat)   
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if
    allocate(thisEAMType%eam_phi_r_pp(eam_n_r),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       if (present(stat)) stat = -1
       return
    end if


  end subroutine allocate_EAMType


  subroutine deallocate_EAMType(thisEAMType)
    type (EAMtype), pointer :: thisEAMType

    ! free Arrays in reverse order of allocation...
    if(associated(thisEAMType%eam_phi_r_pp)) deallocate(thisEAMType%eam_phi_r_pp)      
    if(associated(thisEAMType%eam_rho_r_pp)) deallocate(thisEAMType%eam_rho_r_pp)   
    if(associated(thisEAMType%eam_Z_r_pp)) deallocate(thisEAMType%eam_Z_r_pp)   
    if(associated(thisEAMType%eam_F_rho_pp)) deallocate(thisEAMType%eam_F_rho_pp)   
    if(associated(thisEAMType%eam_phi_r)) deallocate(thisEAMType%eam_phi_r)      
    if(associated(thisEAMType%eam_rho_r)) deallocate(thisEAMType%eam_rho_r)      
    if(associated(thisEAMType%eam_Z_r)) deallocate(thisEAMType%eam_Z_r)   
    if(associated(thisEAMType%eam_F_rho)) deallocate(thisEAMType%eam_F_rho)
    if(associated(thisEAMType%eam_rhovals)) deallocate(thisEAMType%eam_rhovals)
    if(associated(thisEAMType%eam_rvals)) deallocate(thisEAMType%eam_rvals)

  end subroutine deallocate_EAMType

  !! Calculates rho_r
  subroutine calc_eam_prepair_rho(atom1,atom2,d,r,rijsq)
    integer :: atom1,atom2
    real(kind = dp), dimension(3) :: d
    real(kind = dp), intent(inout)               :: r
    real(kind = dp), intent(inout)               :: rijsq
    ! value of electron density rho do to atom i at atom j
    real(kind = dp) :: rho_i_at_j
    ! value of electron density rho do to atom j at atom i
    real(kind = dp) :: rho_j_at_i

    ! we don't use the derivatives, dummy variables
    real( kind = dp) :: drho,d2rho
    integer :: eam_err
    
    integer :: atid1,atid2 ! Global atid    
    integer :: myid_atom1 ! EAM atid
    integer :: myid_atom2 


    ! check to see if we need to be cleaned at the start of a force loop




#ifdef IS_MPI
    Atid1 = Atid_row(Atom1)
    Atid2 = Atid_col(Atom2)
#else
    Atid1 = Atid(Atom1)
    Atid2 = Atid(Atom2)
#endif

    Myid_atom1 = Eamlist%atidtoeamtype(Atid1)
    Myid_atom2 = Eamlist%atidtoeamtype(Atid2)

    if (r.lt.EAMList%EAMParams(myid_atom1)%eam_rcut) then



       call eam_splint(EAMList%EAMParams(myid_atom1)%eam_nr, &
            EAMList%EAMParams(myid_atom1)%eam_rvals, &
            EAMList%EAMParams(myid_atom1)%eam_rho_r, &
            EAMList%EAMParams(myid_atom1)%eam_rho_r_pp, &
            r, rho_i_at_j,drho,d2rho)



#ifdef  IS_MPI
       rho_col(atom2) = rho_col(atom2) + rho_i_at_j
#else
       rho(atom2) = rho(atom2) + rho_i_at_j
#endif
             ! write(*,*) atom1,atom2,r,rho_i_at_j
    endif

    if (r.lt.EAMList%EAMParams(myid_atom2)%eam_rcut) then
       call eam_splint(EAMList%EAMParams(myid_atom2)%eam_nr, &
            EAMList%EAMParams(myid_atom2)%eam_rvals, &
            EAMList%EAMParams(myid_atom2)%eam_rho_r, &
            EAMList%EAMParams(myid_atom2)%eam_rho_r_pp, &
            r, rho_j_at_i,drho,d2rho)




#ifdef  IS_MPI
       rho_row(atom1) = rho_row(atom1) + rho_j_at_i
#else
       rho(atom1) = rho(atom1) + rho_j_at_i
#endif
    endif






  end subroutine calc_eam_prepair_rho




  !! Calculate the functional F(rho) for all local atoms
  subroutine calc_eam_preforce_Frho(nlocal,pot)
    integer :: nlocal
    real(kind=dp) :: pot
    integer :: i,j
    integer :: atom
    real(kind=dp) :: U,U1,U2
    integer :: atype1
    integer :: me,atid1
    integer :: n_rho_points


    cleanme = .true.
    !! Scatter the electron density from  pre-pair calculation back to local atoms
#ifdef IS_MPI
    call scatter(rho_row,rho,plan_atom_row,eam_err)
    if (eam_err /= 0 ) then
       write(errMsg,*) " Error scattering rho_row into rho"
       call handleError(RoutineName,errMesg)
    endif
    call scatter(rho_col,rho_tmp,plan_atom_col,eam_err)
    if (eam_err /= 0 ) then
       write(errMsg,*) " Error scattering rho_col into rho"
       call handleError(RoutineName,errMesg)
    endif

    rho(1:nlocal) = rho(1:nlocal) + rho_tmp(1:nlocal)
#endif



    !! Calculate F(rho) and derivative 
    do atom = 1, nlocal
       atid1 = atid(atom)
       me = eamList%atidToEAMtype(atid1)
       n_rho_points = EAMList%EAMParams(me)%eam_nrho
       !  Check to see that the density is not greater than the larges rho we have calculated
       if (rho(atom) < EAMList%EAMParams(me)%eam_rhovals(n_rho_points)) then
          call eam_splint(n_rho_points, &
               EAMList%EAMParams(me)%eam_rhovals, &
               EAMList%EAMParams(me)%eam_f_rho, &
               EAMList%EAMParams(me)%eam_f_rho_pp, &
               rho(atom), & ! Actual Rho
               u, u1, u2)
       else 
          ! Calculate F(rho with the largest available rho value
          call eam_splint(n_rho_points, &
               EAMList%EAMParams(me)%eam_rhovals, &
               EAMList%EAMParams(me)%eam_f_rho, &
               EAMList%EAMParams(me)%eam_f_rho_pp, &
               EAMList%EAMParams(me)%eam_rhovals(n_rho_points), & ! Largest rho
               u,u1,u2)
       end if


       frho(atom) = u
       dfrhodrho(atom) = u1
       d2frhodrhodrho(atom) = u2
       pot = pot + u

    enddo



#ifdef IS_MPI
    !! communicate f(rho) and derivatives back into row and column arrays
    call gather(frho,frho_row,plan_atom_row, eam_err)
    if (eam_err /=  0) then
       call handleError("cal_eam_forces()","MPI gather frho_row failure")
    endif
    call gather(dfrhodrho,dfrhodrho_row,plan_atom_row, eam_err)
    if (eam_err /=  0) then
       call handleError("cal_eam_forces()","MPI gather dfrhodrho_row failure")
    endif
    call gather(frho,frho_col,plan_atom_col, eam_err)
    if (eam_err /=  0) then
       call handleError("cal_eam_forces()","MPI gather frho_col failure")
    endif
    call gather(dfrhodrho,dfrhodrho_col,plan_atom_col, eam_err)
    if (eam_err /=  0) then
       call handleError("cal_eam_forces()","MPI gather dfrhodrho_col failure")
    endif





    if (nmflag) then
       call gather(d2frhodrhodrho,d2frhodrhodrho_row,plan_atom_row)
       call gather(d2frhodrhodrho,d2frhodrhodrho_col,plan_atom_col)
    endif
#endif


  end subroutine calc_eam_preforce_Frho




  !! Does EAM pairwise Force calculation.  
  subroutine do_eam_pair(atom1, atom2, d, rij, r2, sw, vpair, fpair, &
       pot, f, do_pot)
    !Arguments    
    integer, intent(in) ::  atom1, atom2
    real( kind = dp ), intent(in) :: rij, r2
    real( kind = dp ) :: pot, sw, vpair
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), intent(in), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair

    logical, intent(in) :: do_pot

    real( kind = dp ) :: drdx,drdy,drdz
    real( kind = dp ) :: d2
    real( kind = dp ) :: phab,pha,dvpdr,d2vpdrdr
    real( kind = dp ) :: rha,drha,d2rha, dpha
    real( kind = dp ) :: rhb,drhb,d2rhb, dphb
    real( kind = dp ) :: dudr
    real( kind = dp ) :: rci,rcj
    real( kind = dp ) :: drhoidr,drhojdr
    real( kind = dp ) :: d2rhoidrdr
    real( kind = dp ) :: d2rhojdrdr
    real( kind = dp ) :: Fx,Fy,Fz
    real( kind = dp ) :: r,d2pha,phb,d2phb

    integer :: id1,id2
    integer  :: mytype_atom1
    integer  :: mytype_atom2
    integer  :: atid1,atid2
    !Local Variables

    ! write(*,*) "Frho: ", Frho(atom1)

    phab = 0.0E0_DP
    dvpdr = 0.0E0_DP
    d2vpdrdr = 0.0E0_DP

    if (rij .lt. EAM_rcut) then

#ifdef IS_MPI
       atid1 = atid_row(atom1)
       atid2 = atid_col(atom2)
#else
       atid1 = atid(atom1)
       atid2 = atid(atom2)
#endif

       mytype_atom1 = EAMList%atidToEAMType(atid1)
       mytype_atom2 = EAMList%atidTOEAMType(atid2)


       ! get cutoff for atom 1
       rci = EAMList%EAMParams(mytype_atom1)%eam_rcut
       ! get type specific cutoff for atom 2
       rcj = EAMList%EAMParams(mytype_atom2)%eam_rcut

       drdx = d(1)/rij
       drdy = d(2)/rij
       drdz = d(3)/rij

       if (rij.lt.rci) then
          call eam_splint(EAMList%EAMParams(mytype_atom1)%eam_nr, &
               EAMList%EAMParams(mytype_atom1)%eam_rvals, &
               EAMList%EAMParams(mytype_atom1)%eam_rho_r, &
               EAMList%EAMParams(mytype_atom1)%eam_rho_r_pp, &
               rij, rha,drha,d2rha)
          !! Calculate Phi(r) for atom1.
          call eam_splint(EAMList%EAMParams(mytype_atom1)%eam_nr, &
               EAMList%EAMParams(mytype_atom1)%eam_rvals, &
               EAMList%EAMParams(mytype_atom1)%eam_phi_r, &
               EAMList%EAMParams(mytype_atom1)%eam_phi_r_pp, &
               rij, pha,dpha,d2pha)
       endif

       if (rij.lt.rcj) then      
          ! Calculate rho,drho and d2rho for atom1
          call eam_splint(EAMList%EAMParams(mytype_atom2)%eam_nr, &
               EAMList%EAMParams(mytype_atom2)%eam_rvals, &
               EAMList%EAMParams(mytype_atom2)%eam_rho_r, &
               EAMList%EAMParams(mytype_atom2)%eam_rho_r_pp, &
               rij, rhb,drhb,d2rhb)

          !! Calculate Phi(r) for atom2.
          call eam_splint(EAMList%EAMParams(mytype_atom2)%eam_nr, &
               EAMList%EAMParams(mytype_atom2)%eam_rvals, &
               EAMList%EAMParams(mytype_atom2)%eam_phi_r, &
               EAMList%EAMParams(mytype_atom2)%eam_phi_r_pp, &
               rij, phb,dphb,d2phb)
       endif

       if (rij.lt.rci) then 
          phab = phab + 0.5E0_DP*(rhb/rha)*pha
          dvpdr = dvpdr + 0.5E0_DP*((rhb/rha)*dpha + &
               pha*((drhb/rha) - (rhb*drha/rha/rha)))
          d2vpdrdr = d2vpdrdr + 0.5E0_DP*((rhb/rha)*d2pha + &
               2.0E0_DP*dpha*((drhb/rha) - (rhb*drha/rha/rha)) + &
               pha*((d2rhb/rha) - 2.0E0_DP*(drhb*drha/rha/rha) + &
               (2.0E0_DP*rhb*drha*drha/rha/rha/rha) - (rhb*d2rha/rha/rha)))
       endif

       if (rij.lt.rcj) then
          phab = phab + 0.5E0_DP*(rha/rhb)*phb
          dvpdr = dvpdr + 0.5E0_DP*((rha/rhb)*dphb + &
               phb*((drha/rhb) - (rha*drhb/rhb/rhb)))
          d2vpdrdr = d2vpdrdr + 0.5E0_DP*((rha/rhb)*d2phb + &
               2.0E0_DP*dphb*((drha/rhb) - (rha*drhb/rhb/rhb)) + &
               phb*((d2rha/rhb) - 2.0E0_DP*(drha*drhb/rhb/rhb) + &
               (2.0E0_DP*rha*drhb*drhb/rhb/rhb/rhb) - (rha*d2rhb/rhb/rhb)))
       endif

       drhoidr = drha
       drhojdr = drhb

       d2rhoidrdr = d2rha
       d2rhojdrdr = d2rhb


#ifdef IS_MPI
       dudr = drhojdr*dfrhodrho_row(atom1)+drhoidr*dfrhodrho_col(atom2) &
            + dvpdr

#else
       dudr = drhojdr*dfrhodrho(atom1)+drhoidr*dfrhodrho(atom2) &
            + dvpdr
       ! write(*,*) "Atom1,Atom2, dfrhodrho(atom1) dfrhodrho(atom2): ", atom1,atom2,dfrhodrho(atom1),dfrhodrho(atom2)
#endif

       fx = dudr * drdx
       fy = dudr * drdy
       fz = dudr * drdz


#ifdef IS_MPI
       if (do_pot) then
          pot_Row(METALLIC_POT,atom1) = pot_Row(METALLIC_POT,atom1) + phab*0.5
          pot_Col(METALLIC_POT,atom2) = pot_Col(METALLIC_POT,atom2) + phab*0.5
       end if

       f_Row(1,atom1) = f_Row(1,atom1) + fx
       f_Row(2,atom1) = f_Row(2,atom1) + fy
       f_Row(3,atom1) = f_Row(3,atom1) + fz

       f_Col(1,atom2) = f_Col(1,atom2) - fx
       f_Col(2,atom2) = f_Col(2,atom2) - fy
       f_Col(3,atom2) = f_Col(3,atom2) - fz
#else

       if(do_pot) then
          pot = pot + phab
       end if

       f(1,atom1) = f(1,atom1) + fx
       f(2,atom1) = f(2,atom1) + fy
       f(3,atom1) = f(3,atom1) + fz

       f(1,atom2) = f(1,atom2) - fx
       f(2,atom2) = f(2,atom2) - fy
       f(3,atom2) = f(3,atom2) - fz
#endif

       vpair = vpair + phab
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

       if (nmflag) then

          drhoidr = drha
          drhojdr = drhb
          d2rhoidrdr = d2rha
          d2rhojdrdr = d2rhb

#ifdef IS_MPI
          d2 = d2vpdrdr + &
               d2rhoidrdr*dfrhodrho_col(atom2) + &
               d2rhojdrdr*dfrhodrho_row(atom1) + &
               drhoidr*drhoidr*d2frhodrhodrho_col(atom2) + &
               drhojdr*drhojdr*d2frhodrhodrho_row(atom1)

#else

          d2 = d2vpdrdr + &
               d2rhoidrdr*dfrhodrho(atom2) + &
               d2rhojdrdr*dfrhodrho(atom1) + &
               drhoidr*drhoidr*d2frhodrhodrho(atom2) + &
               drhojdr*drhojdr*d2frhodrhodrho(atom1)
#endif
       end if

    endif
  end subroutine do_eam_pair


  subroutine eam_splint(nx, xa, ya, yppa, x, y, dy, d2y)

    integer :: atype, nx, j
    real( kind = DP ), dimension(:) :: xa
    real( kind = DP ), dimension(:) :: ya
    real( kind = DP ), dimension(:) :: yppa
    real( kind = DP ) :: x, y
    real( kind = DP ) :: dy, d2y
    real( kind = DP ) :: del, h, a, b, c, d
    integer :: pp_arraySize


    ! this spline code assumes that the x points are equally spaced
    ! do not attempt to use this code if they are not.


    ! find the closest point with a value below our own:
    j = FLOOR(real((nx-1),kind=dp) * (x - xa(1)) / (xa(nx) - xa(1))) + 1

    ! check to make sure we're inside the spline range:
    if ((j.gt.nx).or.(j.lt.1)) then
       write(errMSG,*) "EAM_splint: x is outside bounds of spline: ",x,j
       call handleError(routineName,errMSG)
    endif
    ! check to make sure we haven't screwed up the calculation of j:
    if ((x.lt.xa(j)).or.(x.gt.xa(j+1))) then
       if (j.ne.nx) then
          write(errMSG,*) "EAM_splint:",x," x is outside bounding range"
          call handleError(routineName,errMSG)
       endif
    endif

    del = xa(j+1) - x
    h = xa(j+1) - xa(j)

    a = del / h
    b = 1.0E0_DP - a
    c = a*(a*a - 1.0E0_DP)*h*h/6.0E0_DP
    d = b*(b*b - 1.0E0_DP)*h*h/6.0E0_DP

    y = a*ya(j) + b*ya(j+1) + c*yppa(j) + d*yppa(j+1)

    dy = (ya(j+1)-ya(j))/h &
         - (3.0E0_DP*a*a - 1.0E0_DP)*h*yppa(j)/6.0E0_DP &
         + (3.0E0_DP*b*b - 1.0E0_DP)*h*yppa(j+1)/6.0E0_DP


    d2y = a*yppa(j) + b*yppa(j+1)


  end subroutine eam_splint


  subroutine eam_spline(nx, xa, ya, yppa, yp1, ypn, boundary)


    ! yp1 and ypn are the first derivatives of y at the two endpoints
    ! if boundary is 'L' the lower derivative is used
    ! if boundary is 'U' the upper derivative is used
    ! if boundary is 'B' then both derivatives are used
    ! if boundary is anything else, then both derivatives are assumed to be 0

    integer :: nx, i, k, max_array_size

    real( kind = DP ), dimension(:)        :: xa
    real( kind = DP ), dimension(:)        :: ya
    real( kind = DP ), dimension(:)        :: yppa
    real( kind = DP ), dimension(size(xa)) :: u
    real( kind = DP ) :: yp1,ypn,un,qn,sig,p
    character(len=*) :: boundary

    ! make sure the sizes match
    if ((nx /= size(xa)) .or. (nx /= size(ya))) then
       call handleWarning("EAM_SPLINE","Array size mismatch")
    end if

    if ((boundary.eq.'l').or.(boundary.eq.'L').or. &
         (boundary.eq.'b').or.(boundary.eq.'B')) then
       yppa(1) = -0.5E0_DP
       u(1) = (3.0E0_DP/(xa(2)-xa(1)))*((ya(2)-&
            ya(1))/(xa(2)-xa(1))-yp1)
    else
       yppa(1) = 0.0E0_DP
       u(1)  = 0.0E0_DP
    endif

    do i = 2, nx - 1
       sig = (xa(i) - xa(i-1)) / (xa(i+1) - xa(i-1))
       p = sig * yppa(i-1) + 2.0E0_DP
       yppa(i) = (sig - 1.0E0_DP) / p
       u(i) = (6.0E0_DP*((ya(i+1)-ya(i))/(xa(i+1)-xa(i)) - &
            (ya(i)-ya(i-1))/(xa(i)-xa(i-1)))/ &
            (xa(i+1)-xa(i-1)) - sig * u(i-1))/p
    enddo

    if ((boundary.eq.'u').or.(boundary.eq.'U').or. &
         (boundary.eq.'b').or.(boundary.eq.'B')) then
       qn = 0.5E0_DP
       un = (3.0E0_DP/(xa(nx)-xa(nx-1)))* &
            (ypn-(ya(nx)-ya(nx-1))/(xa(nx)-xa(nx-1)))
    else
       qn = 0.0E0_DP
       un = 0.0E0_DP
    endif

    yppa(nx)=(un-qn*u(nx-1))/(qn*yppa(nx-1)+1.0E0_DP)

    do k = nx-1, 1, -1
       yppa(k)=yppa(k)*yppa(k+1)+u(k)
    enddo

  end subroutine eam_spline

end module eam
