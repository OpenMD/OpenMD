! Fortran interface to C entry plug.

module force_globals
  use definitions
#ifdef IS_MPI
  use mpiSimulation
#endif

  implicit none
  PRIVATE

  logical, save :: force_globals_initialized = .false.

#ifdef IS_MPI
  real( kind = dp ), allocatable, dimension(:,:), public :: q_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: q_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: q_group_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: q_group_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: u_l_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: u_l_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: A_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: A_Col
  
  real( kind = dp ), allocatable, dimension(:), public :: pot_Row
  real( kind = dp ), allocatable, dimension(:), public :: pot_Col
  real( kind = dp ), allocatable, dimension(:), public :: pot_Temp
  real( kind = dp ), allocatable, dimension(:,:), public :: f_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: f_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: f_Temp
  real( kind = dp ), allocatable, dimension(:,:), public :: t_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: t_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: t_Temp
  real( kind = dp ), allocatable, dimension(:,:), public :: rf_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: rf_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: rf_Temp

  integer, allocatable, dimension(:), public :: atid_Row
  integer, allocatable, dimension(:), public :: atid_Col
#endif

  integer, allocatable, dimension(:), public :: atid

  real( kind = dp ), allocatable, dimension(:,:), public :: rf
  real(kind = dp), dimension(9), public :: tau_Temp = 0.0_dp
  real(kind = dp), public :: virial_Temp = 0.0_dp
  
  public :: InitializeForceGlobals
  
contains
  
  subroutine InitializeForceGlobals(nlocal, thisStat)
    integer, intent(out) :: thisStat
    integer :: nAtomsInRow, nAtomsInCol
    integer :: nGroupsInRow, nGroupsInCol
    integer :: nlocal
    integer :: ndim = 3
    integer :: alloc_stat
    
    thisStat = 0
    
#ifdef IS_MPI
    nAtomsInRow = getNatomsInRow(plan_atom_row)
    nAtomsInCol = getNatomsInCol(plan_atom_col)
    nGroupsInRow = getNgroupsInRow(plan_group_row)
    nGroupsInCol = getNgroupsInCol(plan_group_col)
    
#endif
    
    call FreeForceGlobals()
    
#ifdef IS_MPI

    allocate(q_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
        
    allocate(q_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(q_group_Row(ndim,nGroupsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
        
    allocate(q_group_Col(ndim,nGroupsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
       
    allocate(u_l_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
     
    allocate(u_l_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
    
    allocate(A_row(9,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
        
    allocate(A_Col(9,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
    
    allocate(pot_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
    
    allocate(pot_Col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(pot_Temp(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
    
    allocate(f_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
    
    allocate(f_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
    
    allocate(f_Temp(ndim,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
    
    allocate(t_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif
    
    allocate(t_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(t_temp(ndim,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(atid(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif


    allocate(atid_Row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(atid_Col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(rf_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(rf_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(rf_Temp(ndim,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

#else

    allocate(atid(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    end if

#endif
    
    allocate(rf(ndim,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    force_globals_initialized = .true.

  end subroutine InitializeForceGlobals
  
  subroutine FreeForceGlobals()
    
    !We free in the opposite order in which we allocate in.
    
    if (allocated(rf))         deallocate(rf)
#ifdef IS_MPI
    if (allocated(rf_Temp))    deallocate(rf_Temp)
    if (allocated(rf_Col))     deallocate(rf_Col)
    if (allocated(rf_Row))     deallocate(rf_Row)    
    if (allocated(atid_Col))   deallocate(atid_Col)
    if (allocated(atid_Row))   deallocate(atid_Row)
    if (allocated(atid))       deallocate(atid)
    if (allocated(t_Temp))     deallocate(t_Temp)
    if (allocated(t_Col))      deallocate(t_Col)
    if (allocated(t_Row))      deallocate(t_Row)
    if (allocated(f_Temp))     deallocate(f_Temp)
    if (allocated(f_Col))      deallocate(f_Col)
    if (allocated(f_Row))      deallocate(f_Row)
    if (allocated(pot_Temp))   deallocate(pot_Temp)
    if (allocated(pot_Col))    deallocate(pot_Col)
    if (allocated(pot_Row))    deallocate(pot_Row)
    if (allocated(A_Col))      deallocate(A_Col)
    if (allocated(A_Row))      deallocate(A_Row)
    if (allocated(u_l_Col))    deallocate(u_l_Col)
    if (allocated(u_l_Row))    deallocate(u_l_Row)
    if (allocated(q_group_Col)) deallocate(q_group_Col)
    if (allocated(q_group_Row)) deallocate(q_group_Row)    
    if (allocated(q_Col))      deallocate(q_Col)
    if (allocated(q_Row))      deallocate(q_Row)    
#else    
    if (allocated(atid))       deallocate(atid)    
#endif
        
  end subroutine FreeForceGlobals
    
end module force_globals
