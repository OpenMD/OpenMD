module reaction_field
  use force_globals
  use definitions
  use atype_module
  use vector_class
  use simulation
  use status
#ifdef IS_MPI
  use mpiSimulation
#endif
  implicit none

  PRIVATE
  
  real(kind=dp), save :: rrf = 1.0_dp
  real(kind=dp), save :: rt
  real(kind=dp), save :: dielect = 1.0_dp
  real(kind=dp), save :: rrfsq = 1.0_dp
  real(kind=dp), save :: pre
  logical, save :: haveCutoffs = .false.
  logical, save :: haveMomentMap = .false.
  logical, save :: haveDielectric = .false.

  type :: MomentList
     real(kind=DP) :: dipole_moment = 0.0_DP
  end type MomentList

  type(MomentList), dimension(:),allocatable :: MomentMap

  PUBLIC::initialize_rf
  PUBLIC::setCutoffsRF
  PUBLIC::accumulate_rf
  PUBLIC::accumulate_self_rf
  PUBLIC::reaction_field_final
  PUBLIC::rf_correct_forces

contains
  
  subroutine initialize_rf(this_dielect)
    real(kind=dp), intent(in) :: this_dielect
    
    dielect = this_dielect

    pre = 14.38362d0*2.0d0*(dielect-1.0d0)/((2.0d0*dielect+1.0d0)*rrfsq*rrf) 
    
    haveDielectric = .true.

    return
  end subroutine initialize_rf

  subroutine setCutoffsRF( this_rrf, this_rt )

    real(kind=dp), intent(in) :: this_rrf, this_rt

    rrf = this_rrf
    rt  = this_rt

    rrfsq = rrf * rrf
    pre = 14.38362d0*2.0d0*(dielect-1.0d0)/((2.0d0*dielect+1.0d0)*rrfsq*rrf)
    
    haveCutoffs = .true.

  end subroutine setCutoffsRF

  subroutine createMomentMap(status)
    integer :: nAtypes
    integer :: status
    integer :: i
    real (kind=DP) :: thisDP
    logical :: thisProperty

    status = 0

    nAtypes = getSize(atypes)
    
    if (nAtypes == 0) then
       status = -1
       return
    end if
    
    if (.not. allocated(MomentMap)) then
       allocate(MomentMap(nAtypes))
    endif

    do i = 1, nAtypes

       call getElementProperty(atypes, i, "is_DP", thisProperty)

       if (thisProperty) then
          call getElementProperty(atypes, i, "dipole_moment", thisDP)
          MomentMap(i)%dipole_moment = thisDP
       endif
       
    end do
    
    haveMomentMap = .true.
    
  end subroutine createMomentMap  

  subroutine accumulate_rf(atom1, atom2, rij, u_l, taper)

    integer, intent(in) :: atom1, atom2
    real (kind = dp), intent(in) :: rij
    real (kind = dp), dimension(3,nLocal) :: u_l    

    integer :: me1, me2
    real (kind = dp), intent(in) :: taper
    real (kind = dp):: mu1, mu2
    real (kind = dp), dimension(3) :: ul1
    real (kind = dp), dimension(3) :: ul2   

    integer :: localError

    if ((.not.haveDielectric).or.(.not.haveCutoffs)) then
       write(default_error,*) 'Reaction field not initialized!'
       return
    endif

    if (.not.haveMomentMap) then
       localError = 0
       call createMomentMap(localError)
       if ( localError .ne. 0 ) then
          call handleError("reaction-field", "MomentMap creation failed!")
          return
       end if
    endif
    
       
#ifdef IS_MPI
    me1 = atid_Row(atom1)
    ul1(1) = u_l_Row(1,atom1)
    ul1(2) = u_l_Row(2,atom1)
    ul1(3) = u_l_Row(3,atom1)
    
    me2 = atid_Col(atom2)
    ul2(1) = u_l_Col(1,atom2)
    ul2(2) = u_l_Col(2,atom2)
    ul2(3) = u_l_Col(3,atom2)
#else
    me1 = atid(atom1)
    ul1(1) = u_l(1,atom1)
    ul1(2) = u_l(2,atom1)
    ul1(3) = u_l(3,atom1)
    
    me2 = atid(atom2)
    ul2(1) = u_l(1,atom2)
    ul2(2) = u_l(2,atom2)
    ul2(3) = u_l(3,atom2)
#endif
    
    mu1 = MomentMap(me1)%dipole_moment
    mu2 = MomentMap(me2)%dipole_moment
    
#ifdef IS_MPI
    rf_Row(1,atom1) = rf_Row(1,atom1) + ul2(1)*mu2*taper
    rf_Row(2,atom1) = rf_Row(2,atom1) + ul2(2)*mu2*taper
    rf_Row(3,atom1) = rf_Row(3,atom1) + ul2(3)*mu2*taper
    
    rf_Col(1,atom2) = rf_Col(1,atom2) + ul1(1)*mu1*taper
    rf_Col(2,atom2) = rf_Col(2,atom2) + ul1(2)*mu1*taper
    rf_Col(3,atom2) = rf_Col(3,atom2) + ul1(3)*mu1*taper
#else
    rf(1,atom1) = rf(1,atom1) + ul2(1)*mu2*taper
    rf(2,atom1) = rf(2,atom1) + ul2(2)*mu2*taper
    rf(3,atom1) = rf(3,atom1) + ul2(3)*mu2*taper
    
    rf(1,atom2) = rf(1,atom2) + ul1(1)*mu1*taper
    rf(2,atom2) = rf(2,atom2) + ul1(2)*mu1*taper
    rf(3,atom2) = rf(3,atom2) + ul1(3)*mu1*taper     
#endif
    
    
    return  
  end subroutine accumulate_rf

  subroutine accumulate_self_rf(atom1, mu1, u_l)
    
    integer, intent(in) :: atom1
    real(kind=dp), intent(in) :: mu1
    real(kind=dp), dimension(3,nLocal) :: u_l
    
    !! should work for both MPI and non-MPI version since this is not pairwise.
    rf(1,atom1) = rf(1,atom1) + u_l(1,atom1)*mu1
    rf(2,atom1) = rf(2,atom1) + u_l(2,atom1)*mu1
    rf(3,atom1) = rf(3,atom1) + u_l(3,atom1)*mu1
        
    return
  end subroutine accumulate_self_rf
  
  subroutine reaction_field_final(a1, mu1, u_l, rfpot, t, do_pot)
            
    integer, intent(in) :: a1
    real (kind=dp), intent(in) :: mu1
    real (kind=dp), intent(inout) :: rfpot
    logical, intent(in) :: do_pot
    real (kind = dp), dimension(3,nLocal) :: u_l    
    real (kind = dp), dimension(3,nLocal) :: t

    integer :: localError

    if ((.not.haveDielectric).or.(.not.haveCutoffs)) then
       write(default_error,*) 'Reaction field not initialized!'
       return
    endif

    if (.not.haveMomentMap) then
       localError = 0
       call createMomentMap(localError)
       if ( localError .ne. 0 ) then
          call handleError("reaction-field", "MomentMap creation failed!")
          return
       end if
    endif

    ! compute torques on dipoles:
    ! pre converts from mu in units of debye to kcal/mol
    
    ! The torque contribution is dipole cross reaction_field   

    t(1,a1) = t(1,a1) + pre*mu1*(u_l(2,a1)*rf(3,a1) - u_l(3,a1)*rf(2,a1))
    t(2,a1) = t(2,a1) + pre*mu1*(u_l(3,a1)*rf(1,a1) - u_l(1,a1)*rf(3,a1))
    t(3,a1) = t(3,a1) + pre*mu1*(u_l(1,a1)*rf(2,a1) - u_l(2,a1)*rf(1,a1))
    
    ! the potential contribution is -1/2 dipole dot reaction_field
    
    if (do_pot) then
       rfpot = rfpot - 0.5d0 * pre * mu1 * &
            (rf(1,a1)*u_l(1,a1) + rf(2,a1)*u_l(2,a1) + rf(3,a1)*u_l(3,a1))
    endif
 
    return
  end subroutine reaction_field_final
  
  subroutine rf_correct_forces(atom1, atom2, d, rij, u_l, taper, f, fpair)
    
    integer, intent(in) :: atom1, atom2
    real(kind=dp), dimension(3), intent(in) :: d
    real(kind=dp), intent(in) :: rij, taper
    real( kind = dp ), dimension(3,nLocal) :: u_l
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), dimension(3), intent(inout) :: fpair
    
    real (kind = dp), dimension(3) :: ul1
    real (kind = dp), dimension(3) :: ul2
    real (kind = dp) :: dtdr
    real (kind = dp) :: dudx, dudy, dudz, u1dotu2
    integer :: me1, me2, id1, id2
    real (kind = dp) :: mu1, mu2
    
    integer :: localError

    if ((.not.haveDielectric).or.(.not.haveCutoffs)) then
       write(default_error,*) 'Reaction field not initialized!'
       return
    endif

    if (.not.haveMomentMap) then
       localError = 0
       call createMomentMap(localError)
       if ( localError .ne. 0 ) then
          call handleError("reaction-field", "MomentMap creation failed!")
          return
       end if
    endif

    if (rij.le.rrf) then
       
       if (rij.lt.rt) then
          dtdr = 0.0d0
       else
 !         write(*,*) 'rf correct in taper region'
          dtdr = 6.0d0*(rij*rij - rij*rt - rij*rrf +rrf*rt)/((rrf-rt)**3)
       endif
       
#ifdef IS_MPI
       me1 = atid_Row(atom1)
       ul1(1) = u_l_Row(1,atom1)
       ul1(2) = u_l_Row(2,atom1)
       ul1(3) = u_l_Row(3,atom1)
       
       me2 = atid_Col(atom2)
       ul2(1) = u_l_Col(1,atom2)
       ul2(2) = u_l_Col(2,atom2)
       ul2(3) = u_l_Col(3,atom2)
#else
       me1 = atid(atom1)
       ul1(1) = u_l(1,atom1)
       ul1(2) = u_l(2,atom1)
       ul1(3) = u_l(3,atom1)
       
       me2 = atid(atom2)
       ul2(1) = u_l(1,atom2)
       ul2(2) = u_l(2,atom2)
       ul2(3) = u_l(3,atom2)
#endif
       
       mu1 = MomentMap(me1)%dipole_moment
       mu2 = MomentMap(me2)%dipole_moment
       
       u1dotu2 = ul1(1)*ul2(1) + ul1(2)*ul2(2) + ul1(3)*ul2(3)
       
       dudx = - pre*mu1*mu2*u1dotu2*dtdr*d(1)/rij
       dudy = - pre*mu1*mu2*u1dotu2*dtdr*d(2)/rij
       dudz = - pre*mu1*mu2*u1dotu2*dtdr*d(3)/rij
       
#ifdef IS_MPI
       f_Row(1,atom1) = f_Row(1,atom1) + dudx
       f_Row(2,atom1) = f_Row(2,atom1) + dudy
       f_Row(3,atom1) = f_Row(3,atom1) + dudz
       
       f_Col(1,atom2) = f_Col(1,atom2) - dudx
       f_Col(2,atom2) = f_Col(2,atom2) - dudy
       f_Col(3,atom2) = f_Col(3,atom2) - dudz
#else
       f(1,atom1) = f(1,atom1) + dudx
       f(2,atom1) = f(2,atom1) + dudy
       f(3,atom1) = f(3,atom1) + dudz
       
       f(1,atom2) = f(1,atom2) - dudx
       f(2,atom2) = f(2,atom2) - dudy
       f(3,atom2) = f(3,atom2) - dudz
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
       
    end if
    return
  end subroutine rf_correct_forces
end module reaction_field
