module dipole_dipole
  
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

  real(kind=dp), parameter :: pre = 14.38362_dp
  logical, save :: haveMomentMap = .false.

  public::do_dipole_pair

  type :: MomentList
     real(kind=DP) :: dipole_moment = 0.0_DP
  end type MomentList

  type(MomentList), dimension(:),allocatable :: MomentMap

contains

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
  
  subroutine do_dipole_pair(atom1, atom2, d, rij, r2, sw, vpair, fpair, &
       pot, u_l, f, t, do_pot)
    
    logical :: do_pot

    integer atom1, atom2, me1, me2, id1, id2
    integer :: localError
    real(kind=dp) :: rij, mu1, mu2
    real(kind=dp) :: dfact1, dfact2, dip2, r2, r3, r5
    real(kind=dp) :: dudx, dudy, dudz, dudu1x, dudu1y, dudu1z
    real(kind=dp) :: dudu2x, dudu2y, dudu2z, rdotu1, rdotu2, u1dotu2
    real(kind=dp) :: sw, vpair, vterm

    real( kind = dp ) :: pot
    real( kind = dp ), dimension(3) :: d, fpair
    real( kind = dp ), dimension(3,nLocal) :: u_l
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), dimension(3,nLocal) :: t
    
    real (kind = dp), dimension(3) :: ul1
    real (kind = dp), dimension(3) :: ul2

    if (.not.haveMomentMap) then
       localError = 0
       call createMomentMap(localError)
       if ( localError .ne. 0 ) then
          call handleError("dipole-dipole", "MomentMap creation failed!")
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
    
    r3 = r2*rij
    r5 = r3*r2
    
    rdotu1 = d(1)*ul1(1) + d(2)*ul1(2) + d(3)*ul1(3)
    rdotu2 = d(1)*ul2(1) + d(2)*ul2(2) + d(3)*ul2(3)
    u1dotu2 = ul1(1)*ul2(1) + ul1(2)*ul2(2) + ul1(3)*ul2(3)
    
    dip2 = pre * mu1 * mu2
    dfact1 = 3.0d0*dip2 / r2
    dfact2 = 3.0d0*dip2 / r5
    
    vterm = dip2*((u1dotu2/r3) - 3.0d0*(rdotu1*rdotu2/r5))
    
    vpair = vpair + vterm
    
    if (do_pot) then
#ifdef IS_MPI 
       pot_row(atom1) = pot_row(atom1) + 0.5d0*vterm*sw
       pot_col(atom2) = pot_col(atom2) + 0.5d0*vterm*sw
#else
       pot = pot + vterm*sw
#endif
    endif
    
    dudx = (-dfact1 * d(1) * ((u1dotu2/r3) - &
         (5.0d0*(rdotu1*rdotu2)/r5)) -  &
            dfact2*(ul1(1)*rdotu2 + ul2(1)*rdotu1))*sw
    
    dudy = (-dfact1 * d(2) * ((u1dotu2/r3) - &
         (5.0d0*(rdotu1*rdotu2)/r5)) -  &
            dfact2*(ul1(2)*rdotu2 + ul2(2)*rdotu1))*sw
    
    dudz = (-dfact1 * d(3) * ((u1dotu2/r3) - &
         (5.0d0*(rdotu1*rdotu2)/r5)) -  &
         dfact2*(ul1(3)*rdotu2 + ul2(3)*rdotu1))*sw
    
    dudu1x = (dip2*((ul2(1)/r3) - (3.0d0*d(1)*rdotu2/r5)))*sw
    dudu1y = (dip2*((ul2(2)/r3) - (3.0d0*d(2)*rdotu2/r5)))*sw
    dudu1z = (dip2*((ul2(3)/r3) - (3.0d0*d(3)*rdotu2/r5)))*sw
    
    dudu2x = (dip2*((ul1(1)/r3) - (3.0d0*d(1)*rdotu1/r5)))*sw
    dudu2y = (dip2*((ul1(2)/r3) - (3.0d0*d(2)*rdotu1/r5)))*sw
    dudu2z = (dip2*((ul1(3)/r3) - (3.0d0*d(3)*rdotu1/r5)))*sw
    
    
#ifdef IS_MPI
    f_Row(1,atom1) = f_Row(1,atom1) + dudx
    f_Row(2,atom1) = f_Row(2,atom1) + dudy
    f_Row(3,atom1) = f_Row(3,atom1) + dudz
    
    f_Col(1,atom2) = f_Col(1,atom2) - dudx
    f_Col(2,atom2) = f_Col(2,atom2) - dudy
    f_Col(3,atom2) = f_Col(3,atom2) - dudz
    
    t_Row(1,atom1) = t_Row(1,atom1) - ul1(2)*dudu1z + ul1(3)*dudu1y
    t_Row(2,atom1) = t_Row(2,atom1) - ul1(3)*dudu1x + ul1(1)*dudu1z
    t_Row(3,atom1) = t_Row(3,atom1) - ul1(1)*dudu1y + ul1(2)*dudu1x
    
    t_Col(1,atom2) = t_Col(1,atom2) - ul2(2)*dudu2z + ul2(3)*dudu2y
    t_Col(2,atom2) = t_Col(2,atom2) - ul2(3)*dudu2x + ul2(1)*dudu2z
    t_Col(3,atom2) = t_Col(3,atom2) - ul2(1)*dudu2y + ul2(2)*dudu2x
#else
    f(1,atom1) = f(1,atom1) + dudx
    f(2,atom1) = f(2,atom1) + dudy
    f(3,atom1) = f(3,atom1) + dudz
    
    f(1,atom2) = f(1,atom2) - dudx
    f(2,atom2) = f(2,atom2) - dudy
    f(3,atom2) = f(3,atom2) - dudz
    
    t(1,atom1) = t(1,atom1) - ul1(2)*dudu1z + ul1(3)*dudu1y
    t(2,atom1) = t(2,atom1) - ul1(3)*dudu1x + ul1(1)*dudu1z
    t(3,atom1) = t(3,atom1) - ul1(1)*dudu1y + ul1(2)*dudu1x
    
    t(1,atom2) = t(1,atom2) - ul2(2)*dudu2z + ul2(3)*dudu2y
    t(2,atom2) = t(2,atom2) - ul2(3)*dudu2x + ul2(1)*dudu2z
    t(3,atom2) = t(3,atom2) - ul2(1)*dudu2y + ul2(2)*dudu2x
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
  end subroutine do_dipole_pair
  
end module dipole_dipole
