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
#ifdef IS_MPI
  use mpiSimulation
#endif
  implicit none

  PRIVATE

  real(kind=dp), parameter :: pre11 = 332.0637778_dp
  real(kind=dp), parameter :: pre12 = 69.13291783_dp
  real(kind=dp), parameter :: pre22 = 14.39289874_dp

  public :: newElectrostaticType
  public :: setCharge
  public :: setDipoleMoment
  public :: setSplitDipoleDistance
  public :: setQuadrupoleMoments
  public :: doElectrostaticPair
  public :: getCharge
  public :: getDipoleMoment

  type :: Electrostatic
     integer :: c_ident
     logical :: is_Charge = .false.
     logical :: is_Dipole = .false.
     logical :: is_SplitDipole = .false.
     logical :: is_Quadrupole = .false.
     real(kind=DP) :: charge = 0.0_DP
     real(kind=DP) :: dipole_moment = 0.0_DP
     real(kind=DP) :: split_dipole_distance = 0.0_DP
     real(kind=DP), dimension(3) :: quadrupole_moments = 0.0_DP
  end type Electrostatic

  type(Electrostatic), dimension(:), allocatable :: ElectrostaticMap

contains

  subroutine newElectrostaticType(c_ident, is_Charge, is_Dipole, &
       is_SplitDipole, is_Quadrupole, status)
    
    integer, intent(in) :: c_ident
    logical, intent(in) :: is_Charge
    logical, intent(in) :: is_Dipole
    logical, intent(in) :: is_SplitDipole
    logical, intent(in) :: is_Quadrupole
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
       
       if (.not. allocated(ElectrostaticMap)) then
          allocate(ElectrostaticMap(nAtypes))
       endif
       
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
    
  end subroutine newElectrostaticType

  subroutine setCharge(c_ident, charge, status)
    integer, intent(in) :: c_ident
    real(kind=dp), intent(in) :: charge
    integer, intent(out) :: status
    integer :: myATID

    status = 0
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    if (.not.allocated(ElectrostaticMap)) then
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

    if (.not.allocated(ElectrostaticMap)) then
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

    if (.not.allocated(ElectrostaticMap)) then
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

    if (.not.allocated(ElectrostaticMap)) then
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
    
    if (.not.allocated(ElectrostaticMap)) then
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
    
    if (.not.allocated(ElectrostaticMap)) then
       call handleError("electrostatic", "no ElectrostaticMap was present before first call of getDipoleMoment!")
       return
    end if
    
    if (.not.ElectrostaticMap(atid)%is_Dipole) then
       call handleError("electrostatic", "getDipoleMoment was called for an atom type that isn't a dipole!")
       return
    endif
    
    dm = ElectrostaticMap(atid)%dipole_moment
  end function getDipoleMoment

  subroutine doElectrostaticPair(atom1, atom2, d, rij, r2, sw, &
       vpair, fpair, pot, eFrame, f, t, do_pot)
    
    logical, intent(in) :: do_pot
    
    integer, intent(in) :: atom1, atom2
    integer :: localError

    real(kind=dp), intent(in) :: rij, r2, sw
    real(kind=dp), intent(in), dimension(3) :: d
    real(kind=dp), intent(inout) :: vpair
    real(kind=dp), intent(inout), dimension(3) :: fpair

    real( kind = dp ) :: pot
    real( kind = dp ), dimension(9,nLocal) :: eFrame
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), dimension(3,nLocal) :: t
    
    real (kind = dp), dimension(3) :: ul_i
    real (kind = dp), dimension(3) :: ul_j

    logical :: i_is_Charge, i_is_Dipole, i_is_SplitDipole, i_is_Quadrupole
    logical :: j_is_Charge, j_is_Dipole, j_is_SplitDipole, j_is_Quadrupole
    integer :: me1, me2, id1, id2
    real (kind=dp) :: q_i, q_j, mu_i, mu_j, d_i, d_j
    real (kind=dp) :: ct_i, ct_j, ct_ij, a1
    real (kind=dp) :: riji, ri2, ri3, ri4
    real (kind=dp) :: pref, vterm, epot, dudr    
    real (kind=dp) :: dudx, dudy, dudz
    real (kind=dp) :: drdxj, drdyj, drdzj
    real (kind=dp) :: duduix, duduiy, duduiz, dudujx, dudujy, dudujz


    if (.not.allocated(ElectrostaticMap)) then
       call handleError("electrostatic", "no ElectrostaticMap was present before first call of do_electrostatic_pair!")
       return
    end if

#ifdef IS_MPI
    me1 = atid_Row(atom1)
    me2 = atid_Col(atom2)
#else
    me1 = atid(atom1)
    me2 = atid(atom2)
#endif

    !! some variables we'll need independent of electrostatic type:

    riji = 1.0d0 / rij

    !! these are also useful as the unit vector of \vec{r} 
    !! \hat{r} = \vec{r} / r =   {(x_j-x_i) / r, (y_j-y_i)/r, (z_j-z_i)/r}

    drdxj = d(1) * riji
    drdyj = d(2) * riji
    drdzj = d(3) * riji

    !! logicals

    i_is_Charge = ElectrostaticMap(me1)%is_Charge
    i_is_Dipole = ElectrostaticMap(me1)%is_Dipole
    i_is_SplitDipole = ElectrostaticMap(me1)%is_SplitDipole
    i_is_Quadrupole = ElectrostaticMap(me1)%is_Quadrupole

    j_is_Charge = ElectrostaticMap(me2)%is_Charge
    j_is_Dipole = ElectrostaticMap(me2)%is_Dipole
    j_is_SplitDipole = ElectrostaticMap(me2)%is_SplitDipole
    j_is_Quadrupole = ElectrostaticMap(me2)%is_Quadrupole

    if (i_is_Charge) then
       q_i = ElectrostaticMap(me1)%charge      
    endif
    
    if (i_is_Dipole) then
       mu_i = ElectrostaticMap(me1)%dipole_moment
#ifdef IS_MPI
       ul_i(1) = eFrame_Row(3,atom1)
       ul_i(2) = eFrame_Row(6,atom1)
       ul_i(3) = eFrame_Row(9,atom1)
#else
       ul_i(1) = eFrame(3,atom1)
       ul_i(2) = eFrame(6,atom1)
       ul_i(3) = eFrame(9,atom1)
#endif
       ct_i = ul_i(1)*drdxj + ul_i(2)*drdyj + ul_i(3)*drdzj

       if (i_is_SplitDipole) then
          d_i = ElectrostaticMap(me1)%split_dipole_distance
       endif
       
    endif

    if (j_is_Charge) then
       q_j = ElectrostaticMap(me2)%charge      
    endif
    
    if (j_is_Dipole) then
       mu_j = ElectrostaticMap(me2)%dipole_moment
#ifdef IS_MPI
       ul_j(1) = eFrame_Col(3,atom2)
       ul_j(2) = eFrame_Col(6,atom2)
       ul_j(3) = eFrame_Col(9,atom2)
#else
       ul_j(1) = eFrame(3,atom2)
       ul_j(2) = eFrame(6,atom2)
       ul_j(3) = eFrame(9,atom2)
#endif
       ct_j = ul_j(1)*drdxj + ul_j(2)*drdyj + ul_j(3)*drdzj

       if (j_is_SplitDipole) then
          d_j = ElectrostaticMap(me2)%split_dipole_distance
       endif
    endif

    epot = 0.0_dp
    dudx = 0.0_dp
    dudy = 0.0_dp
    dudz = 0.0_dp

    duduix = 0.0_dp
    duduiy = 0.0_dp
    duduiz = 0.0_dp

    dudujx = 0.0_dp
    dudujy = 0.0_dp
    dudujz = 0.0_dp

    if (i_is_Charge) then

       if (j_is_Charge) then
          
          vterm = pre11 * q_i * q_j * riji
          vpair = vpair + vterm
          epot = epot + sw*vterm

          dudr  = - sw * vterm * riji

          dudx = dudx + dudr * drdxj
          dudy = dudy + dudr * drdyj
          dudz = dudz + dudr * drdzj
       
       endif

       if (j_is_Dipole) then

          ri2 = riji * riji
          ri3 = ri2 * riji

          pref = pre12 * q_i * mu_j
          vterm = pref * ct_j * riji * riji
          vpair = vpair + vterm
          epot = epot + sw * vterm

          dudx = dudx + pref * sw * ri3 * ( ul_j(1) + 3.0d0 * ct_j * drdxj)
          dudy = dudy + pref * sw * ri3 * ( ul_j(2) + 3.0d0 * ct_j * drdyj)
          dudz = dudz + pref * sw * ri3 * ( ul_j(3) + 3.0d0 * ct_j * drdzj)

          dudujx = dudujx - pref * sw * ri2 * drdxj
          dudujy = dudujy - pref * sw * ri2 * drdyj
          dudujz = dudujz - pref * sw * ri2 * drdzj
          
       endif
    endif
  
    if (i_is_Dipole) then 
       
       if (j_is_Charge) then

          ri2 = riji * riji
          ri3 = ri2 * riji

          pref = pre12 * q_j * mu_i
          vterm = pref * ct_i * riji * riji
          vpair = vpair + vterm
          epot = epot + sw * vterm

          dudx = dudx + pref * sw * ri3 * ( ul_i(1) - 3.0d0 * ct_i * drdxj)
          dudy = dudy + pref * sw * ri3 * ( ul_i(2) - 3.0d0 * ct_i * drdyj)
          dudz = dudz + pref * sw * ri3 * ( ul_i(3) - 3.0d0 * ct_i * drdzj)

          duduix = duduix + pref * sw * ri2 * drdxj
          duduiy = duduiy + pref * sw * ri2 * drdyj
          duduiz = duduiz + pref * sw * ri2 * drdzj
       endif

       if (j_is_Dipole) then

          ct_ij = ul_i(1)*ul_j(1) + ul_i(2)*ul_j(2) + ul_i(3)*ul_j(3)
          ri2 = riji * riji
          ri3 = ri2 * riji
          ri4 = ri2 * ri2

          pref = pre22 * mu_i * mu_j
          vterm = pref * ri3 * (ct_ij - 3.0d0 * ct_i * ct_j)
          vpair = vpair + vterm
          epot = epot + sw * vterm
          
          a1 = 5.0d0 * ct_i * ct_j - ct_ij

          dudx = dudx + pref*sw*3.0d0*ri4*(a1*drdxj-ct_i*ul_j(1)-ct_j*ul_i(1))
          dudy = dudy + pref*sw*3.0d0*ri4*(a1*drdyj-ct_i*ul_j(2)-ct_j*ul_i(2))
          dudz = dudz + pref*sw*3.0d0*ri4*(a1*drdzj-ct_i*ul_j(3)-ct_j*ul_i(3))

          duduix = duduix + pref*sw*ri3*(ul_j(1) - 3.0d0*ct_j*drdxj)
          duduiy = duduiy + pref*sw*ri3*(ul_j(2) - 3.0d0*ct_j*drdyj)
          duduiz = duduiz + pref*sw*ri3*(ul_j(3) - 3.0d0*ct_j*drdzj)

          dudujx = dudujx + pref*sw*ri3*(ul_i(1) - 3.0d0*ct_i*drdxj)
          dudujy = dudujy + pref*sw*ri3*(ul_i(2) - 3.0d0*ct_i*drdyj)
          dudujz = dudujz + pref*sw*ri3*(ul_i(3) - 3.0d0*ct_i*drdzj)
       endif

    endif
    
    if (do_pot) then
#ifdef IS_MPI 
       pot_row(atom1) = pot_row(atom1) + 0.5d0*epot
       pot_col(atom2) = pot_col(atom2) + 0.5d0*epot
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
       t_Row(1,atom1) = t_Row(1,atom1) - ul_i(2)*duduiz + ul_i(3)*duduiy
       t_Row(2,atom1) = t_Row(2,atom1) - ul_i(3)*duduix + ul_i(1)*duduiz
       t_Row(3,atom1) = t_Row(3,atom1) - ul_i(1)*duduiy + ul_i(2)*duduix
    endif

    if (j_is_Dipole .or. j_is_Quadrupole) then
       t_Col(1,atom2) = t_Col(1,atom2) - ul_j(2)*dudujz + ul_j(3)*dudujy
       t_Col(2,atom2) = t_Col(2,atom2) - ul_j(3)*dudujx + ul_j(1)*dudujz
       t_Col(3,atom2) = t_Col(3,atom2) - ul_j(1)*dudujy + ul_j(2)*dudujx
    endif

#else
    f(1,atom1) = f(1,atom1) + dudx
    f(2,atom1) = f(2,atom1) + dudy
    f(3,atom1) = f(3,atom1) + dudz
    
    f(1,atom2) = f(1,atom2) - dudx
    f(2,atom2) = f(2,atom2) - dudy
    f(3,atom2) = f(3,atom2) - dudz
    
    if (i_is_Dipole .or. i_is_Quadrupole) then
       t(1,atom1) = t(1,atom1) - ul_i(2)*duduiz + ul_i(3)*duduiy
       t(2,atom1) = t(2,atom1) - ul_i(3)*duduix + ul_i(1)*duduiz
       t(3,atom1) = t(3,atom1) - ul_i(1)*duduiy + ul_i(2)*duduix
    endif
       
    if (j_is_Dipole .or. j_is_Quadrupole) then
       t(1,atom2) = t(1,atom2) - ul_j(2)*dudujz + ul_j(3)*dudujy
       t(2,atom2) = t(2,atom2) - ul_j(3)*dudujx + ul_j(1)*dudujz
       t(3,atom2) = t(3,atom2) - ul_j(1)*dudujy + ul_j(2)*dudujx
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
  
end module electrostatic_module

