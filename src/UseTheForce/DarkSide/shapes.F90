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


module shapes

  use force_globals
  use definitions
  use atype_module
  use vector_class
  use simulation
  use status
  use lj
#ifdef IS_MPI
  use mpiSimulation
#endif
  implicit none

  PRIVATE
  
  INTEGER, PARAMETER:: CHEBYSHEV_TN = 1
  INTEGER, PARAMETER:: CHEBYSHEV_UN = 2
  INTEGER, PARAMETER:: LAGUERRE     = 3
  INTEGER, PARAMETER:: HERMITE      = 4
  INTEGER, PARAMETER:: SH_COS       = 0
  INTEGER, PARAMETER:: SH_SIN       = 1

  logical, save :: haveShapeMap = .false.

  public :: do_shape_pair
  public :: newShapeType
  public :: complete_Shape_FF


  type, private :: Shape
     integer :: atid
     integer :: nContactFuncs 
     integer :: nRangeFuncs 
     integer :: nStrengthFuncs 
     integer :: bigL
     integer :: bigM
     integer, pointer, dimension(:) :: ContactFuncLValue             => null()
     integer, pointer, dimension(:) :: ContactFuncMValue             => null()
     integer, pointer, dimension(:) :: ContactFunctionType           => null()
     real(kind=dp), pointer, dimension(:) :: ContactFuncCoefficient  => null()
     integer, pointer, dimension(:) :: RangeFuncLValue               => null()
     integer, pointer, dimension(:) :: RangeFuncMValue               => null()
     integer, pointer, dimension(:) :: RangeFunctionType             => null()
     real(kind=dp), pointer, dimension(:) :: RangeFuncCoefficient    => null()
     integer, pointer, dimension(:) :: StrengthFuncLValue            => null()
     integer, pointer, dimension(:) :: StrengthFuncMValue            => null()
     integer, pointer, dimension(:) :: StrengthFunctionType          => null()
     real(kind=dp), pointer, dimension(:) :: StrengthFuncCoefficient => null()
     logical :: isLJ 
     real ( kind = dp )  :: epsilon 
     real ( kind = dp )  :: sigma 
  end type Shape
  
  type, private :: ShapeList
     integer :: n_shapes = 0
     integer :: currentShape = 0
     type (Shape), pointer :: Shapes(:)      => null()
     integer, pointer      :: atidToShape(:) => null()
  end type ShapeList
  
  type(ShapeList), save :: ShapeMap

  integer :: lmax

contains  

  subroutine newShapeType(nContactFuncs, ContactFuncLValue, &
       ContactFuncMValue, ContactFunctionType, ContactFuncCoefficient, &
       nRangeFuncs, RangeFuncLValue, RangeFuncMValue, RangeFunctionType, &
       RangeFuncCoefficient, nStrengthFuncs, StrengthFuncLValue, &
       StrengthFuncMValue, StrengthFunctionType, StrengthFuncCoefficient, &
       myATID, status)

    integer :: nContactFuncs 
    integer :: nRangeFuncs 
    integer :: nStrengthFuncs 
    integer :: shape_ident
    integer :: status
    integer :: myATID
    integer :: bigL
    integer :: bigM
    integer :: j, me, nShapeTypes, nLJTypes, ntypes, current, alloc_stat
    integer, pointer :: MatchList(:) => null()

    integer, dimension(nContactFuncs) :: ContactFuncLValue           
    integer, dimension(nContactFuncs) :: ContactFuncMValue           
    integer, dimension(nContactFuncs) :: ContactFunctionType         
    real(kind=dp), dimension(nContactFuncs) :: ContactFuncCoefficient
    integer, dimension(nRangeFuncs) :: RangeFuncLValue             
    integer, dimension(nRangeFuncs) :: RangeFuncMValue             
    integer, dimension(nRangeFuncs) :: RangeFunctionType           
    real(kind=dp), dimension(nRangeFuncs) :: RangeFuncCoefficient  
    integer, dimension(nStrengthFuncs) :: StrengthFuncLValue          
    integer, dimension(nStrengthFuncs) :: StrengthFuncMValue          
    integer, dimension(nStrengthFuncs) :: StrengthFunctionType        
    real(kind=dp), dimension(nStrengthFuncs) :: StrengthFuncCoefficient

    status = 0
    ! check to see if this is the first time into this routine...
    if (.not.associated(ShapeMap%Shapes)) then

       call getMatchingElementList(atypes, "is_Shape", .true., &
            nShapeTypes, MatchList)
       
       call getMatchingElementList(atypes, "is_LennardJones", .true., &
            nLJTypes, MatchList)
       
       ShapeMap%n_shapes = nShapeTypes + nLJTypes
       
       allocate(ShapeMap%Shapes(nShapeTypes + nLJTypes))
       
       ntypes = getSize(atypes)
       
       allocate(ShapeMap%atidToShape(0:ntypes))
    end if
    
    ShapeMap%currentShape = ShapeMap%currentShape + 1
    current = ShapeMap%currentShape

    call allocateShape(nContactFuncs, nRangeFuncs, nStrengthFuncs, &
         ShapeMap%Shapes(current), stat=alloc_stat)
    if (alloc_stat .ne. 0) then
       status = -1
       return
    endif

    call getElementProperty(atypes, myATID, 'c_ident', me)

    ShapeMap%atidToShape(me)                         = current
    ShapeMap%Shapes(current)%atid                    = me
    ShapeMap%Shapes(current)%nContactFuncs           = nContactFuncs
    ShapeMap%Shapes(current)%nRangeFuncs             = nRangeFuncs
    ShapeMap%Shapes(current)%nStrengthFuncs          = nStrengthFuncs
    ShapeMap%Shapes(current)%ContactFuncLValue       = ContactFuncLValue
    ShapeMap%Shapes(current)%ContactFuncMValue       = ContactFuncMValue
    ShapeMap%Shapes(current)%ContactFunctionType     = ContactFunctionType
    ShapeMap%Shapes(current)%ContactFuncCoefficient  = ContactFuncCoefficient
    ShapeMap%Shapes(current)%RangeFuncLValue         = RangeFuncLValue
    ShapeMap%Shapes(current)%RangeFuncMValue         = RangeFuncMValue
    ShapeMap%Shapes(current)%RangeFunctionType       = RangeFunctionType
    ShapeMap%Shapes(current)%RangeFuncCoefficient    = RangeFuncCoefficient
    ShapeMap%Shapes(current)%StrengthFuncLValue      = StrengthFuncLValue
    ShapeMap%Shapes(current)%StrengthFuncMValue      = StrengthFuncMValue
    ShapeMap%Shapes(current)%StrengthFunctionType    = StrengthFunctionType
    ShapeMap%Shapes(current)%StrengthFuncCoefficient = StrengthFuncCoefficient

    bigL = -1
    bigM = -1
    
    do j = 1, ShapeMap%Shapes(current)%nContactFuncs
       if (ShapeMap%Shapes(current)%ContactFuncLValue(j) .gt. bigL) then
          bigL = ShapeMap%Shapes(current)%ContactFuncLValue(j)
       endif
       if (ShapeMap%Shapes(current)%ContactFuncMValue(j) .gt. bigM) then
          bigM = ShapeMap%Shapes(current)%ContactFuncMValue(j)
       endif
    enddo
    do j = 1, ShapeMap%Shapes(current)%nRangeFuncs
       if (ShapeMap%Shapes(current)%RangeFuncLValue(j) .gt. bigL) then
          bigL = ShapeMap%Shapes(current)%RangeFuncLValue(j)
       endif
       if (ShapeMap%Shapes(current)%RangeFuncMValue(j) .gt. bigM) then
          bigM = ShapeMap%Shapes(current)%RangeFuncMValue(j)
       endif
    enddo
    do j = 1, ShapeMap%Shapes(current)%nStrengthFuncs
       if (ShapeMap%Shapes(current)%StrengthFuncLValue(j) .gt. bigL) then
          bigL = ShapeMap%Shapes(current)%StrengthFuncLValue(j)
       endif
       if (ShapeMap%Shapes(current)%StrengthFuncMValue(j) .gt. bigM) then
          bigM = ShapeMap%Shapes(current)%StrengthFuncMValue(j)
       endif
    enddo

    ShapeMap%Shapes(current)%bigL                    = bigL
    ShapeMap%Shapes(current)%bigM                    = bigM

  end subroutine newShapeType

  subroutine allocateShape(nContactFuncs, nRangeFuncs, nStrengthFuncs, &
       myShape, stat)

    integer, intent(in) :: nContactFuncs, nRangeFuncs, nStrengthFuncs
    type(Shape), intent(inout) :: myShape
    integer, intent(out) :: stat
    integer :: alloc_stat
 
    stat = 0
    if (associated(myShape%contactFuncLValue)) then
       deallocate(myShape%contactFuncLValue)
    endif
    allocate(myShape%contactFuncLValue(nContactFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%contactFuncMValue)) then
       deallocate(myShape%contactFuncMValue)
    endif
    allocate(myShape%contactFuncMValue(nContactFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%contactFunctionType)) then
       deallocate(myShape%contactFunctionType)
    endif
    allocate(myShape%contactFunctionType(nContactFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%contactFuncCoefficient)) then
       deallocate(myShape%contactFuncCoefficient)
    endif
    allocate(myShape%contactFuncCoefficient(nContactFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif

    if (associated(myShape%rangeFuncLValue)) then
       deallocate(myShape%rangeFuncLValue)
    endif
    allocate(myShape%rangeFuncLValue(nRangeFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%rangeFuncMValue)) then
       deallocate(myShape%rangeFuncMValue)
    endif
    allocate(myShape%rangeFuncMValue(nRangeFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%rangeFunctionType)) then
       deallocate(myShape%rangeFunctionType)
    endif
    allocate(myShape%rangeFunctionType(nRangeFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%rangeFuncCoefficient)) then
       deallocate(myShape%rangeFuncCoefficient)
    endif
    allocate(myShape%rangeFuncCoefficient(nRangeFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif

    if (associated(myShape%strengthFuncLValue)) then
       deallocate(myShape%strengthFuncLValue)
    endif
    allocate(myShape%strengthFuncLValue(nStrengthFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%strengthFuncMValue)) then
       deallocate(myShape%strengthFuncMValue)
    endif
    allocate(myShape%strengthFuncMValue(nStrengthFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%strengthFunctionType)) then
       deallocate(myShape%strengthFunctionType)
    endif
    allocate(myShape%strengthFunctionType(nStrengthFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%strengthFuncCoefficient)) then
       deallocate(myShape%strengthFuncCoefficient)
    endif
    allocate(myShape%strengthFuncCoefficient(nStrengthFuncs), stat=alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif

    return

  end subroutine allocateShape
    
  subroutine complete_Shape_FF(status)
    integer :: status
    integer :: i, j, l, m, lm, function_type
    real(kind=dp) :: thisDP, sigma
    integer :: alloc_stat, iTheta, iPhi, nSteps, nAtypes, thisIP, current
    logical :: thisProperty

    status = 0
    if (ShapeMap%currentShape == 0) then
       call handleError("init_Shape_FF", "No members in ShapeMap")
       status = -1
       return
    end if
    
    nAtypes = getSize(atypes)

    if (nAtypes == 0) then
       status = -1
       return
    end if

    ! atypes comes from c side
    do i = 0, nAtypes
       
       call getElementProperty(atypes, i, "is_LennardJones", thisProperty)
       
       if (thisProperty) then
          
          ShapeMap%currentShape = ShapeMap%currentShape + 1
          current = ShapeMap%currentShape

          call getElementProperty(atypes, i, "c_ident",  thisIP)
          ShapeMap%atidToShape(thisIP) = current
          ShapeMap%Shapes(current)%atid = thisIP

          ShapeMap%Shapes(current)%isLJ = .true.

          ShapeMap%Shapes(current)%epsilon = getEpsilon(thisIP)
          ShapeMap%Shapes(current)%sigma = getSigma(thisIP)
          
       endif
       
    end do

    haveShapeMap = .true.
    
  end subroutine complete_Shape_FF
    
  subroutine do_shape_pair(atom1, atom2, d, rij, r2, sw, vpair, fpair, &
       pot, A, f, t, do_pot)
    
    INTEGER, PARAMETER:: LMAX         = 64
    INTEGER, PARAMETER:: MMAX         = 64

    integer, intent(in) :: atom1, atom2
    real (kind=dp), intent(inout) :: rij, r2
    real (kind=dp), dimension(3), intent(in) :: d
    real (kind=dp), dimension(3), intent(inout) :: fpair
    real (kind=dp) :: pot, vpair, sw, dswdr
    real (kind=dp), dimension(9,nLocal) :: A
    real (kind=dp), dimension(3,nLocal) :: f
    real (kind=dp), dimension(3,nLocal) :: t
    logical, intent(in) :: do_pot

    real (kind=dp) :: r3, r5, rt2, rt3, rt5, rt6, rt11, rt12, rt126
    integer :: atid1, atid2, st1, st2
    integer :: l, m, lm, id1, id2, localError, function_type
    real (kind=dp) :: sigma_i, s_i, eps_i, sigma_j, s_j, eps_j
    real (kind=dp) :: coeff
    real (kind=dp) :: pot_temp

    real (kind=dp) :: dsigmaidx, dsigmaidy, dsigmaidz
    real (kind=dp) :: dsigmaidux, dsigmaiduy, dsigmaiduz
    real (kind=dp) :: dsigmajdx, dsigmajdy, dsigmajdz
    real (kind=dp) :: dsigmajdux, dsigmajduy, dsigmajduz

    real (kind=dp) :: dsidx, dsidy, dsidz
    real (kind=dp) :: dsidux, dsiduy, dsiduz
    real (kind=dp) :: dsjdx, dsjdy, dsjdz
    real (kind=dp) :: dsjdux, dsjduy, dsjduz

    real (kind=dp) :: depsidx, depsidy, depsidz
    real (kind=dp) :: depsidux, depsiduy, depsiduz
    real (kind=dp) :: depsjdx, depsjdy, depsjdz
    real (kind=dp) :: depsjdux, depsjduy, depsjduz

    real (kind=dp) :: xi, yi, zi, xj, yj, zj, xi2, yi2, zi2, xj2, yj2, zj2

    real (kind=dp) :: sti2, stj2

    real (kind=dp) :: proji, proji3, projj, projj3
    real (kind=dp) :: cti, ctj, cpi, cpj, spi, spj
    real (kind=dp) :: Phunc, sigma, s, eps, rtdenom, rt

    real (kind=dp) :: dctidx, dctidy, dctidz
    real (kind=dp) :: dctidux, dctiduy, dctiduz
    real (kind=dp) :: dctjdx, dctjdy, dctjdz
    real (kind=dp) :: dctjdux, dctjduy, dctjduz

    real (kind=dp) :: dcpidx, dcpidy, dcpidz
    real (kind=dp) :: dcpidux, dcpiduy, dcpiduz
    real (kind=dp) :: dcpjdx, dcpjdy, dcpjdz
    real (kind=dp) :: dcpjdux, dcpjduy, dcpjduz

    real (kind=dp) :: dspidx, dspidy, dspidz
    real (kind=dp) :: dspidux, dspiduy, dspiduz
    real (kind=dp) :: dspjdx, dspjdy, dspjdz
    real (kind=dp) :: dspjdux, dspjduy, dspjduz

    real (kind=dp) :: dPhuncdX, dPhuncdY, dPhuncdZ
    real (kind=dp) :: dPhuncdUx, dPhuncdUy, dPhuncdUz

    real (kind=dp) :: dsigmadxi, dsigmadyi, dsigmadzi
    real (kind=dp) :: dsigmaduxi, dsigmaduyi, dsigmaduzi
    real (kind=dp) :: dsigmadxj, dsigmadyj, dsigmadzj
    real (kind=dp) :: dsigmaduxj, dsigmaduyj, dsigmaduzj

    real (kind=dp) :: dsdxi, dsdyi, dsdzi
    real (kind=dp) :: dsduxi, dsduyi, dsduzi
    real (kind=dp) :: dsdxj, dsdyj, dsdzj
    real (kind=dp) :: dsduxj, dsduyj, dsduzj
    
    real (kind=dp) :: depsdxi, depsdyi, depsdzi
    real (kind=dp) :: depsduxi, depsduyi, depsduzi
    real (kind=dp) :: depsdxj, depsdyj, depsdzj
    real (kind=dp) :: depsduxj, depsduyj, depsduzj

    real (kind=dp) :: drtdxi, drtdyi, drtdzi
    real (kind=dp) :: drtduxi, drtduyi, drtduzi
    real (kind=dp) :: drtdxj, drtdyj, drtdzj
    real (kind=dp) :: drtduxj, drtduyj, drtduzj

    real (kind=dp) :: drdxi, drdyi, drdzi
    real (kind=dp) :: drduxi, drduyi, drduzi
    real (kind=dp) :: drdxj, drdyj, drdzj
    real (kind=dp) :: drduxj, drduyj, drduzj

    real (kind=dp) :: dvdxi, dvdyi, dvdzi
    real (kind=dp) :: dvduxi, dvduyi, dvduzi
    real (kind=dp) :: dvdxj, dvdyj, dvdzj
    real (kind=dp) :: dvduxj, dvduyj, dvduzj  

    real (kind=dp) :: fxi, fyi, fzi, fxj, fyj, fzj
    real (kind=dp) :: txi, tyi, tzi, txj, tyj, tzj
    real (kind=dp) :: fxii, fyii, fzii, fxij, fyij, fzij
    real (kind=dp) :: fxji, fyji, fzji, fxjj, fyjj, fzjj
    real (kind=dp) :: fxradial, fyradial, fzradial

    real (kind=dp) :: plm_i(0:LMAX,0:MMAX), dlm_i(0:LMAX,0:MMAX)
    real (kind=dp) :: plm_j(0:LMAX,0:MMAX), dlm_j(0:LMAX,0:MMAX)
    real (kind=dp) :: tm_i(0:MMAX), dtm_i(0:MMAX), um_i(0:MMAX), dum_i(0:MMAX)
    real (kind=dp) :: tm_j(0:MMAX), dtm_j(0:MMAX), um_j(0:MMAX), dum_j(0:MMAX)

    if (.not.haveShapeMap) then
       call handleError("calc_shape", "NO SHAPEMAP!!!!")
       return       
    endif
    
    !! We assume that the rotation matrices have already been calculated
    !! and placed in the A array.

    r3 = r2*rij
    r5 = r3*r2
    
    drdxi = -d(1) / rij
    drdyi = -d(2) / rij
    drdzi = -d(3) / rij

    drdxj = d(1) / rij
    drdyj = d(2) / rij
    drdzj = d(3) / rij
    
    ! find the atom type id (atid) for each atom:
#ifdef IS_MPI
    atid1 = atid_Row(atom1)
    atid2 = atid_Col(atom2)
#else
    atid1 = atid(atom1)
    atid2 = atid(atom2)
#endif

    ! use the atid to find the shape type (st) for each atom:
    st1 = ShapeMap%atidToShape(atid1)
    st2 = ShapeMap%atidToShape(atid2)

    if (ShapeMap%Shapes(st1)%isLJ) then

       sigma_i = ShapeMap%Shapes(st1)%sigma
       s_i = ShapeMap%Shapes(st1)%sigma
       eps_i = ShapeMap%Shapes(st1)%epsilon
       dsigmaidx = 0.0d0
       dsigmaidy = 0.0d0
       dsigmaidz = 0.0d0
       dsigmaidux = 0.0d0
       dsigmaiduy = 0.0d0
       dsigmaiduz = 0.0d0
       dsidx = 0.0d0
       dsidy = 0.0d0
       dsidz = 0.0d0
       dsidux = 0.0d0
       dsiduy = 0.0d0
       dsiduz = 0.0d0
       depsidx = 0.0d0
       depsidy = 0.0d0
       depsidz = 0.0d0
       depsidux = 0.0d0
       depsiduy = 0.0d0
       depsiduz = 0.0d0
    else

#ifdef IS_MPI
       ! rotate the inter-particle separation into the two different
       ! body-fixed coordinate systems:
       
       xi = A_row(1,atom1)*d(1) + A_row(2,atom1)*d(2) + A_row(3,atom1)*d(3)
       yi = A_row(4,atom1)*d(1) + A_row(5,atom1)*d(2) + A_row(6,atom1)*d(3)
       zi = A_row(7,atom1)*d(1) + A_row(8,atom1)*d(2) + A_row(9,atom1)*d(3)
       
#else
       ! rotate the inter-particle separation into the two different
       ! body-fixed coordinate systems:
       
       xi = a(1,atom1)*d(1) + a(2,atom1)*d(2) + a(3,atom1)*d(3)
       yi = a(4,atom1)*d(1) + a(5,atom1)*d(2) + a(6,atom1)*d(3)
       zi = a(7,atom1)*d(1) + a(8,atom1)*d(2) + a(9,atom1)*d(3)
       
#endif

       xi2 = xi*xi
       yi2 = yi*yi
       zi2 = zi*zi             
       cti = zi / rij

       if (cti .gt. 1.0_dp) cti = 1.0_dp
       if (cti .lt. -1.0_dp) cti = -1.0_dp

       dctidx = - zi * xi / r3
       dctidy = - zi * yi / r3
       dctidz = 1.0d0 / rij - zi2 / r3
       dctidux = - (zi * xi2) / r3
       dctiduy = - (zi * yi2) / r3
       dctiduz = zi / rij - (zi2 * zi) / r3

       ! this is an attempt to try to truncate the singularity when
       ! sin(theta) is near 0.0:

       sti2 = 1.0_dp - cti*cti
       if (dabs(sti2) .lt. 1.0d-12) then
          proji = sqrt(rij * 1.0d-12)
          dcpidx = 1.0d0 / proji
          dcpidy = 0.0d0
          dcpidux = xi / proji
          dcpiduy = 0.0d0
          dspidx = 0.0d0
          dspidy = 1.0d0 / proji
          dspidux = 0.0d0
          dspiduy = yi / proji
       else
          proji = sqrt(xi2 + yi2)
          proji3 = proji*proji*proji
          dcpidx = 1.0d0 / proji - xi2 / proji3
          dcpidy = - xi * yi / proji3
          dcpidux = xi / proji - (xi2 * xi) / proji3
          dcpiduy = - (xi * yi2) / proji3
          dspidx = - xi * yi / proji3
          dspidy = 1.0d0 / proji - yi2 / proji3
          dspidux = - (yi * xi2) / proji3
          dspiduy = yi / proji - (yi2 * yi) / proji3
       endif
       
       cpi = xi / proji
       dcpidz = 0.0d0
       dcpiduz = 0.0d0
       
       spi = yi / proji
       dspidz = 0.0d0
       dspiduz = 0.0d0

       call Associated_Legendre(cti, ShapeMap%Shapes(st1)%bigM, &
            ShapeMap%Shapes(st1)%bigL, LMAX, &
            plm_i, dlm_i)

       call Orthogonal_Polynomial(cpi, ShapeMap%Shapes(st1)%bigM, MMAX, &
            CHEBYSHEV_TN, tm_i, dtm_i)
       call Orthogonal_Polynomial(cpi, ShapeMap%Shapes(st1)%bigM, MMAX, &
            CHEBYSHEV_UN, um_i, dum_i)
       
       sigma_i = 0.0d0
       s_i = 0.0d0
       eps_i = 0.0d0
       dsigmaidx = 0.0d0
       dsigmaidy = 0.0d0
       dsigmaidz = 0.0d0
       dsigmaidux = 0.0d0
       dsigmaiduy = 0.0d0
       dsigmaiduz = 0.0d0
       dsidx = 0.0d0
       dsidy = 0.0d0
       dsidz = 0.0d0
       dsidux = 0.0d0
       dsiduy = 0.0d0
       dsiduz = 0.0d0
       depsidx = 0.0d0
       depsidy = 0.0d0
       depsidz = 0.0d0
       depsidux = 0.0d0
       depsiduy = 0.0d0
       depsiduz = 0.0d0

       do lm = 1, ShapeMap%Shapes(st1)%nContactFuncs
          l = ShapeMap%Shapes(st1)%ContactFuncLValue(lm)
          m = ShapeMap%Shapes(st1)%ContactFuncMValue(lm)
          coeff = ShapeMap%Shapes(st1)%ContactFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st1)%ContactFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_i(m)
             dPhuncdX = coeff * dtm_i(m) * dcpidx
             dPhuncdY = coeff * dtm_i(m) * dcpidy
             dPhuncdZ = coeff * dtm_i(m) * dcpidz
             dPhuncdUz = coeff * dtm_i(m) * dcpidux
             dPhuncdUy = coeff * dtm_i(m) * dcpiduy
             dPhuncdUz = coeff * dtm_i(m) * dcpiduz
          else
             Phunc = coeff * spi * um_i(m-1)
             dPhuncdX = coeff * (spi * dum_i(m-1) * dcpidx + dspidx *um_i(m-1))
             dPhuncdY = coeff * (spi * dum_i(m-1) * dcpidy + dspidy *um_i(m-1))
             dPhuncdZ = coeff * (spi * dum_i(m-1) * dcpidz + dspidz *um_i(m-1))
             dPhuncdUx = coeff*(spi * dum_i(m-1)*dcpidux + dspidux *um_i(m-1))
             dPhuncdUy = coeff*(spi * dum_i(m-1)*dcpiduy + dspiduy *um_i(m-1))
             dPhuncdUz = coeff*(spi * dum_i(m-1)*dcpiduz + dspiduz *um_i(m-1))
          endif

          sigma_i = sigma_i + plm_i(m,l)*Phunc

          dsigmaidx = dsigmaidx + plm_i(m,l)*dPhuncdX + &
               Phunc * dlm_i(m,l) * dctidx
          dsigmaidy = dsigmaidy + plm_i(m,l)*dPhuncdY + &
               Phunc * dlm_i(m,l) * dctidy
          dsigmaidz = dsigmaidz + plm_i(m,l)*dPhuncdZ + &
               Phunc * dlm_i(m,l) * dctidz
          
          dsigmaidux = dsigmaidux + plm_i(m,l)* dPhuncdUx + &
               Phunc * dlm_i(m,l) * dctidux
          dsigmaiduy = dsigmaiduy + plm_i(m,l)* dPhuncdUy + &
               Phunc * dlm_i(m,l) * dctiduy
          dsigmaiduz = dsigmaiduz + plm_i(m,l)* dPhuncdUz + &
               Phunc * dlm_i(m,l) * dctiduz

       end do

       do lm = 1, ShapeMap%Shapes(st1)%nRangeFuncs
          l = ShapeMap%Shapes(st1)%RangeFuncLValue(lm)
          m = ShapeMap%Shapes(st1)%RangeFuncMValue(lm)
          coeff = ShapeMap%Shapes(st1)%RangeFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st1)%RangeFunctionType(lm)
          
          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_i(m)
             dPhuncdX = coeff * dtm_i(m) * dcpidx
             dPhuncdY = coeff * dtm_i(m) * dcpidy
             dPhuncdZ = coeff * dtm_i(m) * dcpidz
             dPhuncdUz = coeff * dtm_i(m) * dcpidux
             dPhuncdUy = coeff * dtm_i(m) * dcpiduy
             dPhuncdUz = coeff * dtm_i(m) * dcpiduz
          else
             Phunc = coeff * spi * um_i(m-1)
             dPhuncdX = coeff * (spi * dum_i(m-1) * dcpidx + dspidx *um_i(m-1))
             dPhuncdY = coeff * (spi * dum_i(m-1) * dcpidy + dspidy *um_i(m-1))
             dPhuncdZ = coeff * (spi * dum_i(m-1) * dcpidz + dspidz *um_i(m-1))
             dPhuncdUx = coeff*(spi * dum_i(m-1)*dcpidux + dspidux *um_i(m-1))
             dPhuncdUy = coeff*(spi * dum_i(m-1)*dcpiduy + dspiduy *um_i(m-1))
             dPhuncdUz = coeff*(spi * dum_i(m-1)*dcpiduz + dspiduz *um_i(m-1))
          endif

          s_i = s_i + plm_i(m,l)*Phunc
          
          dsidx = dsidx + plm_i(m,l)*dPhuncdX + &
               Phunc * dlm_i(m,l) * dctidx
          dsidy = dsidy + plm_i(m,l)*dPhuncdY + &
               Phunc * dlm_i(m,l) * dctidy
          dsidz = dsidz + plm_i(m,l)*dPhuncdZ + &
               Phunc * dlm_i(m,l) * dctidz
          
          dsidux = dsidux + plm_i(m,l)* dPhuncdUx + &
               Phunc * dlm_i(m,l) * dctidux
          dsiduy = dsiduy + plm_i(m,l)* dPhuncdUy + &
               Phunc * dlm_i(m,l) * dctiduy
          dsiduz = dsiduz + plm_i(m,l)* dPhuncdUz + &
               Phunc * dlm_i(m,l) * dctiduz      

       end do
              
       do lm = 1, ShapeMap%Shapes(st1)%nStrengthFuncs
          l = ShapeMap%Shapes(st1)%StrengthFuncLValue(lm)
          m = ShapeMap%Shapes(st1)%StrengthFuncMValue(lm)
          coeff = ShapeMap%Shapes(st1)%StrengthFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st1)%StrengthFunctionType(lm)
          
          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_i(m)
             dPhuncdX = coeff * dtm_i(m) * dcpidx
             dPhuncdY = coeff * dtm_i(m) * dcpidy
             dPhuncdZ = coeff * dtm_i(m) * dcpidz
             dPhuncdUz = coeff * dtm_i(m) * dcpidux
             dPhuncdUy = coeff * dtm_i(m) * dcpiduy
             dPhuncdUz = coeff * dtm_i(m) * dcpiduz
          else
             Phunc = coeff * spi * um_i(m-1)
             dPhuncdX = coeff * (spi * dum_i(m-1) * dcpidx + dspidx *um_i(m-1))
             dPhuncdY = coeff * (spi * dum_i(m-1) * dcpidy + dspidy *um_i(m-1))
             dPhuncdZ = coeff * (spi * dum_i(m-1) * dcpidz + dspidz *um_i(m-1))
             dPhuncdUx = coeff*(spi * dum_i(m-1)*dcpidux + dspidux *um_i(m-1))
             dPhuncdUy = coeff*(spi * dum_i(m-1)*dcpiduy + dspiduy *um_i(m-1))
             dPhuncdUz = coeff*(spi * dum_i(m-1)*dcpiduz + dspiduz *um_i(m-1))
          endif

          eps_i = eps_i + plm_i(m,l)*Phunc
          
          depsidx = depsidx + plm_i(m,l)*dPhuncdX + &
               Phunc * dlm_i(m,l) * dctidx
          depsidy = depsidy + plm_i(m,l)*dPhuncdY + &
               Phunc * dlm_i(m,l) * dctidy
          depsidz = depsidz + plm_i(m,l)*dPhuncdZ + &
               Phunc * dlm_i(m,l) * dctidz
          
          depsidux = depsidux + plm_i(m,l)* dPhuncdUx + &
               Phunc * dlm_i(m,l) * dctidux
          depsiduy = depsiduy + plm_i(m,l)* dPhuncdUy + &
               Phunc * dlm_i(m,l) * dctiduy
          depsiduz = depsiduz + plm_i(m,l)* dPhuncdUz + &
               Phunc * dlm_i(m,l) * dctiduz      

       end do

    endif
       
       ! now do j:

    if (ShapeMap%Shapes(st2)%isLJ) then
       sigma_j = ShapeMap%Shapes(st2)%sigma
       s_j = ShapeMap%Shapes(st2)%sigma
       eps_j = ShapeMap%Shapes(st2)%epsilon
       dsigmajdx = 0.0d0
       dsigmajdy = 0.0d0
       dsigmajdz = 0.0d0
       dsigmajdux = 0.0d0
       dsigmajduy = 0.0d0
       dsigmajduz = 0.0d0
       dsjdx = 0.0d0
       dsjdy = 0.0d0
       dsjdz = 0.0d0
       dsjdux = 0.0d0
       dsjduy = 0.0d0
       dsjduz = 0.0d0
       depsjdx = 0.0d0
       depsjdy = 0.0d0
       depsjdz = 0.0d0
       depsjdux = 0.0d0
       depsjduy = 0.0d0
       depsjduz = 0.0d0
    else
       
#ifdef IS_MPI
       ! rotate the inter-particle separation into the two different
       ! body-fixed coordinate systems:
       ! negative sign because this is the vector from j to i:
       
       xj = -(A_Col(1,atom2)*d(1) + A_Col(2,atom2)*d(2) + A_Col(3,atom2)*d(3))
       yj = -(A_Col(4,atom2)*d(1) + A_Col(5,atom2)*d(2) + A_Col(6,atom2)*d(3))
       zj = -(A_Col(7,atom2)*d(1) + A_Col(8,atom2)*d(2) + A_Col(9,atom2)*d(3))
#else
       ! rotate the inter-particle separation into the two different
       ! body-fixed coordinate systems:
       ! negative sign because this is the vector from j to i:
       
       xj = -(a(1,atom2)*d(1) + a(2,atom2)*d(2) + a(3,atom2)*d(3))
       yj = -(a(4,atom2)*d(1) + a(5,atom2)*d(2) + a(6,atom2)*d(3))
       zj = -(a(7,atom2)*d(1) + a(8,atom2)*d(2) + a(9,atom2)*d(3))
#endif
       
       xj2 = xj*xj
       yj2 = yj*yj
       zj2 = zj*zj
       ctj = zj / rij
       
       if (ctj .gt. 1.0_dp) ctj = 1.0_dp
       if (ctj .lt. -1.0_dp) ctj = -1.0_dp

       dctjdx = - zj * xj / r3
       dctjdy = - zj * yj / r3
       dctjdz = 1.0d0 / rij - zj2 / r3
       dctjdux = - (zi * xj2) / r3
       dctjduy = - (zj * yj2) / r3
       dctjduz = zj / rij - (zj2 * zj) / r3
       
       ! this is an attempt to try to truncate the singularity when
       ! sin(theta) is near 0.0:

       stj2 = 1.0_dp - ctj*ctj
       if (dabs(stj2) .lt. 1.0d-12) then
          projj = sqrt(rij * 1.0d-12)
          dcpjdx = 1.0d0 / projj 
          dcpjdy = 0.0d0
          dcpjdux = xj / projj
          dcpjduy = 0.0d0
          dspjdx = 0.0d0
          dspjdy = 1.0d0 / projj
          dspjdux = 0.0d0
          dspjduy = yj / projj
       else
          projj = sqrt(xj2 + yj2)
          projj3 = projj*projj*projj
          dcpjdx = 1.0d0 / projj - xj2 / projj3
          dcpjdy = - xj * yj / projj3
          dcpjdux = xj / projj - (xj2 * xj) / projj3
          dcpjduy = - (xj * yj2) / projj3
          dspjdx = - xj * yj / projj3
          dspjdy = 1.0d0 / projj - yj2 / projj3
          dspjdux = - (yj * xj2) / projj3
          dspjduy = yj / projj - (yj2 * yj) / projj3
       endif

       cpj = xj / projj
       dcpjdz = 0.0d0
       dcpjduz = 0.0d0
       
       spj = yj / projj
       dspjdz = 0.0d0
       dspjduz = 0.0d0


       write(*,*) 'dcpdu = ' ,dcpidux, dcpiduy, dcpiduz
       write(*,*) 'dcpdu = ' ,dcpjdux, dcpjduy, dcpjduz
       call Associated_Legendre(ctj, ShapeMap%Shapes(st2)%bigM, &
            ShapeMap%Shapes(st2)%bigL, LMAX, &
            plm_j, dlm_j)
       
       call Orthogonal_Polynomial(cpj, ShapeMap%Shapes(st2)%bigM, MMAX, &
            CHEBYSHEV_TN, tm_j, dtm_j)
       call Orthogonal_Polynomial(cpj, ShapeMap%Shapes(st2)%bigM, MMAX, &
            CHEBYSHEV_UN, um_j, dum_j)
       
       sigma_j = 0.0d0
       s_j = 0.0d0
       eps_j = 0.0d0
       dsigmajdx = 0.0d0
       dsigmajdy = 0.0d0
       dsigmajdz = 0.0d0
       dsigmajdux = 0.0d0
       dsigmajduy = 0.0d0
       dsigmajduz = 0.0d0
       dsjdx = 0.0d0
       dsjdy = 0.0d0
       dsjdz = 0.0d0
       dsjdux = 0.0d0
       dsjduy = 0.0d0
       dsjduz = 0.0d0
       depsjdx = 0.0d0
       depsjdy = 0.0d0
       depsjdz = 0.0d0
       depsjdux = 0.0d0
       depsjduy = 0.0d0
       depsjduz = 0.0d0

       do lm = 1, ShapeMap%Shapes(st2)%nContactFuncs
          l = ShapeMap%Shapes(st2)%ContactFuncLValue(lm)
          m = ShapeMap%Shapes(st2)%ContactFuncMValue(lm)
          coeff = ShapeMap%Shapes(st2)%ContactFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st2)%ContactFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_j(m)
             dPhuncdX = coeff * dtm_j(m) * dcpjdx
             dPhuncdY = coeff * dtm_j(m) * dcpjdy
             dPhuncdZ = coeff * dtm_j(m) * dcpjdz
             dPhuncdUz = coeff * dtm_j(m) * dcpjdux
             dPhuncdUy = coeff * dtm_j(m) * dcpjduy
             dPhuncdUz = coeff * dtm_j(m) * dcpjduz
          else
             Phunc = coeff * spj * um_j(m-1)
             dPhuncdX = coeff * (spj * dum_j(m-1) * dcpjdx + dspjdx *um_j(m-1))
             dPhuncdY = coeff * (spj * dum_j(m-1) * dcpjdy + dspjdy *um_j(m-1))
             dPhuncdZ = coeff * (spj * dum_j(m-1) * dcpjdz + dspjdz *um_j(m-1))
             dPhuncdUx = coeff*(spj * dum_j(m-1)*dcpjdux + dspjdux *um_j(m-1))
             dPhuncdUy = coeff*(spj * dum_j(m-1)*dcpjduy + dspjduy *um_j(m-1))
             dPhuncdUz = coeff*(spj * dum_j(m-1)*dcpjduz + dspjduz *um_j(m-1))
          endif
 
          sigma_j = sigma_j + plm_j(m,l)*Phunc
          
          dsigmajdx = dsigmajdx + plm_j(m,l)*dPhuncdX + &
               Phunc * dlm_j(m,l) * dctjdx
          dsigmajdy = dsigmajdy + plm_j(m,l)*dPhuncdY + &
               Phunc * dlm_j(m,l) * dctjdy
          dsigmajdz = dsigmajdz + plm_j(m,l)*dPhuncdZ + &
               Phunc * dlm_j(m,l) * dctjdz
          
          dsigmajdux = dsigmajdux + plm_j(m,l)* dPhuncdUx + &
               Phunc * dlm_j(m,l) * dctjdux
          dsigmajduy = dsigmajduy + plm_j(m,l)* dPhuncdUy + &
               Phunc * dlm_j(m,l) * dctjduy
          dsigmajduz = dsigmajduz + plm_j(m,l)* dPhuncdUz + &
               Phunc * dlm_j(m,l) * dctjduz

       end do

       do lm = 1, ShapeMap%Shapes(st2)%nRangeFuncs
          l = ShapeMap%Shapes(st2)%RangeFuncLValue(lm)
          m = ShapeMap%Shapes(st2)%RangeFuncMValue(lm)
          coeff = ShapeMap%Shapes(st2)%RangeFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st2)%RangeFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_j(m)
             dPhuncdX = coeff * dtm_j(m) * dcpjdx
             dPhuncdY = coeff * dtm_j(m) * dcpjdy
             dPhuncdZ = coeff * dtm_j(m) * dcpjdz
             dPhuncdUz = coeff * dtm_j(m) * dcpjdux
             dPhuncdUy = coeff * dtm_j(m) * dcpjduy
             dPhuncdUz = coeff * dtm_j(m) * dcpjduz
          else
             Phunc = coeff * spj * um_j(m-1)
             dPhuncdX = coeff * (spj * dum_j(m-1) * dcpjdx + dspjdx *um_j(m-1))
             dPhuncdY = coeff * (spj * dum_j(m-1) * dcpjdy + dspjdy *um_j(m-1))
             dPhuncdZ = coeff * (spj * dum_j(m-1) * dcpjdz + dspjdz *um_j(m-1))
             dPhuncdUx = coeff*(spj * dum_j(m-1)*dcpjdux + dspjdux *um_j(m-1))
             dPhuncdUy = coeff*(spj * dum_j(m-1)*dcpjduy + dspjduy *um_j(m-1))
             dPhuncdUz = coeff*(spj * dum_j(m-1)*dcpjduz + dspjduz *um_j(m-1))
          endif

          s_j = s_j + plm_j(m,l)*Phunc
          
          dsjdx = dsjdx + plm_j(m,l)*dPhuncdX + &
               Phunc * dlm_j(m,l) * dctjdx
          dsjdy = dsjdy + plm_j(m,l)*dPhuncdY + &
               Phunc * dlm_j(m,l) * dctjdy
          dsjdz = dsjdz + plm_j(m,l)*dPhuncdZ + &
               Phunc * dlm_j(m,l) * dctjdz
          
          dsjdux = dsjdux + plm_j(m,l)* dPhuncdUx + &
               Phunc * dlm_j(m,l) * dctjdux
          dsjduy = dsjduy + plm_j(m,l)* dPhuncdUy + &
               Phunc * dlm_j(m,l) * dctjduy
          dsjduz = dsjduz + plm_j(m,l)* dPhuncdUz + &
               Phunc * dlm_j(m,l) * dctjduz

       end do

       do lm = 1, ShapeMap%Shapes(st2)%nStrengthFuncs
          l = ShapeMap%Shapes(st2)%StrengthFuncLValue(lm)
          m = ShapeMap%Shapes(st2)%StrengthFuncMValue(lm)
          coeff = ShapeMap%Shapes(st2)%StrengthFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st2)%StrengthFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_j(m)
             dPhuncdX = coeff * dtm_j(m) * dcpjdx
             dPhuncdY = coeff * dtm_j(m) * dcpjdy
             dPhuncdZ = coeff * dtm_j(m) * dcpjdz
             dPhuncdUz = coeff * dtm_j(m) * dcpjdux
             dPhuncdUy = coeff * dtm_j(m) * dcpjduy
             dPhuncdUz = coeff * dtm_j(m) * dcpjduz
          else
             Phunc = coeff * spj * um_j(m-1)
             dPhuncdX = coeff * (spj * dum_j(m-1) * dcpjdx + dspjdx *um_j(m-1))
             dPhuncdY = coeff * (spj * dum_j(m-1) * dcpjdy + dspjdy *um_j(m-1))
             dPhuncdZ = coeff * (spj * dum_j(m-1) * dcpjdz + dspjdz *um_j(m-1))
             dPhuncdUx = coeff*(spj * dum_j(m-1)*dcpjdux + dspjdux *um_j(m-1))
             dPhuncdUy = coeff*(spj * dum_j(m-1)*dcpjduy + dspjduy *um_j(m-1))
             dPhuncdUz = coeff*(spj * dum_j(m-1)*dcpjduz + dspjduz *um_j(m-1))
          endif

          write(*,*) 'l,m = ', l, m, coeff, dPhuncdUx, dPhuncdUy, dPhuncdUz

          eps_j = eps_j + plm_j(m,l)*Phunc
          
          depsjdx = depsjdx + plm_j(m,l)*dPhuncdX + &
               Phunc * dlm_j(m,l) * dctjdx
          depsjdy = depsjdy + plm_j(m,l)*dPhuncdY + &
               Phunc * dlm_j(m,l) * dctjdy
          depsjdz = depsjdz + plm_j(m,l)*dPhuncdZ + &
               Phunc * dlm_j(m,l) * dctjdz
          
          depsjdux = depsjdux + plm_j(m,l)* dPhuncdUx + &
               Phunc * dlm_j(m,l) * dctjdux
          depsjduy = depsjduy + plm_j(m,l)* dPhuncdUy + &
               Phunc * dlm_j(m,l) * dctjduy
          depsjduz = depsjduz + plm_j(m,l)* dPhuncdUz + &
               Phunc * dlm_j(m,l) * dctjduz

       end do

    endif

    ! phew, now let's assemble the potential energy:

    sigma = 0.5*(sigma_i + sigma_j)

    dsigmadxi = 0.5*dsigmaidx
    dsigmadyi = 0.5*dsigmaidy
    dsigmadzi = 0.5*dsigmaidz
    dsigmaduxi = 0.5*dsigmaidux
    dsigmaduyi = 0.5*dsigmaiduy
    dsigmaduzi = 0.5*dsigmaiduz

    dsigmadxj = 0.5*dsigmajdx
    dsigmadyj = 0.5*dsigmajdy
    dsigmadzj = 0.5*dsigmajdz
    dsigmaduxj = 0.5*dsigmajdux
    dsigmaduyj = 0.5*dsigmajduy
    dsigmaduzj = 0.5*dsigmajduz

    s = 0.5*(s_i + s_j)

    dsdxi = 0.5*dsidx
    dsdyi = 0.5*dsidy
    dsdzi = 0.5*dsidz
    dsduxi = 0.5*dsidux
    dsduyi = 0.5*dsiduy
    dsduzi = 0.5*dsiduz

    dsdxj = 0.5*dsjdx
    dsdyj = 0.5*dsjdy
    dsdzj = 0.5*dsjdz
    dsduxj = 0.5*dsjdux
    dsduyj = 0.5*dsjduy
    dsduzj = 0.5*dsjduz

    eps = sqrt(eps_i * eps_j)

    depsdxi = eps_j * depsidx / (2.0d0 * eps)
    depsdyi = eps_j * depsidy / (2.0d0 * eps)
    depsdzi = eps_j * depsidz / (2.0d0 * eps)
    depsduxi = eps_j * depsidux / (2.0d0 * eps)
    depsduyi = eps_j * depsiduy / (2.0d0 * eps)
    depsduzi = eps_j * depsiduz / (2.0d0 * eps)

    depsdxj = eps_i * depsjdx / (2.0d0 * eps)
    depsdyj = eps_i * depsjdy / (2.0d0 * eps)
    depsdzj = eps_i * depsjdz / (2.0d0 * eps)
    depsduxj = eps_i * depsjdux / (2.0d0 * eps)
    depsduyj = eps_i * depsjduy / (2.0d0 * eps)
    depsduzj = eps_i * depsjduz / (2.0d0 * eps)
    
!!$    write(*,*) 'depsidu = ', depsidux, depsiduy, depsiduz
!!$    write(*,*) 'depsjdu = ', depsjdux, depsjduy, depsjduz
!!$
!!$    write(*,*) 'depsdui = ', depsduxi, depsduyi, depsduzi
!!$    write(*,*) 'depsduj = ', depsduxj, depsduyj, depsduzj
!!$
!!$    write(*,*) 's, sig, eps = ', s, sigma, eps

    rtdenom = rij-sigma+s
    rt = s / rtdenom

    drtdxi = (dsdxi + rt * (drdxi - dsigmadxi + dsdxi)) / rtdenom
    drtdyi = (dsdyi + rt * (drdyi - dsigmadyi + dsdyi)) / rtdenom
    drtdzi = (dsdzi + rt * (drdzi - dsigmadzi + dsdzi)) / rtdenom
    drtduxi = (dsduxi + rt * (drduxi - dsigmaduxi + dsduxi)) / rtdenom
    drtduyi = (dsduyi + rt * (drduyi - dsigmaduyi + dsduyi)) / rtdenom
    drtduzi = (dsduzi + rt * (drduzi - dsigmaduzi + dsduzi)) / rtdenom
    drtdxj = (dsdxj + rt * (drdxj - dsigmadxj + dsdxj)) / rtdenom
    drtdyj = (dsdyj + rt * (drdyj - dsigmadyj + dsdyj)) / rtdenom
    drtdzj = (dsdzj + rt * (drdzj - dsigmadzj + dsdzj)) / rtdenom
    drtduxj = (dsduxj + rt * (drduxj - dsigmaduxj + dsduxj)) / rtdenom
    drtduyj = (dsduyj + rt * (drduyj - dsigmaduyj + dsduyj)) / rtdenom
    drtduzj = (dsduzj + rt * (drduzj - dsigmaduzj + dsduzj)) / rtdenom
    
    rt2 = rt*rt
    rt3 = rt2*rt
    rt5 = rt2*rt3
    rt6 = rt3*rt3
    rt11 = rt5*rt6
    rt12 = rt6*rt6
    rt126 = rt12 - rt6

    pot_temp = 4.0d0 * eps * rt126

    vpair = vpair + pot_temp
    if (do_pot) then
#ifdef IS_MPI
       pot_row(atom1) = pot_row(atom1) + 0.5d0*pot_temp*sw
       pot_col(atom2) = pot_col(atom2) + 0.5d0*pot_temp*sw
#else
       pot = pot + pot_temp*sw
#endif
    endif

!!$    write(*,*) 'drtdu, depsdu = ', drtduxi, depsduxi
    
    dvdxi = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtdxi + 4.0d0*depsdxi*rt126
    dvdyi = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtdyi + 4.0d0*depsdyi*rt126
    dvdzi = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtdzi + 4.0d0*depsdzi*rt126
    dvduxi = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtduxi + 4.0d0*depsduxi*rt126
    dvduyi = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtduyi + 4.0d0*depsduyi*rt126
    dvduzi = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtduzi + 4.0d0*depsduzi*rt126

    dvdxj = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtdxj + 4.0d0*depsdxj*rt126
    dvdyj = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtdyj + 4.0d0*depsdyj*rt126
    dvdzj = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtdzj + 4.0d0*depsdzj*rt126
    dvduxj = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtduxj + 4.0d0*depsduxj*rt126
    dvduyj = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtduyj + 4.0d0*depsduyj*rt126
    dvduzj = 24.0d0*eps*(2.0d0*rt11 - rt5)*drtduzj + 4.0d0*depsduzj*rt126

    ! do the torques first since they are easy:
    ! remember that these are still in the body fixed axes


!!$    write(*,*) 'sw = ', sw
!!$    write(*,*) 'dvdu1 = ', dvduxi, dvduyi, dvduzi
!!$    write(*,*) 'dvdu2 = ', dvduxj, dvduyj, dvduzj
!!$
    txi =  (dvduzi - dvduyi) * sw
    tyi =  (dvduxi - dvduzi) * sw
    tzi =  (dvduyi - dvduxi) * sw

    txj = (dvduzj - dvduyj) * sw
    tyj = (dvduxj - dvduzj) * sw
    tzj = (dvduyj - dvduxj) * sw

!!$    txi = -dvduxi * sw
!!$    tyi = -dvduyi * sw
!!$    tzi = -dvduzi * sw
!!$
!!$    txj = dvduxj * sw
!!$    tyj = dvduyj * sw
!!$    tzj = dvduzj * sw
 
    write(*,*) 't1 = ', txi, tyi, tzi
    write(*,*) 't2 = ', txj, tyj, tzj

    ! go back to lab frame using transpose of rotation matrix:
    
#ifdef IS_MPI
    t_Row(1,atom1) = t_Row(1,atom1) + a_Row(1,atom1)*txi + &
         a_Row(4,atom1)*tyi + a_Row(7,atom1)*tzi
    t_Row(2,atom1) = t_Row(2,atom1) + a_Row(2,atom1)*txi + &
         a_Row(5,atom1)*tyi + a_Row(8,atom1)*tzi
    t_Row(3,atom1) = t_Row(3,atom1) + a_Row(3,atom1)*txi + &
         a_Row(6,atom1)*tyi + a_Row(9,atom1)*tzi
    
    t_Col(1,atom2) = t_Col(1,atom2) + a_Col(1,atom2)*txj + &
         a_Col(4,atom2)*tyj + a_Col(7,atom2)*tzj
    t_Col(2,atom2) = t_Col(2,atom2) + a_Col(2,atom2)*txj + &
            a_Col(5,atom2)*tyj + a_Col(8,atom2)*tzj
    t_Col(3,atom2) = t_Col(3,atom2) + a_Col(3,atom2)*txj + &
         a_Col(6,atom2)*tyj + a_Col(9,atom2)*tzj
#else
    t(1,atom1) = t(1,atom1) + a(1,atom1)*txi + a(4,atom1)*tyi + a(7,atom1)*tzi
    t(2,atom1) = t(2,atom1) + a(2,atom1)*txi + a(5,atom1)*tyi + a(8,atom1)*tzi
    t(3,atom1) = t(3,atom1) + a(3,atom1)*txi + a(6,atom1)*tyi + a(9,atom1)*tzi
    
    t(1,atom2) = t(1,atom2) + a(1,atom2)*txj + a(4,atom2)*tyj + a(7,atom2)*tzj
    t(2,atom2) = t(2,atom2) + a(2,atom2)*txj + a(5,atom2)*tyj + a(8,atom2)*tzj
    t(3,atom2) = t(3,atom2) + a(3,atom2)*txj + a(6,atom2)*tyj + a(9,atom2)*tzj
#endif
    ! Now, on to the forces:
    
    ! first rotate the i terms back into the lab frame:
    
    fxi = dvdxi * sw
    fyi = dvdyi * sw
    fzi = dvdzi * sw

    fxj = dvdxj * sw
    fyj = dvdyj * sw
    fzj = dvdzj * sw

#ifdef IS_MPI
    fxii = a_Row(1,atom1)*fxi + a_Row(4,atom1)*fyi + a_Row(7,atom1)*fzi
    fyii = a_Row(2,atom1)*fxi + a_Row(5,atom1)*fyi + a_Row(8,atom1)*fzi
    fzii = a_Row(3,atom1)*fxi + a_Row(6,atom1)*fyi + a_Row(9,atom1)*fzi

    fxjj = a_Col(1,atom2)*fxj + a_Col(4,atom2)*fyj + a_Col(7,atom2)*fzj
    fyjj = a_Col(2,atom2)*fxj + a_Col(5,atom2)*fyj + a_Col(8,atom2)*fzj
    fzjj = a_Col(3,atom2)*fxj + a_Col(6,atom2)*fyj + a_Col(9,atom2)*fzj
#else
    fxii = a(1,atom1)*fxi + a(4,atom1)*fyi + a(7,atom1)*fzi
    fyii = a(2,atom1)*fxi + a(5,atom1)*fyi + a(8,atom1)*fzi
    fzii = a(3,atom1)*fxi + a(6,atom1)*fyi + a(9,atom1)*fzi
    
    fxjj = a(1,atom2)*fxj + a(4,atom2)*fyj + a(7,atom2)*fzj
    fyjj = a(2,atom2)*fxj + a(5,atom2)*fyj + a(8,atom2)*fzj
    fzjj = a(3,atom2)*fxj + a(6,atom2)*fyj + a(9,atom2)*fzj
#endif

    fxij = -fxii
    fyij = -fyii
    fzij = -fzii
    
    fxji = -fxjj
    fyji = -fyjj
    fzji = -fzjj

    fxradial = 0.5_dp * (fxii + fxji)
    fyradial = 0.5_dp * (fyii + fyji)
    fzradial = 0.5_dp * (fzii + fzji)

#ifdef IS_MPI
    f_Row(1,atom1) = f_Row(1,atom1) + fxradial
    f_Row(2,atom1) = f_Row(2,atom1) + fyradial
    f_Row(3,atom1) = f_Row(3,atom1) + fzradial
    
    f_Col(1,atom2) = f_Col(1,atom2) - fxradial
    f_Col(2,atom2) = f_Col(2,atom2) - fyradial
    f_Col(3,atom2) = f_Col(3,atom2) - fzradial
#else
    f(1,atom1) = f(1,atom1) + fxradial
    f(2,atom1) = f(2,atom1) + fyradial
    f(3,atom1) = f(3,atom1) + fzradial
    
    f(1,atom2) = f(1,atom2) - fxradial
    f(2,atom2) = f(2,atom2) - fyradial
    f(3,atom2) = f(3,atom2) - fzradial
#endif

#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif
    
    if (molMembershipList(id1) .ne. molMembershipList(id2)) then
       
       fpair(1) = fpair(1) + fxradial
       fpair(2) = fpair(2) + fyradial
       fpair(3) = fpair(3) + fzradial
       
    endif

  end subroutine do_shape_pair
    
  SUBROUTINE Associated_Legendre(x, l, m, lmax, plm, dlm)        

    ! Purpose: Compute the associated Legendre functions 
    !          Plm(x) and their derivatives Plm'(x)
    ! Input :  x  --- Argument of Plm(x)
    !          l  --- Order of Plm(x),  l = 0,1,2,...,n
    !          m  --- Degree of Plm(x), m = 0,1,2,...,N
    !          lmax --- Physical dimension of PLM and DLM
    ! Output:  PLM(l,m) --- Plm(x)
    !          DLM(l,m) --- Plm'(x)
    !
    ! adapted from the routines in 
    ! COMPUTATION OF SPECIAL FUNCTIONS by Shanjie Zhang and Jianming Jin
    ! ISBN 0-471-11963-6
    !
    ! The original Fortran77 codes can be found here:
    ! http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
    
    real (kind=dp), intent(in) :: x
    integer, intent(in) :: l, m, lmax 
    real (kind=dp), dimension(0:lmax,0:m), intent(out) :: PLM, DLM
    integer :: i, j, ls
    real (kind=dp) :: xq, xs

    ! zero out both arrays:
    DO I = 0, m
       DO J = 0, l
          PLM(J,I) = 0.0_dp
          DLM(J,I) = 0.0_dp
       end DO
    end DO

    ! start with 0,0:
    PLM(0,0) = 1.0D0
  
    ! x = +/- 1 functions are easy:
    IF (abs(X).EQ.1.0D0) THEN
       DO I = 1, m
          PLM(0, I) = X**I
          DLM(0, I) = 0.5D0*I*(I+1.0D0)*X**(I+1)
       end DO
       DO J = 1, m
          DO I = 1, l
             IF (I.EQ.1) THEN
                DLM(I, J) = 1.0D+300
             ELSE IF (I.EQ.2) THEN
                DLM(I, J) = -0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
             ENDIF
          end DO
       end DO
       RETURN
    ENDIF

    LS = 1
    IF (abs(X).GT.1.0D0) LS = -1
    XQ = sqrt(LS*(1.0D0-X*X))
    XS = LS*(1.0D0-X*X)

    DO I = 1, l
       PLM(I, I) = -LS*(2.0D0*I-1.0D0)*XQ*PLM(I-1, I-1)
    enddo

    DO I = 0, l
       PLM(I, I+1)=(2.0D0*I+1.0D0)*X*PLM(I, I)
    enddo

    DO I = 0, l
       DO J = I+2, m
          PLM(I, J)=((2.0D0*J-1.0D0)*X*PLM(I,J-1) - &
               (I+J-1.0D0)*PLM(I,J-2))/(J-I)
       end DO
    end DO

    DLM(0, 0)=0.0D0
    DO J = 1, m
       DLM(0, J)=LS*J*(PLM(0,J-1)-X*PLM(0,J))/XS
    end DO

    DO I = 1, l
       DO J = I, m
          DLM(I,J) = LS*I*X*PLM(I, J)/XS + (J+I)*(J-I+1.0D0)/XQ*PLM(I-1, J)
       end DO
    end DO

    RETURN
  END SUBROUTINE Associated_Legendre


  subroutine Orthogonal_Polynomial(x, m, mmax, function_type, pl, dpl)
  
    ! Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
    !          or Ln(x) or Hn(x), and their derivatives
    ! Input :  function_type --- Function code
    !                 =1 for Chebyshev polynomial Tn(x)
    !                 =2 for Chebyshev polynomial Un(x)
    !                 =3 for Laguerre polynomial Ln(x)
    !                 =4 for Hermite polynomial Hn(x)
    !          n ---  Order of orthogonal polynomials
    !          x ---  Argument of orthogonal polynomials
    ! Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
    !          DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
    !
    ! adapted from the routines in 
    ! COMPUTATION OF SPECIAL FUNCTIONS by Shanjie Zhang and Jianming Jin
    ! ISBN 0-471-11963-6
    !
    ! The original Fortran77 codes can be found here:
    ! http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
  
    real(kind=8), intent(in) :: x
    integer, intent(in):: m, mmax
    integer, intent(in):: function_type
    real(kind=8), dimension(0:mmax), intent(inout) :: pl, dpl
  
    real(kind=8) :: a, b, c, y0, y1, dy0, dy1, yn, dyn
    integer :: k

    A = 2.0D0
    B = 0.0D0
    C = 1.0D0
    Y0 = 1.0D0
    Y1 = 2.0D0*X
    DY0 = 0.0D0
    DY1 = 2.0D0
    PL(0) = 1.0D0
    PL(1) = 2.0D0*X
    DPL(0) = 0.0D0
    DPL(1) = 2.0D0
    IF (function_type.EQ.CHEBYSHEV_TN) THEN
       Y1 = X
       DY1 = 1.0D0
       PL(1) = X
       DPL(1) = 1.0D0
    ELSE IF (function_type.EQ.LAGUERRE) THEN
       Y1 = 1.0D0-X
       DY1 = -1.0D0
       PL(1) = 1.0D0-X
       DPL(1) = -1.0D0
    ENDIF
    DO K = 2, m
       IF (function_type.EQ.LAGUERRE) THEN
          A = -1.0D0/K
          B = 2.0D0+A
          C = 1.0D0+A
       ELSE IF (function_type.EQ.HERMITE) THEN
          C = 2.0D0*(K-1.0D0)
       ENDIF
       YN = (A*X+B)*Y1-C*Y0
       DYN = A*Y1+(A*X+B)*DY1-C*DY0
       PL(K) = YN
       DPL(K) = DYN
       Y0 = Y1
       Y1 = YN
       DY0 = DY1
       DY1 = DYN
    end DO


    RETURN
    
  end subroutine Orthogonal_Polynomial
  
end module shapes
