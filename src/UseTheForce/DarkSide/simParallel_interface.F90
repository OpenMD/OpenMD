#ifdef IS_MPI
subroutine setFsimParallel(thisComponentPlan, nAtomTags, atomTags, &
     nGroupTags, groupTags, status)

  use mpiSimulation
  
  !! Passed Arguments
  !! mpiComponentPlan struct from C
  type (mpiComponentPlan), intent(inout) :: thisComponentPlan
  !! Number of tags passed
  integer, intent(in) :: nAtomTags, nGroupTags
  !! Result status, 0 = normal, -1 = error
  integer, intent(out) :: status
  integer :: localStatus
  !! Global reference tag for local particles
  integer, dimension(nAtomTags), intent(inout) :: atomTags
  integer, dimension(nGroupTags), intent(inout) :: groupTags

  call setupSimParallel(thisComponentPlan, nAtomTags, atomTags, &
       nGroupTags, groupTags, status)  
  
end subroutine setFsimParallel 

#else

!! Dummy routine so that we don't have an empty compilation unit:

subroutine setFsimParallel(status)

  integer, intent(out) :: status
  status = 0
  return

end subroutine setFsimParallel

#endif
