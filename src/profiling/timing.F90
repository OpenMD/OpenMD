subroutine getTimes(forceTime,commTime)
      use doForces
      use definitions, ONLY: dp
#ifdef IS_MPI
      use mpiSimulation
#endif
      implicit none

      real(kind=dp) :: forceTime
      real(kind=dp) :: commTime
      
      forceTime = 0.0
      commTime  = 0.0
#ifdef PROFILE
      forceTime = getForceTime()

#ifdef IS_MPI
      commTime  = getCommTime()
#endif
#endif

end subroutine getTimes
