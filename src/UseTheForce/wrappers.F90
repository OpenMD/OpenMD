
subroutine wrapForceField(infoWrapper)

  use atype_module, ONLY: new_atype
  use doForces, ONLY: init_FF, do_force_loop
  use sticky_pair, ONLY: set_sticky_params
  use gb_pair, ONLY: set_gb_pair_params
  use eam, ONLY: newEAMtype
  external infoWrapper

  call infoWrapper(new_atype,init_FF,do_force_loop, &
       set_sticky_params,set_gb_pair_params,newEAMtype)
  
end subroutine wrapForceField


subroutine wrapSimMod(infoWrapper)
  use simulation, only: simulationSetup, setBox
  use notifyCutoffs, only: cutoffNotify
  
  external infoWrapper

  call infoWrapper( simulationSetup, setBox, cutoffNotify )

end subroutine wrapSimMod


#ifdef IS_MPI

subroutine wrapSimParallelMod(infoWrapper)
  use mpiSimulation, only: setupSimParallel
  
  external infoWrapper
  
  call infoWrapper( setupSimParallel )

end subroutine wrapSimParallelMod

#endif
