subroutine newEAMtype(lattice_constant,eam_nrho,eam_drho,eam_nr,&
     eam_dr,rcut,eam_Z_r,eam_rho_r,eam_F_rho,&
     eam_ident,status)
  use definitions, ONLY : dp
  use eam, ONLY : module_newEAMtype => newEAMtype
  real (kind = dp )                      :: lattice_constant
  integer                                :: eam_nrho
  real (kind = dp )                      :: eam_drho
  integer                                :: eam_nr
  real (kind = dp )                     :: eam_dr
  real (kind = dp )                      :: rcut
  real (kind = dp ), dimension(eam_nr)    :: eam_Z_r
  real (kind = dp ), dimension(eam_nr)   :: eam_rho_r
  real (kind = dp ), dimension(eam_nrho) :: eam_F_rho
  integer                                :: eam_ident
  integer                                :: status


  call module_newEAMtype(lattice_constant,eam_nrho,eam_drho,eam_nr,&
       eam_dr,rcut,eam_Z_r,eam_rho_r,eam_F_rho,&
       eam_ident,status)

end subroutine newEAMtype

subroutine destroyEAMTypes()
  use eam, ONLY: module_destroyEAMTypes =>destroyEAMTypes

  call module_destroyEAMTypes()  

end subroutine destroyEAMTypes
