!! Machine and compiler dependent definitions
!! Charles F. Vardeman II 2/26/02
!! PUBLIC VARIABLES
!! DEFAULT_INPUT    default standard input
!! DEFAULT_OUTPUT   default standard output
!! DEFAULT_ERROR    default standard error
!! SP               single precision type
!! DP               double precision type
!! MAX_UNITS        system dependend maximum number of open units

module definitions
  IMPLICIT NONE
  PUBLIC

!! Machine dependent input and output (fortran 2000 will fix this standard)
  INTEGER, PARAMETER :: DEFAULT_INPUT  = 5
  INTEGER, PARAMETER :: DEFAULT_OUTPUT = 6
  INTEGER, PARAMETER :: DEFAULT_ERROR  = 0

!! Various precision parameters
  
  INTEGER, PARAMETER :: SP = selected_real_kind(4)
  INTEGER, PARAMETER :: DP = selected_real_kind(8)


!! Maximum number of fortran streams...
  INTEGER, PARAMETER :: MAX_UNITS = 100

!! number of dimensions in simulation
  INTEGER, PARAMETER :: ndim = 3
!! Default Size parameter of nlists
  INTEGER, PARAMETER :: nlistPrefactor = 80

end module definitions
