module status
  implicit none
  PRIVATE

#define __FORTRAN90
#include "simError.h"

  character(len=1), parameter :: nullchar = char(0)
  character(len=1), parameter :: newline = char(10)
  character(len=1), parameter :: tab = char(9)
  INTEGER, PARAMETER:: statusMsgSize = MAX_SIM_ERROR_MSG_LENGTH 

!!$interface
!!$   
!!$   subroutine c_simError(painCave)
!!$     type(errorStruct), pointer :: painCave
!!$   end subroutine c_simError
!!$
!!$end interface

public :: handleInfo
public :: handleError
public :: handleWarning
public :: statusMsgSize
public :: nullchar
public :: newline
public :: tab

contains

  subroutine handleInfo(myRoutine, myMessage)
    character(len=*), intent(in) :: myRoutine
    character(len=*), intent(in) :: myMessage
   
    painCave%errMsg = "Location: " // trim(myRoutine) // newline // &
         tab // trim(myMessage) // newline // nullchar

    painCave%severity = OOPSE_INFO
    painCave%isFatal = .false.

    call c_simError(painCave)
  end subroutine handleInfo


  subroutine handleError(myRoutine, myMessage)
    character(len=*), intent(in) :: myRoutine
    character(len=*), intent(in) :: myMessage

    painCave%errMsg = "Location: " // trim(myRoutine) // newline // &
         tab // trim(myMessage) // newline // nullchar

    painCave%severity = OOPSE_ERROR
    painCave%isFatal = .true.

    call c_simError(painCave)
    
  end subroutine handleError

  subroutine handleWarning(myRoutine, myMessage)
    character(len=*), intent(in) :: myRoutine
    character(len=*), intent(in) :: myMessage

    painCave%errMsg = "Location: " // trim(myRoutine) // newline // &
         tab // trim(myMessage) // newline // nullchar

    painCave%severity = OOPSE_WARNING
    painCave%isFatal = .false.

    call c_simError(painCave)

  end subroutine handleWarning
  
end module status
