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

!! Wrapper interfaces for mpi.
!! This wrapper should work independently of mpi implimentation.
!! We only provide an explicit interface for routines we use. 
!! Access to other functions is external
module oopseMPI
  use definitions, ONLY : dp
#ifdef IS_MPI
  implicit none
  PUBLIC !WARNING everything in this module is public
  

  include "mpif.h"

  private  :: dp


! interfaces for things that we use 
! These routines are not overloaded and only include one argument type.
interface
  SUBROUTINE MPI_INIT(IERROR)
    INTEGER IERROR
  END SUBROUTINE MPI_INIT
  
  SUBROUTINE MPI_FINALIZE(IERROR)
    INTEGER IERROR
  END SUBROUTINE MPI_FINALIZE
  
  SUBROUTINE MPI_BARRIER(COMM, IERROR) 
    INTEGER COMM, IERROR
  END SUBROUTINE MPI_BARRIER
  
  SUBROUTINE MPI_COMM_RANK(COMM, RANK, IERROR)
    INTEGER COMM, RANK, IERROR
  END SUBROUTINE MPI_COMM_RANK
  
  SUBROUTINE MPI_COMM_SIZE(COMM, SIZE, IERROR)
    INTEGER COMM, SIZE, IERROR
  END SUBROUTINE MPI_COMM_SIZE
  !
  SUBROUTINE MPI_COMM_SPLIT(COMM, COLOR, KEY, NEWCOMM, IERROR)
    INTEGER COMM, COLOR, KEY, NEWCOMM, IERROR
  END SUBROUTINE MPI_COMM_SPLIT
  
  SUBROUTINE MPI_GET_PROCESSOR_NAME( NAME, RESULTLEN, IERROR)
    CHARACTER(len=*) :: NAME
    INTEGER RESULTLEN,IERROR
  END SUBROUTINE MPI_GET_PROCESSOR_NAME

!  FUNCTION MPI_WTICK()
!    DOUBLE PRECISION MPI_WTICK
!  END FUNCTION MPI_WTICK
end interface

!! These routines are overloaded and require multiple argument types
  interface mpi_allreduce
     module procedure mpi_allreduce_int
     module procedure mpi_allreduce_int_1d
     module procedure mpi_allreduce_int_2d
     module procedure mpi_allreduce_dp     
     module procedure mpi_allreduce_dp_1d     
     module procedure mpi_allreduce_dp_2d     
  end interface

!  interface mpi_reduce
!     module procedure mpi_reduce_int
!     module procedure mpi_reduce_int_1d
!     module procedure mpi_reduce_int_2d
!     module procedure mpi_reduce_dp     
!     module procedure mpi_reduce_dp_1d     
!     module procedure mpi_reduce_dp_2d     
!  end interface
 
  interface mpi_reduce_scatter
     module procedure mpi_reduce_scatter_int
     module procedure mpi_reduce_scatter_int_1d
     module procedure mpi_reduce_scatter_int_2d
     module procedure mpi_reduce_scatter_dp     
     module procedure mpi_reduce_scatter_dp_1d     
     module procedure mpi_reduce_scatter_dp_2d    
  end interface

  interface mpi_allgatherv
     module procedure mpi_allgatherv_int
     module procedure mpi_allgatherv_int_1d
     module procedure mpi_allgatherv_int_2d
     module procedure mpi_allgatherv_dp     
     module procedure mpi_allgatherv_dp_1d     
     module procedure mpi_allgatherv_dp_2d    
  end interface

  interface mpi_allgather
     module procedure mpi_allgather_int
     module procedure mpi_allgather_int_1d
     module procedure mpi_allgather_int_2d
     module procedure mpi_allgather_dp     
     module procedure mpi_allgather_dp_1d     
     module procedure mpi_allgather_dp_2d    
  end interface

  interface mpi_send
     module procedure mpi_send_dp
     module procedure mpi_send_dp_1d
     module procedure mpi_send_dp_2d
     module procedure mpi_send_int
     module procedure mpi_send_int_1d
     module procedure mpi_send_int_2d
     module procedure mpi_send_logical
     module procedure mpi_send_logical_1d
     module procedure mpi_send_t
     module procedure mpi_send_char
  end interface

  interface mpi_bcast
     module procedure mpi_bcast_dp
     module procedure mpi_bcast_dp_1d
     module procedure mpi_bcast_dp_2d
     module procedure mpi_bcast_int
     module procedure mpi_bcast_int_1d
     module procedure mpi_bcast_int_2d
     module procedure mpi_bcast_t
     module procedure mpi_bcast_char
     module procedure mpi_bcast_logical
     module procedure mpi_bcast_logical_1d
  end interface 

  interface mpi_recv
     module procedure mpi_recv_int
     module procedure mpi_recv_int_1d
     module procedure mpi_recv_int_2d
     module procedure mpi_recv_dp
     module procedure mpi_recv_dp_1d
     module procedure mpi_recv_dp_2d
     module procedure mpi_recv_t
     module procedure mpi_recv_char
     module procedure mpi_recv_logical
     module procedure mpi_recv_logical_1d
 end interface



contains



!! MPI BCAST FUNCTIONS
  subroutine mpi_bcast_t(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
        character (len=*),dimension(:) :: BUFFER
        integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
        external MPI_BCAST
        call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_t

  subroutine mpi_bcast_char(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
        character (len=*) :: BUFFER
        integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
        external MPI_BCAST
        call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_char

  subroutine mpi_bcast_int(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
        integer :: BUFFER
        integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
        external MPI_BCAST
        call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_int

  subroutine mpi_bcast_int_1d(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
        integer, dimension(:) :: BUFFER
        integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
        external MPI_BCAST
        call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_int_1d

  subroutine mpi_bcast_int_2d(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
        integer, dimension(:,:) :: BUFFER
        integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
        external MPI_BCAST
        call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_int_2d

  subroutine mpi_bcast_dp(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
        real(kind = dp) :: BUFFER
        integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
        external MPI_BCAST
        call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_dp

  subroutine mpi_bcast_dp_1d(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
        real(kind = dp),dimension(:) :: BUFFER
        integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
        external MPI_BCAST
        call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_dp_1d

  subroutine mpi_bcast_dp_2d(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
        real(kind = dp),dimension(:,:) :: BUFFER
        integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
        external MPI_BCAST
        call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_dp_2d

  subroutine mpi_bcast_logical(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
    logical :: BUFFER
    integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
    external MPI_BCAST
    call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_logical

  subroutine mpi_bcast_logical_1d(BUFFER, COUNT, DATATYPE, ROOT, COMM,    &
       IERROR)
    logical,dimension(:) :: BUFFER
    integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
    external MPI_BCAST
    call MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
  end subroutine mpi_bcast_logical_1d



!---------------------END MPIBCAST---------------------------------


!--------------------MPISEND-------------------------------------
  SUBROUTINE MPI_SEND_T(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    character(len=*), dimension(:) ::  BUF
    INTEGER  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_T
  SUBROUTINE MPI_SEND_CHAR(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    character(len=*) ::  BUF
    INTEGER  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_CHAR

  SUBROUTINE MPI_SEND_DP(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    real(kind=dp) ::  BUF
    INTEGER  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_DP

  SUBROUTINE MPI_SEND_DP_1D(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    real(kind=dp), dimension(:) ::  BUF
    INTEGER  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_DP_1D

  SUBROUTINE MPI_SEND_DP_2D(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    real(kind=dp), dimension(:,:) ::  BUF
    INTEGER  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_DP_2D

  SUBROUTINE MPI_SEND_INT(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    INTEGER ::  BUF
    INTEGER ::  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_INT

  SUBROUTINE MPI_SEND_INT_1D(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    INTEGER, dimension(:) ::  BUF
    INTEGER  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_INT_1D

  SUBROUTINE MPI_SEND_INT_2D(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    INTEGER, dimension(:,:) ::  BUF
    INTEGER  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_INT_2D
 
  SUBROUTINE MPI_SEND_LOGICAL(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    LOGICAL ::  BUF
    INTEGER  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_LOGICAL

  SUBROUTINE MPI_SEND_LOGICAL_1D(BUF, COUNT, DATATYPE, DEST, TAG, COMM,   &
       IERROR) 
    LOGICAL,dimension(:) ::  BUF
    INTEGER  COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    EXTERNAL MPI_SEND
    CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  END SUBROUTINE MPI_SEND_LOGICAL_1D
! ----------------END MPISEND------------------------------>

!------------------BEGIN MPIRECV-------------------------->

  subroutine mpi_recv_T(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
         STATUS, IERROR)
        character(len=*), dimension(:) :: BUF
        INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
         STATUS(MPI_STATUS_SIZE), IERROR 
        EXTERNAL MPI_RECV
        CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
         IERROR) 
  end subroutine mpi_recv_T

  subroutine mpi_recv_char(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
         STATUS, IERROR)
        character(len=*) :: BUF
        INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
         STATUS(MPI_STATUS_SIZE), IERROR 
        EXTERNAL MPI_RECV
        CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
         IERROR) 
  end subroutine mpi_recv_char

  subroutine mpi_recv_int(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
         STATUS, IERROR)
        INTEGER :: BUF
        INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
         STATUS(MPI_STATUS_SIZE), IERROR 
        EXTERNAL MPI_RECV
        CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
         IERROR) 
  end subroutine mpi_recv_int

  subroutine mpi_recv_int_1d(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
         STATUS, IERROR)
        INTEGER, dimension(:) :: BUF
        INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
         STATUS(MPI_STATUS_SIZE), IERROR 
        EXTERNAL MPI_RECV
        CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
         IERROR) 
  end subroutine mpi_recv_int_1d
  subroutine mpi_recv_int_2d(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
       STATUS, IERROR)
    INTEGER, dimension(:,:) :: BUF
    INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
            STATUS(MPI_STATUS_SIZE), IERROR 
    EXTERNAL MPI_RECV
    CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
            IERROR) 
  end subroutine mpi_recv_int_2d

  subroutine mpi_recv_dp(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
       STATUS, IERROR)
    real(kind=dp) :: BUF
    INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
         STATUS(MPI_STATUS_SIZE), IERROR 
    EXTERNAL MPI_RECV
    CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
         IERROR) 
  end subroutine mpi_recv_dp

  subroutine mpi_recv_dp_1d(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
       STATUS, IERROR)
    real(kind=dp), dimension(:) :: BUF
    INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
         STATUS(MPI_STATUS_SIZE), IERROR 
    EXTERNAL MPI_RECV
    CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
         IERROR) 
  end subroutine mpi_recv_dp_1d

  subroutine mpi_recv_dp_2d(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
       STATUS, IERROR)
    real(kind=dp), dimension(:,:) :: BUF
    INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
         STATUS(MPI_STATUS_SIZE), IERROR 
    EXTERNAL MPI_RECV
    CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
         IERROR) 
  end subroutine mpi_recv_dp_2d

  subroutine mpi_recv_logical_1d(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
       STATUS, IERROR)
    logical, dimension(:) :: BUF
    INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
         STATUS(MPI_STATUS_SIZE), IERROR 
    EXTERNAL MPI_RECV
    CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
         IERROR) 
  end subroutine mpi_recv_logical_1d

  subroutine mpi_recv_logical(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, &
       STATUS, IERROR)
    logical :: BUF
    INTEGER  COUNT, DATATYPE, SOURCE, TAG, COMM,                   &
         STATUS(MPI_STATUS_SIZE), IERROR 
    EXTERNAL MPI_RECV
    CALL MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, &
         IERROR) 
  end subroutine mpi_recv_logical
!-------------------------END MPIRECV------------------------------

!-------------------------MPI_ALLREDUCE----------------------------

  SUBROUTINE MPI_ALLREDUCE_int(SENDBUF, RECVBUF, COUNT, DATATYPE,  &
          OP, COMM, IERROR) 
    INTEGER :: SENDBUF, RECVBUF 
    INTEGER COUNT, DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_ALLREDUCE
    CALL MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP,      &
            COMM, IERROR) 
  END SUBROUTINE MPI_ALLREDUCE_INT

  SUBROUTINE MPI_ALLREDUCE_INT_1d(SENDBUF, RECVBUF, COUNT, DATATYPE,  &
          OP, COMM, IERROR) 
    INTEGER,dimension(:) :: SENDBUF, RECVBUF 
    INTEGER COUNT, DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_ALLREDUCE
    CALL MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP,      &
            COMM, IERROR) 
  END SUBROUTINE MPI_ALLREDUCE_INT_1D

  SUBROUTINE MPI_ALLREDUCE_INT_2D(SENDBUF, RECVBUF, COUNT, DATATYPE,  &
          OP, COMM, IERROR) 
    integer,dimension(:,:) :: SENDBUF, RECVBUF 
    INTEGER COUNT, DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_ALLREDUCE
    CALL MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP,      &
            COMM, IERROR) 
  END SUBROUTINE MPI_ALLREDUCE_INT_2D

  SUBROUTINE MPI_ALLREDUCE_DP(SENDBUF, RECVBUF, COUNT, DATATYPE,  &
          OP, COMM, IERROR) 
    REAL(kind=dp) :: SENDBUF, RECVBUF 
    INTEGER COUNT, DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_ALLREDUCE
    CALL MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP,      &
            COMM, IERROR) 
  END SUBROUTINE MPI_ALLREDUCE_DP

  SUBROUTINE MPI_ALLREDUCE_DP_1d(SENDBUF, RECVBUF, COUNT, DATATYPE,  &
          OP, COMM, IERROR) 
    REAL(kind=dp),dimension(:) :: SENDBUF, RECVBUF 
    INTEGER COUNT, DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_ALLREDUCE
    CALL MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP,      &
            COMM, IERROR) 
  END SUBROUTINE MPI_ALLREDUCE_DP_1D

  SUBROUTINE MPI_ALLREDUCE_DP_2D(SENDBUF, RECVBUF, COUNT, DATATYPE,  &
          OP, COMM, IERROR) 
    real(kind=dp),dimension(:,:) :: SENDBUF, RECVBUF 
    INTEGER COUNT, DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_ALLREDUCE
    CALL MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP,      &
            COMM, IERROR) 
  END SUBROUTINE MPI_ALLREDUCE_DP_2D
!-----------------END MPI_ALLREDUCE-------------------------->

!----------------BEGIN MPI_REDUCE_SCATTER
  SUBROUTINE MPI_REDUCE_SCATTER_DP(SENDBUF, RECVBUF, RECVCOUNTS,  &
       DATATYPE, OP, COMM, IERROR) 
    real(kind=dp) :: SENDBUF, RECVBUF 
    INTEGER RECVCOUNTS(*), DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_REDUCE_SCATTER
    CALL MPI_REDUCE_SCATTER(SENDBUF, RECVBUF, RECVCOUNTS,          &
         DATATYPE, OP, COMM, IERROR) 
  END SUBROUTINE MPI_REDUCE_SCATTER_DP

  SUBROUTINE MPI_REDUCE_SCATTER_DP_1D(SENDBUF, RECVBUF, RECVCOUNTS,  &
       DATATYPE, OP, COMM, IERROR) 
    real(kind=dp),dimension(:) :: SENDBUF, RECVBUF 
    INTEGER RECVCOUNTS(*), DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_REDUCE_SCATTER
    CALL MPI_REDUCE_SCATTER(SENDBUF, RECVBUF, RECVCOUNTS,          &
         DATATYPE, OP, COMM, IERROR) 
  END SUBROUTINE MPI_REDUCE_SCATTER_DP_1D

  SUBROUTINE MPI_REDUCE_SCATTER_DP_2D(SENDBUF, RECVBUF, RECVCOUNTS,  &
       DATATYPE, OP, COMM, IERROR) 
    real(kind=dp),dimension(:,:) :: SENDBUF, RECVBUF 
    INTEGER RECVCOUNTS(*), DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_REDUCE_SCATTER
    CALL MPI_REDUCE_SCATTER(SENDBUF, RECVBUF, RECVCOUNTS,          &
         DATATYPE, OP, COMM, IERROR) 
  END SUBROUTINE MPI_REDUCE_SCATTER_DP_2D

  SUBROUTINE MPI_REDUCE_SCATTER_INT(SENDBUF, RECVBUF, RECVCOUNTS,  &
       DATATYPE, OP, COMM, IERROR) 
    integer :: SENDBUF, RECVBUF 
    INTEGER RECVCOUNTS(*), DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_REDUCE_SCATTER
    CALL MPI_REDUCE_SCATTER(SENDBUF, RECVBUF, RECVCOUNTS,          &
         DATATYPE, OP, COMM, IERROR) 
  END SUBROUTINE MPI_REDUCE_SCATTER_INT

  SUBROUTINE MPI_REDUCE_SCATTER_INT_1D(SENDBUF, RECVBUF, RECVCOUNTS,  &
       DATATYPE, OP, COMM, IERROR) 
    integer,dimension(:) :: SENDBUF, RECVBUF 
    INTEGER RECVCOUNTS(*), DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_REDUCE_SCATTER
    CALL MPI_REDUCE_SCATTER(SENDBUF, RECVBUF, RECVCOUNTS,          &
         DATATYPE, OP, COMM, IERROR) 
  END SUBROUTINE MPI_REDUCE_SCATTER_INT_1D

  SUBROUTINE MPI_REDUCE_SCATTER_INT_2D(SENDBUF, RECVBUF, RECVCOUNTS,  &
       DATATYPE, OP, COMM, IERROR) 
    integer,dimension(:,:) :: SENDBUF, RECVBUF 
    INTEGER RECVCOUNTS(*), DATATYPE, OP, COMM, IERROR
    EXTERNAL MPI_REDUCE_SCATTER
    CALL MPI_REDUCE_SCATTER(SENDBUF, RECVBUF, RECVCOUNTS,          &
         DATATYPE, OP, COMM, IERROR) 
  END SUBROUTINE MPI_REDUCE_SCATTER_INT_2D

!end ---------------------MPI_REDUCE_SCATTER----------------->

!BEGIN------------------- MPI_ALLGATHERV--------------------->
  SUBROUTINE MPI_ALLGATHERV_INT(SENDBUF, SENDCOUNT, SENDTYPE,         &
       RECVBUF, RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
    INTEGER :: SENDBUF
    INTEGER,dimension(:) :: RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNTS(*), DISPLS(*),         &
         RECVTYPE, COMM, IERROR 
    EXTERNAL MPI_ALLGATHERV
    CALL MPI_ALLGATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,        &
         RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHERV_INT

  SUBROUTINE MPI_ALLGATHERV_INT_1D(SENDBUF, SENDCOUNT, SENDTYPE,         &
       RECVBUF, RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
    INTEGER,dimension(:) :: SENDBUF, RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNTS(*), DISPLS(*),         &
         RECVTYPE, COMM, IERROR 
    EXTERNAL MPI_ALLGATHERV
    CALL MPI_ALLGATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,        &
         RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_AlLGATHERV_INT_1D

  SUBROUTINE MPI_ALLGATHERV_INT_2D(SENDBUF, SENDCOUNT, SENDTYPE,         &
       RECVBUF, RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
    INTEGER,dimension(:,:) :: SENDBUF, RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNTS(*), DISPLS(*),         &
         RECVTYPE, COMM, IERROR 
    EXTERNAL MPI_ALLGATHERV
    CALL MPI_ALLGATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,        &
         RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHERV_INT_2D

 SUBROUTINE MPI_ALLGATHERV_DP(SENDBUF, SENDCOUNT, SENDTYPE,         &
       RECVBUF, RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
    real(kind=dp) :: SENDBUF
    real(kind=dp),dimension(:) :: RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNTS(*), DISPLS(*),         &
         RECVTYPE, ROOT, COMM, IERROR 
    EXTERNAL MPI_ALLGATHERV
    CALL MPI_ALLGATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,        &
         RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHERV_DP

  SUBROUTINE MPI_ALLGATHERV_DP_1D(SENDBUF, SENDCOUNT, SENDTYPE,         &
       RECVBUF, RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
    real(kind=dp),dimension(:) :: SENDBUF, RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNTS(*), DISPLS(*),         &
         RECVTYPE, COMM, IERROR 
    EXTERNAL MPI_ALLGATHERV
    CALL MPI_ALLGATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,        &
         RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_AlLGATHERV_dp_1D

  SUBROUTINE MPI_ALLGATHERV_dp_2D(SENDBUF, SENDCOUNT, SENDTYPE,         &
       RECVBUF, RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
    REAL(kind=dp),dimension(:,:) :: SENDBUF, RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNTS(*), DISPLS(*),         &
         RECVTYPE, COMM, IERROR 
    EXTERNAL MPI_ALLGATHERV
    CALL MPI_ALLGATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,        &
         RECVCOUNTS, DISPLS, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHERV_DP_2D
!--------------------------end MPI_ALLGATHERV----------------------->

  SUBROUTINE MPI_ALLGATHER_DP(SENDBUF, SENDCOUNT, SENDTYPE,       &
         RECVBUF, RECVCOUNT, RECVTYPE, COMM, IERROR) 
    real(kind=dp) :: SENDBUF
    real(kind=dp), dimension(:) :: RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM, IERROR
    EXTERNAL MPI_ALLGATHER
    CALL MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,      &
         RECVCOUNT, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHER_DP

  SUBROUTINE MPI_ALLGATHER_DP_1D(SENDBUF, SENDCOUNT, SENDTYPE,       &
         RECVBUF, RECVCOUNT, RECVTYPE, COMM, IERROR) 
    real(kind=dp),dimension(:) :: SENDBUF, RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM, IERROR
    EXTERNAL MPI_ALLGATHER
    CALL MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,      &
         RECVCOUNT, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHER_DP_1D

  SUBROUTINE MPI_ALLGATHER_DP_2D(SENDBUF, SENDCOUNT, SENDTYPE,       &
         RECVBUF, RECVCOUNT, RECVTYPE, COMM, IERROR) 
    real(kind=dp),dimension(:,:) :: SENDBUF, RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM, IERROR
    EXTERNAL MPI_ALLGATHER
    CALL MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,      &
         RECVCOUNT, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHER_DP_2D

  SUBROUTINE MPI_ALLGATHER_INT(SENDBUF, SENDCOUNT, SENDTYPE,       &
         RECVBUF, RECVCOUNT, RECVTYPE, COMM, IERROR) 
    integer :: SENDBUF 
    integer,dimension(:) :: RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM, IERROR
    EXTERNAL MPI_ALLGATHER
    CALL MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,      &
         RECVCOUNT, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHER_INT

  SUBROUTINE MPI_ALLGATHER_INT_1D(SENDBUF, SENDCOUNT, SENDTYPE,       &
         RECVBUF, RECVCOUNT, RECVTYPE, COMM, IERROR) 
    integer,dimension(:) :: SENDBUF, RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM, IERROR
    EXTERNAL MPI_ALLGATHER
    CALL MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,      &
         RECVCOUNT, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHER_INT_1D

  SUBROUTINE MPI_ALLGATHER_INT_2D(SENDBUF, SENDCOUNT, SENDTYPE,       &
         RECVBUF, RECVCOUNT, RECVTYPE, COMM, IERROR) 
    integer,dimension(:,:) :: SENDBUF, RECVBUF 
    INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM, IERROR
    EXTERNAL MPI_ALLGATHER
    CALL MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF,      &
         RECVCOUNT, RECVTYPE, COMM, IERROR) 
  END SUBROUTINE MPI_ALLGATHER_INT_2D
!-----------------------END MPI_ALLGATHER---------------------------

#endif
end module oopseMPI
