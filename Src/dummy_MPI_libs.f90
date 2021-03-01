!---------------------------------------------------------------------
!     *******  The following subroutines and functions are     *******
!     *******  to be used in the sequential environment.       *******
!---------------------------------------------------------------------


subroutine MPI_INIT(ierr)
  implicit none
  integer(kind=4) :: ierr
  return
end subroutine MPI_INIT
 
subroutine MPI_COMM_RANK(icomm, iproc, ierr)
  implicit none
  integer(kind=4) :: icomm, iproc, ierr
!-----------------------------------------
  iproc = 0
  return
end
 
subroutine MPI_COMM_SIZE(icomm, nproc, ierr)
  implicit none
  integer(kind=4) :: icomm, nproc, ierr
!-----------------------------------------
  nproc = 1
  return
end
 
subroutine MPI_Barrier(icomm,ierr)
  implicit none
  integer(kind=4) :: icomm, ierr
!-----------------------------------------
  return
end
 
 
subroutine MPI_Finalize(ierr)
  implicit none
  integer(kind=4) :: ierr
!-----------------------------------------
  return
end

subroutine MPI_Bcast()
  return
end subroutine MPI_Bcast
 
subroutine MPI_Reduce()
  return
end subroutine MPI_Reduce
 
subroutine MPI_Allreduce()
  return
end subroutine MPI_Allreduce
 
subroutine MPI_Abort()
  call exit
  return
end subroutine MPI_Abort

subroutine MPI_INFO_CREATE()
  return
end subroutine MPI_INFO_CREATE

subroutine MPI_INFO_SET()
  return
end subroutine MPI_INFO_SET

subroutine MPI_FILE_OPEN()
  return
end subroutine MPI_FILE_OPEN

subroutine MPI_FILE_SET_VIEW()
  return
end subroutine MPI_FILE_SET_VIEW

subroutine MPI_FILE_READ_AT()
  return
end subroutine MPI_FILE_READ_AT

subroutine MPI_FILE_WRITE_AT()
  return
end subroutine MPI_FILE_WRITE_AT

subroutine MPI_COMM_CREATE()
  return
end subroutine MPI_COMM_CREATE

subroutine MPI_COMM_GROUP()
  return
end subroutine MPI_COMM_GROUP

subroutine MPI_GROUP_EXCL()
  return
end subroutine MPI_GROUP_EXCL

subroutine MPI_GROUP_INCL()
  return
end subroutine MPI_GROUP_INCL

subroutine MPI_GROUP_RANK()
  return
end subroutine MPI_GROUP_RANK

subroutine MPI_GROUP_SIZE()
  return
end subroutine MPI_GROUP_SIZE


