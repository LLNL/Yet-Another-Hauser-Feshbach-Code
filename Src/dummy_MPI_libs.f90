!---------------------------------------------------------------------
!     *******  The following subroutines and functions are     *******
!     *******  to be used in the sequential environment.       *******
!     *******  and are "dummy" subroutines for MPI             *******
!---------------------------------------------------------------------
!
!*******************************************************************************
!
subroutine MPI_INIT(ierr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine for MPI initialization
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  implicit none
  integer(kind=4) :: ierr
  ierr = 0
  return
end subroutine MPI_INIT
!
!*******************************************************************************
!
subroutine MPI_COMM_RANK(icomm, iproc, ierr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subrtoutine to return the MPI rank
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  implicit none
  integer(kind=4) :: icomm, iproc, ierr
!-----------------------------------------
  icomm = 0
  iproc = 0
  ierr = 0
  return
end
!
!*******************************************************************************
!
subroutine MPI_COMM_SIZE(icomm, nproc, ierr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine to return the MPI size
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  implicit none
  integer(kind=4) :: icomm, nproc, ierr
!-----------------------------------------
  icomm = 0
  nproc = 1
  ierr = 0
  return
end
!
!*******************************************************************************
!
subroutine MPI_Barrier(icomm,ierr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine for setting an MPI barrier across all nodes
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  implicit none
  integer(kind=4) :: icomm, ierr
  icomm = 0
  ierr = 0
!-----------------------------------------
  return
end
!
!*******************************************************************************
!
subroutine MPI_Finalize(ierr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine to finalize an MPI run
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  implicit none
  integer(kind=4) :: ierr
!-----------------------------------------
  ierr = 0
  return
end
!
!*******************************************************************************
!
subroutine MPI_Bcast()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine for an MPI broadcast
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_Bcast
!
!*******************************************************************************
!
subroutine MPI_Reduce()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine for an MPI reduce
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_Reduce
!
!*******************************************************************************
!
subroutine MPI_Allreduce()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine to perform an all reduce over all MPI nodes
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
 return
end subroutine MPI_Allreduce
!
!*******************************************************************************
!
subroutine MPI_Abort()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a substitute subroutine for MPI_Abort that will terminate execution
!    when called 
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  call exit
  return
end subroutine MPI_Abort
!
!*******************************************************************************
!
subroutine MPI_INFO_CREATE()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy MPI subroutine
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_INFO_CREATE
!
!*******************************************************************************
!
subroutine MPI_INFO_SET()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy MPI subroutine
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_INFO_SET
!
!*******************************************************************************
!
subroutine MPI_FILE_OPEN()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy MPI subroutine for opening a file
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_FILE_OPEN
!
!*******************************************************************************
!
subroutine MPI_FILE_SET_VIEW()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy MPI subroutine that sets the view for the file
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_FILE_SET_VIEW
!
!*******************************************************************************
!
subroutine MPI_FILE_READ_AT()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy MPI subroutine to read at a specific point
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_FILE_READ_AT
!
!*******************************************************************************
!
subroutine MPI_FILE_WRITE_AT()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy MPI subroutine to write a file at a specific point
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_FILE_WRITE_AT
!
!*******************************************************************************
!
subroutine MPI_COMM_CREATE()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine to mimic creating an MPI communication group
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_COMM_CREATE
!
!*******************************************************************************
!
subroutine MPI_COMM_GROUP()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine 
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_COMM_GROUP
!
!*******************************************************************************
!
subroutine MPI_GROUP_EXCL()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine 
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_GROUP_EXCL
!
!*******************************************************************************
!
subroutine MPI_GROUP_INCL()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine 
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_GROUP_INCL
!
!*******************************************************************************
!
subroutine MPI_GROUP_RANK()
!
!*******************************************************************************
!
!  Discussion:
!
!    This function calculates well depth
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_GROUP_RANK
!
!*******************************************************************************
!
subroutine MPI_GROUP_SIZE()
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a dummy subroutine returning the group size
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
  return
end subroutine MPI_GROUP_SIZE


