!
!*******************************************************************************
!
subroutine read_saved_parameters(data_path,len_path,icomp)                
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine to set up read previously saved nuclear parameters for 
!    nucleus icomp in the calculation. The parameters are written to the file
!    data_path(1:len_path)//'Saved-Parameters.txt' and under the header
!    start Z A. Parameters are read from this file using the parse_command
!    subroutine. The saved parameters may be structure or fission parameters
!    and are entered as standard YAHFC commands. Only nucleus specific parameters
!    can be entered in this manner. All global commands must be entered 
!    in the script file being executed.
! 
!    The saved parameters must be in a block between different saved
!    nuclei. The blocks are defined between different start elements. The file
!    may contain comments using the has '#' or bang '!'.
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        nodeinfo
!        constants
!
!     Subroutines:
!
!       parse_string
!       lower_case_word
!       extract_ZA_data
!       parse_command
!       exit_YAHFC
!
!     External functions:
!
!       integer(kind=4) :: rank_commands
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
!    28 June 2021
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   use variable_kinds
   use options
   use nuclei
   use nodeinfo
   use constants
   implicit none
   character(len=200), intent(in) :: data_path
   integer(kind=4), intent(in) :: len_path
   integer(kind=4), intent(in) :: icomp
!-------------------------------------------------------------------
   logical :: finish, finished, reading
   character(len=132) :: char_temp
   character(len=132) :: command_buffer
   character(len=132) :: command(100)    !   maximum of 100 commands - shouldn't be more
   integer(kind=4) :: command_rank(100)  !   being a bit lazy not wanting to run twice to count
   integer(kind=4) :: num_commands
   integer(kind=4) :: numw
   integer(kind=4) :: startw(66), stopw(66)
   integer(kind=4) :: read_err
   integer(kind=4) :: nchar
   integer(kind=4) :: itemp
   integer(kind=4) :: i, j

   logical :: read_error
   integer(kind=4) :: ndat
   integer(kind=4) :: nw
   real(kind=8) :: X
   integer(kind=4) :: izz, iaa
!------------------------------------------------------------------------------
   integer(kind=4) :: rank_commands
!------------------------------------------------------------------------------
!------   Check the file 'Saved-Fission-Parameters.txt'
   open(unit=153,file=data_path(1:len_path)//'Saved-Parameters.txt',status='old')
   finished = .false.
   reading = .false.
   num_commands = 0
   do while(.not. finished)
      read(153,'(a)', iostat = read_err)command_buffer
      if(read_err /= 0)finished = .true.                                     !  Hit end of file
      if(command_buffer(1:1) == '!' .or. command_buffer(1:1) == '#')cycle                  !  it is a comment line
      call parse_string(command_buffer,numw,startw,stopw)
      nchar = stopw(1)

      call lower_case_word(nchar,command_buffer(1:stopw(1)))                     !  convert to lower case
      if(command_buffer(startw(1):stopw(1)) == 'end')finished = .true.
      if(.not. reading)then
         if(command_buffer(startw(1):stopw(1)) == 'start')then
            ndat = 0
            call extract_ZA_data(command_buffer, numw, startw, stopw, ndat,         &
                                 iZZ, iAA, X, nw, read_error)

            if(read_error)then
               if(iproc == 0)then
                   write(6,*)'Error specifying nucleus in "Saved-Parameters.txt"'
                   write(6,*)'Either isotope or Z and A are not set correctly'
                   write(6,*)'Unable to determine how to read data. Fix, and restart'
                end if
                call exit_YAHFC(501)
            end if
            if(izz == nucleus(icomp)%Z .and. iaa == nucleus(icomp)%A)reading = .true.
         end if
      else
         if(command_buffer(startw(1):stopw(1)) == 'start')then
            exit                      !  hit the next nucleus, so we are done
         end if
         num_commands = num_commands + 1
         if(num_commands > 100)then
            if(iproc == 0)then
               write(6,*)'Too many commands in saved file > 100'
               write(6,*)'Increase dimension of command and command_rank in subroutine'
               write(6,*)'read_saved_parameters'
            end if
            call exit_YAHFC(501)
         end if
         command(num_commands) = command_buffer
      end if  
   end do
   close(unit=153)

   if(num_commands == 0)return
!---   Rank them so that they will be executed in the correct order
   do i = 1, num_commands
      command_rank(i) = rank_commands(command(i))
   end do  
!---------------------------------------------------------------------------+
!-------   order the commands by rank                                       +
!---------------------------------------------------------------------------+

   do i = 1, num_commands - 1
      do j = i + 1, num_commands
         if(command_rank(j) < command_rank(i))then
            itemp = command_rank(i)
            char_temp = command(i)
            command_rank(i) = command_rank(j)
            command(i) = command(j)
            command_rank(j) = itemp
            command(j) = char_temp
         end if
      end do
   end do

   if(print_me)write(6,*)'Reading and executing saved parameters from $YAHFC_DATA/Saved-Parameters.txt'

!---------------------------------------------------------------------------+
!-------   Now execute them using parse_commad                              +
!---------------------------------------------------------------------------+
   nucleus(icomp)%reading_param = .true.
   do i = 1, num_commands
      if(print_me)write(6,*)'executing command ',command(i)
      call parse_command(command(i),finish)           !  parse data and put in proper place
   end do
   nucleus(icomp)%reading_param = .false.
   nucleus(icomp)%param_read = .true.

   return

end subroutine read_saved_parameters              
