program x4_plot
   implicit none
   character*132 :: header(246)
   character*132 :: data_header(49)
   character*3 :: temp_char
   integer(kind=4) :: max_sets
   parameter (max_sets = 210)
   real(kind=4) :: x(max_sets,50000), y(max_sets,50000), dx(max_sets,50000), dy(max_sets,50000)
   real(kind=4) :: xx
   logical :: file_test
   integer(kind=4) :: num_sets, num_data, ndat(max_sets), data_type(max_sets), temp_type
   integer(kind=4) :: symbol(0:max_sets), color(0:max_sets), pattern(0:max_sets)
   character*100 :: data_label(max_sets), temp_label
   character*132 input, output
   character*132 output_file
   character*132 input_file
   integer(kind=4) :: num_words, start_word(66), stop_word(66)
   integer(kind=4) :: i, j, k, l
   integer(kind=4) :: ilast, istart, istop, ibegin, iend, length
   integer(kind=4) :: num_files
   real(kind=4) :: xmin, xmax
   real(kind=4) :: xplot_min, xplot_max, yplot_min, yplot_max
!---------------   internal functions ----------------------
   integer(kind=4) :: last_char
!---------------    Start program   ------------------------
   call set_header(header)
   call set_data_header(data_header)
   read(5,'(a)')output_file
   ilast = index(output_file,' ') - 1
   open(unit=20, file = output_file(1:ilast)//'.agr', status = 'unknown')
   read(5,*)num_files
   if(num_files > max_sets) stop 'Too many files, cannot process'
   read(5,*)xmin
   read(5,*)xmax
   if(xmax < 0.0) xmax = 1.0e10
   num_sets = 0
   symbol(0) = 1
   color(0) = 0
   pattern(0) = 0
   do i = 1, num_files
      input(1:132)= ' '
      read(5,'(a)')input
      call get_file_name(input,input_file,length)
      open(unit=10, file = input_file(1:length), status = 'old')
 2    read(10,'(a)')input
      temp_type = 0
      if(input(1:1) == '#')then
         call parse_string(input, num_words, start_word, stop_word)    !   turn string into num_words words
         if(num_words == 1) goto 2
         if(num_words >= 2 .and. input(start_word(2):stop_word(2)) == 'SubEntry:')then
            temp_label = input(start_word(3):stop_word(3))
         end if
         if(num_words >= 2 .and. input(start_word(2):stop_word(2)) == 'Energy')then
            if(num_words == 4)then
               if(input(start_word(4):stop_word(4)) == 'd(Data)') temp_type = 1
               if(input(start_word(4):stop_word(4)) == 'd(Energy)') temp_type = 3
            elseif(num_words == 5)then
               if(input(start_word(4):stop_word(4)) == 'd(Energy)' .and.          &
                  input(start_word(5):stop_word(5)) == 'd(Data)') temp_type = 3
            end if
            input(1:132)=' '
            read(10,'(a)')input
            call parse_string(input, num_words, start_word, stop_word)    !   turn string into num_words words
            if(input(start_word(3):stop_word(3)) /= 'barns')cycle      !   This is not a cross section file
         else
            goto 2
         end if
      end if
!-------------    scan Energy (x) column to see if data falls within min and max, otherwise throw out
      file_test = .false.
  3   read(10,*,end=4)xx
      if(xx >= xmin .and. xx <= xmax)file_test = .true.
!  write(6,*)xx
      goto 3
  4   continue
      if(.not. file_test)cycle       !    data doesn't fall within limits so, skip
      rewind(10)
      num_sets = num_sets + 1
      data_type(num_sets) = temp_type
      num_data = 0
  5   read(10,'(a)',end=6)input
      if(input(1:1) /= '#')then     !   this is a data entry
         call parse_string(input, num_words, start_word, stop_word)    !   turn string into num_words words
         if(input(start_word(1):stop_word(1)) == 'None')then
            goto 5
         else
            read(input(start_word(1):stop_word(1)),*)xx
         end if
         if(xx >= xmin .and. xx <= xmax)then
            num_data = num_data + 1
            x(num_sets, num_data) = xx
            if(input(start_word(2):stop_word(2)) == 'None')then
               goto 5
            else
               read(input(start_word(2):stop_word(2)),*)y(num_sets, num_data)
            end if
            if(temp_type == 1)then
               dy(num_sets, num_data) = 0.0
               if(input(start_word(3):stop_word(3)) /= 'None')read(input(start_word(3):stop_word(3)),*)dy(num_sets, num_data)
            end if
            if(temp_type == 2)then
               dx(num_sets, num_data) = 0.0
               if(input(start_word(3):stop_word(3)) /= 'None')read(input(start_word(3):stop_word(3)),*)dx(num_sets, num_data)
               dy(num_sets, num_data) = 0.0
               if(input(start_word(4):stop_word(4)) /= 'None')read(input(start_word(4):stop_word(4)),*)dy(num_sets, num_data)
            end if
            if(temp_type == 3)then
               dx(num_sets, num_data) = 0.0
               if(input(start_word(3):stop_word(3)) /= 'None')read(input(start_word(3):stop_word(3)),*)dx(num_sets, num_data)
            end if
         end if
      end if     
      goto 5
  6   continue
      ndat(num_sets) = num_data
      symbol(num_sets) = symbol(num_sets - 1)
      color(num_sets) = color(num_sets - 1) + 1
      pattern(num_sets) = pattern(num_sets - 1)
      data_label(num_sets) = temp_label
      if(color(num_sets) > 15)then
         color(num_sets) = 1
         symbol(num_sets) = symbol(num_sets) + 1
         if(symbol(num_sets) > 7)then
            symbol(num_sets) = 1
            pattern(num_sets) = pattern(num_sets) + 1
         end if
      end if
      close(unit=10)
   end do

! ------    Scan data to get max x and max y to set world size
   xplot_min = 1.0e10
   xplot_max = 0.0
   yplot_min = 1.0e10
   yplot_max = 0.0
   do i = 1, num_sets
      do j = 1, ndat(i)
         if(x(i,j) > xplot_max) xplot_max = x(i,j)
         if(y(i,j) > yplot_max) yplot_max = y(i,j)
      end do
   end do
   write(6,*)xplot_max, yplot_max
   k = xplot_max/5.0
   xplot_max = real(k+1)*5.0
   k = yplot_max/0.5
   yplot_max = real(k+1)*0.5

!-----------    Set up Header information to make agr file
   do j = 1, 246
      write(output,'(a)')header(j)
      iend = last_char(output)
!-------    Special as it specifies the world
      if(j == 104)then
         write(output(iend+1:132),'(3(1x,f6.1,",")1x,f6.1)')0.0, 0.0, xplot_max, yplot_max
         iend = last_char(output)
      end if
      write(20,'(a)')output(1:iend)
   end do
!--------   Now the information for each data set
   do i = 1, num_sets
      do j = 1, 49
         output(1:6) = '@    s'
         if(i < 11)then
            write(output(7:7),'(i1)')i-1
            istart = 8
         elseif(i < 101)then
            write(output(7:8),'(i2)')i-1
            istart = 9
         elseif(i <= max_sets)then
            write(output(7:9),'(i3)')i-1
            istart = 10
         end if
         output = output(1:istart -1)//' '//data_header(j)
         iend = last_char(output)
         if(j == 2)then
            if(data_type(i) == 0)then
               output = output(1:iend)//' xy'
            end if
            if(data_type(i) == 1)then
               output = output(1:iend)//' xydy'
            end if
            if(data_type(i) == 2)then
               output = output(1:iend)//' xydxdy'
            end if
            if(data_type(i) == 3)then
               output = output(1:iend)//' xydx'
            end if
         end if
         iend = last_char(output)
         temp_char = '   '
         if(j == 3)then
            write(temp_char(1:2),'(1x,i1)')symbol(i)
            output = output(1:iend)//temp_char
         end if
         if(j == 5 .or. j == 7 .or. j == 17 .or. j == 24 .or. j ==39)then
            if(color(i) < 10)then
               write(temp_char(1:2),'(1x,i1)')color(i)
            else
               write(temp_char(1:3),'(1x,i2)')color(i)
            end if
            output = output(1:iend)//temp_char
         end if
         if(j == 8)then
            write(temp_char(1:2),'(1x,i1)')pattern(i)
            output = output(1:iend)//temp_char
         end if
!         if( j == 48)then
!            output = output(1:iend)//' "'//output_file(1:ilast)//'.agr'
!            iend = last_char(output)
!            iend = iend + 1
!            output(iend:iend) = '"'
!         end if
         if( j == 48)then
            output = output(1:iend)//' "'//data_label(i)
            iend = last_char(output)
            iend = iend + 1
            output(iend:iend) = '"'
         end if
         iend = last_char(output)
         write(20,'(a)')output(1:iend)
      end do
   end do
  

! -------     Now write out data for plotting
   do i = 1, num_sets
       if(i < 11)then
          write(20,'(''@target G0.S'',i1)')i-1
       elseif(i < 101)then
          write(20,'(''@target G0.S'',i2)')i-1
       elseif(i <= max_sets)then        
          write(20,'(''@target G0.S'',i3)')i-1
       end if
       if(data_type(i) == 0)write(20,'(''@type xy'')')
       if(data_type(i) == 1)write(20,'(''@type xydy'')')
       if(data_type(i) == 2)write(20,'(''@type xydxdy'')')
       if(data_type(i) == 3)write(20,'(''@type xydx'')')
       do j = 1, ndat(i)
          if(data_type(i) == 0)write(20,*)x(i,j),y(i,j)
          if(data_type(i) == 1)write(20,*)x(i,j),y(i,j),dy(i,j)
          if(data_type(i) == 2)write(20,*)x(i,j),y(i,j),dx(i,j),dy(i,j)
          if(data_type(i) == 3)write(20,*)x(i,j),y(i,j),dx(i,j)          
       end do
       write(20,'(''&'')')
   end do

end program x4_plot

integer(kind=4) function last_char(string)
   implicit none
   character*132 string
   integer(kind=4) :: i
   do i = 132, 1, -1
      if(string(i:i) /= ' ')then
         last_char = i
          exit
      end if
   end do
   return
end function last_char


subroutine set_header(header)
   character*132 :: header(246)
    header(  1) = '# Grace project file '
    header(  2) = '# '
    header(  3) = '@version 50125 '
    header(  4) = '@page size 792, 612 '
    header(  5) = '@page scroll 5% '
    header(  6) = '@page inout 5% '
    header(  7) = '@link page off '
    header(  8) = '@map font 0 to "Times-Roman", "Times-Roman" '
    header(  9) = '@map font 1 to "Times-Italic", "Times-Italic" '
    header( 10) = '@map font 2 to "Times-Bold", "Times-Bold" '
    header( 11) = '@map font 3 to "Times-BoldItalic", "Times-BoldItalic" '
    header( 12) = '@map font 4 to "Helvetica", "Helvetica" '
    header( 13) = '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique" '
    header( 14) = '@map font 6 to "Helvetica-Bold", "Helvetica-Bold" '
    header( 15) = '@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique" '
    header( 16) = '@map font 8 to "Courier", "Courier" '
    header( 17) = '@map font 9 to "Courier-Oblique", "Courier-Oblique" '
    header( 18) = '@map font 10 to "Courier-Bold", "Courier-Bold" '
    header( 19) = '@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique" '
    header( 20) = '@map font 12 to "Symbol", "Symbol" '
    header( 21) = '@map font 13 to "ZapfDingbats", "ZapfDingbats" '
    header( 22) = '@map color 0 to (255, 255, 255), "white" '
    header( 23) = '@map color 1 to (0, 0, 0), "black" '
    header( 24) = '@map color 2 to (255, 0, 0), "red" '
    header( 25) = '@map color 3 to (0, 255, 0), "green" '
    header( 26) = '@map color 4 to (0, 0, 255), "blue" '
    header( 27) = '@map color 5 to (255, 255, 0), "yellow" '
    header( 28) = '@map color 6 to (188, 143, 143), "brown" '
    header( 29) = '@map color 7 to (220, 220, 220), "grey" '
    header( 30) = '@map color 8 to (148, 0, 211), "violet" '
    header( 31) = '@map color 9 to (0, 255, 255), "cyan" '
    header( 32) = '@map color 10 to (255, 0, 255), "magenta" '
    header( 33) = '@map color 11 to (255, 165, 0), "orange" '
    header( 34) = '@map color 12 to (114, 33, 188), "indigo" '
    header( 35) = '@map color 13 to (103, 7, 72), "maroon" '
    header( 36) = '@map color 14 to (64, 224, 208), "turquoise" '
    header( 37) = '@map color 15 to (0, 139, 0), "green4" '
    header( 38) = '@reference date 0 '
    header( 39) = '@date wrap off '
    header( 40) = '@date wrap year 1950 '
    header( 41) = '@default linewidth 1.0 '
    header( 42) = '@default linestyle 1 '
    header( 43) = '@default color 1 '
    header( 44) = '@default pattern 1 '
    header( 45) = '@default font 0 '
    header( 46) = '@default char size 1.000000 '
    header( 47) = '@default symbol size 1.000000 '
    header( 48) = '@default sformat "%.8g" '
    header( 49) = '@background color 0 '
    header( 50) = '@page background fill on '
    header( 51) = '@timestamp off '
    header( 52) = '@timestamp 0.03, 0.03 '
    header( 53) = '@timestamp color 1 '
    header( 54) = '@timestamp rot 0 '
    header( 55) = '@timestamp font 0 '
    header( 56) = '@timestamp char size 1.000000 '
    header( 57) = '@timestamp def "Tue Feb  7 17:40:47 2017" '
    header( 58) = '@r0 off '
    header( 59) = '@link r0 to g0 '
    header( 60) = '@r0 type above '
    header( 61) = '@r0 linestyle 1 '
    header( 62) = '@r0 linewidth 1.0 '
    header( 63) = '@r0 color 1 '
    header( 64) = '@r0 line 0, 0, 0, 0 '
    header( 65) = '@r1 off '
    header( 66) = '@link r1 to g0 '
    header( 67) = '@r1 type above '
    header( 68) = '@r1 linestyle 1 '
    header( 69) = '@r1 linewidth 1.0 '
    header( 70) = '@r1 color 1 '
    header( 71) = '@r1 line 0, 0, 0, 0 '
    header( 72) = '@r2 off '
    header( 73) = '@link r2 to g0 '
    header( 74) = '@r2 type above '
    header( 75) = '@r2 linestyle 1 '
    header( 76) = '@r2 linewidth 1.0 '
    header( 77) = '@r2 color 1 '
    header( 78) = '@r2 line 0, 0, 0, 0 '
    header( 79) = '@r3 off '
    header( 80) = '@link r3 to g0 '
    header( 81) = '@r3 type above '
    header( 82) = '@r3 linestyle 1 '
    header( 83) = '@r3 linewidth 1.0 '
    header( 84) = '@r3 color 1 '
    header( 85) = '@r3 line 0, 0, 0, 0 '
    header( 86) = '@r4 off '
    header( 87) = '@link r4 to g0 '
    header( 88) = '@r4 type above '
    header( 89) = '@r4 linestyle 1 '
    header( 90) = '@r4 linewidth 1.0 '
    header( 91) = '@r4 color 1 '
    header( 92) = '@r4 line 0, 0, 0, 0 '
    header( 93) = '@g0 on '
    header( 94) = '@g0 hidden false '
    header( 95) = '@g0 type XY '
    header( 96) = '@g0 stacked false '
    header( 97) = '@g0 bar hgap 0.000000 '
    header( 98) = '@g0 fixedpoint off '
    header( 99) = '@g0 fixedpoint type 0 '
    header(100) = '@g0 fixedpoint xy 0.000000, 0.000000 '
    header(101) = '@g0 fixedpoint format general general '
    header(102) = '@g0 fixedpoint prec 6, 6 '
    header(103) = '@with g0 '
    header(104) = '@    world '
    header(105) = '@    stack world 0, 0, 0, 0 '
    header(106) = '@    znorm 1 '
    header(107) = '@    view 0.150000, 0.150000, 1.150000, 0.850000 '
    header(108) = '@    title "" '
    header(109) = '@    title font 0 '
    header(110) = '@    title size 1.500000 '
    header(111) = '@    title color 1 '
    header(112) = '@    subtitle "" '
    header(113) = '@    subtitle font 0 '
    header(114) = '@    subtitle size 1.000000 '
    header(115) = '@    subtitle color 1 '
    header(116) = '@    xaxes scale Normal '
    header(117) = '@    yaxes scale Normal '
    header(118) = '@    xaxes invert off '
    header(119) = '@    yaxes invert off '
    header(120) = '@    xaxis  on '
    header(121) = '@    xaxis  type zero false '
    header(122) = '@    xaxis  offset 0.000000 , 0.000000 '
    header(123) = '@    xaxis  bar on '
    header(124) = '@    xaxis  bar color 1 '
    header(125) = '@    xaxis  bar linestyle 1 '
    header(126) = '@    xaxis  bar linewidth 1.0 '
    header(127) = '@    xaxis  label "E\sn\N (MeV)" '
    header(128) = '@    xaxis  label layout para '
    header(129) = '@    xaxis  label place auto '
    header(130) = '@    xaxis  label char size 1.000000 '
    header(131) = '@    xaxis  label font 0 '
    header(132) = '@    xaxis  label color 1 '
    header(133) = '@    xaxis  label place normal '
    header(134) = '@    xaxis  tick on '
    header(135) = '@    xaxis  tick major 5.0 '
    header(136) = '@    xaxis  tick minor ticks 4 '
    header(137) = '@    xaxis  tick default 6 '
    header(138) = '@    xaxis  tick place rounded true '
    header(139) = '@    xaxis  tick in '
    header(140) = '@    xaxis  tick major size 1.000000 '
    header(141) = '@    xaxis  tick major color 1 '
    header(142) = '@    xaxis  tick major linewidth 1.0 '
    header(143) = '@    xaxis  tick major linestyle 1 '
    header(144) = '@    xaxis  tick major grid off '
    header(145) = '@    xaxis  tick minor color 1 '
    header(146) = '@    xaxis  tick minor linewidth 1.0 '
    header(147) = '@    xaxis  tick minor linestyle 1 '
    header(148) = '@    xaxis  tick minor grid off '
    header(149) = '@    xaxis  tick minor size 0.500000 '
    header(150) = '@    xaxis  ticklabel on '
    header(151) = '@    xaxis  ticklabel format general '
    header(152) = '@    xaxis  ticklabel prec 5 '
    header(153) = '@    xaxis  ticklabel formula "" '
    header(154) = '@    xaxis  ticklabel append "" '
    header(155) = '@    xaxis  ticklabel prepend "" '
    header(156) = '@    xaxis  ticklabel angle 0 '
    header(157) = '@    xaxis  ticklabel skip 0 '
    header(158) = '@    xaxis  ticklabel stagger 0 '
    header(159) = '@    xaxis  ticklabel place normal '
    header(160) = '@    xaxis  ticklabel offset auto '
    header(161) = '@    xaxis  ticklabel offset 0.000000 , 0.010000 '
    header(162) = '@    xaxis  ticklabel start type auto '
    header(163) = '@    xaxis  ticklabel start 0.000000 '
    header(164) = '@    xaxis  ticklabel stop type auto '
    header(165) = '@    xaxis  ticklabel stop 0.000000 '
    header(166) = '@    xaxis  ticklabel char size 1.000000 '
    header(167) = '@    xaxis  ticklabel font 0 '
    header(168) = '@    xaxis  ticklabel color 1 '
    header(169) = '@    xaxis  tick place both '
    header(170) = '@    xaxis  tick spec type none '
    header(171) = '@    yaxis  on '
    header(172) = '@    yaxis  type zero false '
    header(173) = '@    yaxis  offset 0.000000 , 0.000000 '
    header(174) = '@    yaxis  bar on '
    header(175) = '@    yaxis  bar color 1 '
    header(176) = '@    yaxis  bar linestyle 1 '
    header(177) = '@    yaxis  bar linewidth 1.0 '
    header(178) = '@    yaxis  label "\f{12}s\f{0} (b)" '
    header(179) = '@    yaxis  label layout para '
    header(180) = '@    yaxis  label place auto '
    header(181) = '@    yaxis  label char size 1.000000 '
    header(182) = '@    yaxis  label font 0 '
    header(183) = '@    yaxis  label color 1 '
    header(184) = '@    yaxis  label place normal '
    header(185) = '@    yaxis  tick on '
    header(186) = '@    yaxis  tick major  0.5 '
    header(187) = '@    yaxis  tick minor ticks 4 '
    header(188) = '@    yaxis  tick default 6 '
    header(189) = '@    yaxis  tick place rounded true '
    header(190) = '@    yaxis  tick in '
    header(191) = '@    yaxis  tick major size 1.000000 '
    header(192) = '@    yaxis  tick major color 1 '
    header(193) = '@    yaxis  tick major linewidth 1.0 '
    header(194) = '@    yaxis  tick major linestyle 1 '
    header(195) = '@    yaxis  tick major grid off '
    header(196) = '@    yaxis  tick minor color 1 '
    header(197) = '@    yaxis  tick minor linewidth 1.0 '
    header(198) = '@    yaxis  tick minor linestyle 1 '
    header(199) = '@    yaxis  tick minor grid off '
    header(200) = '@    yaxis  tick minor size 0.500000 '
    header(201) = '@    yaxis  ticklabel on '
    header(202) = '@    yaxis  ticklabel format general '
    header(203) = '@    yaxis  ticklabel prec 5 '
    header(204) = '@    yaxis  ticklabel formula "" '
    header(205) = '@    yaxis  ticklabel append "" '
    header(206) = '@    yaxis  ticklabel prepend "" '
    header(207) = '@    yaxis  ticklabel angle 0 '
    header(208) = '@    yaxis  ticklabel skip 0 '
    header(209) = '@    yaxis  ticklabel stagger 0 '
    header(210) = '@    yaxis  ticklabel place normal '
    header(211) = '@    yaxis  ticklabel offset auto '
    header(212) = '@    yaxis  ticklabel offset 0.000000 , 0.010000 '
    header(213) = '@    yaxis  ticklabel start type auto '
    header(214) = '@    yaxis  ticklabel start 0.000000 '
    header(215) = '@    yaxis  ticklabel stop type auto '
    header(216) = '@    yaxis  ticklabel stop 0.000000 '
    header(217) = '@    yaxis  ticklabel char size 1.000000 '
    header(218) = '@    yaxis  ticklabel font 0 '
    header(219) = '@    yaxis  ticklabel color 1 '
    header(220) = '@    yaxis  tick place both '
    header(221) = '@    yaxis  tick spec type none '
    header(222) = '@    altxaxis  off '
    header(223) = '@    altyaxis  off '
    header(224) = '@    legend on '
    header(225) = '@    legend loctype view '
    header(226) = '@    legend 0.85, 0.8 '
    header(227) = '@    legend box color 1 '
    header(228) = '@    legend box pattern 1 '
    header(229) = '@    legend box linewidth 1.0 '
    header(230) = '@    legend box linestyle 1 '
    header(231) = '@    legend box fill color 0 '
    header(232) = '@    legend box fill pattern 1 '
    header(233) = '@    legend font 0 '
    header(234) = '@    legend char size 1.000000 '
    header(235) = '@    legend color 1 '
    header(236) = '@    legend length 4 '
    header(237) = '@    legend vgap 1 '
    header(238) = '@    legend hgap 1 '
    header(239) = '@    legend invert false '
    header(240) = '@    frame type 0 '
    header(241) = '@    frame linestyle 1 '
    header(242) = '@    frame linewidth 1.0 '
    header(243) = '@    frame color 1 '
    header(244) = '@    frame pattern 1 '
    header(245) = '@    frame background color 0 '
    header(246) = '@    frame background pattern 0 '
    return
end subroutine set_header

subroutine set_data_header(data_header)
    character*132 data_header(49)
     data_header(1) = 'hidden false'
     data_header(2) = 'type'
     data_header(3) = 'symbol' 
     data_header(4) = 'symbol size 0.300000'
     data_header(5) = 'symbol color' 
     data_header(6) = 'symbol pattern 1' 
     data_header(7) = 'symbol fill color' 
     data_header(8) = 'symbol fill pattern'
     data_header(9) = 'symbol linewidth 1.0'
     data_header(10) = 'symbol linestyle 1'
     data_header(11) = 'symbol char 65'
     data_header(12) = 'symbol char font 0'
     data_header(13) = 'symbol skip 0'
     data_header(14) = 'line type 1'
     data_header(15) = 'line linestyle 0'
     data_header(16) = 'line linewidth 1.0'
     data_header(17) = 'line color'
     data_header(18) = 'line pattern 1'
     data_header(19) = 'baseline type 0'
     data_header(20) = 'baseline off'
     data_header(21) = 'dropline off'
     data_header(22) = 'fill type 0'
     data_header(23) = 'fill rule 0'
     data_header(24) = 'fill color'
     data_header(25) = 'fill pattern 1'
     data_header(26) = 'avalue off'
     data_header(27) = 'avalue type 2'
     data_header(28) = 'avalue char size 1.000000'
     data_header(29) = 'avalue font 0'
     data_header(30) = 'avalue color 1'
     data_header(31) = 'avalue rot 0'
     data_header(32) = 'avalue format general'
     data_header(33) = 'avalue prec 3'
     data_header(34) = 'avalue prepend ""'
     data_header(35) = 'avalue append ""'
     data_header(36) = 'avalue offset 0.000000 , 0.000000'
     data_header(37) = 'errorbar on'
     data_header(38) = 'errorbar place both'
     data_header(39) = 'errorbar color'
     data_header(40) = 'errorbar pattern 1'
     data_header(41) = 'errorbar size 0.000000'
     data_header(42) = 'errorbar linewidth 1.0'
     data_header(43) = 'errorbar linestyle 1'
     data_header(44) = 'errorbar riser linewidth 1.0'
     data_header(45) = 'errorbar riser linestyle 1'
     data_header(46) = 'errorbar riser clip off'
     data_header(47) = 'errorbar riser clip length 0.100000'
     data_header(48) = 'comment '
     data_header(49) = 'legend ""'
    return
end subroutine set_data_header

subroutine get_file_name(input,file,length)
   implicit none
   character*132 input, file
   integer(kind=4) :: length
   integer(kind=4) :: istart,iend, i
   integer(kind=4) :: last_char
   write(6,'(a)')input
   iend = last_char(input)
   do i = iend-1, 1, -1
      if(input(i:i) == ' ')then
         istart = i + 1
         exit
      end if
   end do
   file = input(istart:iend)
   length = last_char(file)
   return
end subroutine get_file_name

subroutine parse_string(input,numw,startw,stopw)
   implicit none
   character*132 input
   integer(kind=4) :: numw
   integer(kind=4) :: startw(66), stopw(66)
   integer i, j, k, n, istart, jstart
   logical word
   startw(1:66) = 0
   stopw(1:66) = 0
   numw = 0
   word = .false.
   k = 1
 1 if(input(k:k) /= ' ')then
      if(.not.word)then
         numw = numw + 1
         word = .true.
         startw(numw) = k
         if(k == 132)then
            stopw(numw) = k
            goto 4
         end if
      end if
   elseif(input(k:k) == ' ')then
      if(word)then
         word = .false.
         stopw(numw) = k - 1
      end if
      if(k == 132) goto 4
   end if
   k = k + 1
   if(k > 132)goto 4
   goto 1
 4 continue

!   write(6,'(a)')input
!   write(6,*)numw
!   do k = 1, numw
!!      write(6,*)startw(k),stopw(k)
!      write(6,'(a)')input(startw(k):stopw(k))
!   end do

   return
end subroutine parse_string

