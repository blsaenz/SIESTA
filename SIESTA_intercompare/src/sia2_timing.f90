! file: SIA2_timing.F90
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================
! ======================================================================


! ======================================================================
! Subroutine: sia2_getmonth
! Purpose: what numerical (1-12) month are we in? mo is returned
! ======================================================================
      SUBROUTINE sia2_getmonth(hour,leap_year,mo)
      double precision :: hour
      integer :: mo
      logical :: leap_year

      if (hour .lt. 744.) then
          mo = 1
          return
      elseif (leap_year) then
          if (hour .lt. 1440.) then
              mo = 2
              return
          elseif (hour .lt. 2184) then
              mo = 3
              return
          elseif (hour .lt. 2904) then
              mo = 4
              return
          elseif (hour .lt. 3648) then
              mo = 5
              return
          elseif (hour .lt. 4368) then
              mo = 6
              return
          elseif (hour .lt. 5112) then
              mo = 7
              return
          elseif (hour .lt. 5856) then
              mo = 8
              return
          elseif (hour .lt. 6576) then
              mo = 9
              return
          elseif (hour .lt. 7320) then
              mo = 10
              return
          elseif (hour .lt. 8040) then
              mo = 11
              return
          else
              mo = 12
              return
          endif
      else
         if (hour .lt. 1416.) then
              mo = 2
              return
          elseif (hour .lt. 2160) then
              mo = 3
              return
          elseif (hour .lt. 2880) then
              mo = 4
              return
          elseif (hour .lt. 3624) then
              mo = 5
              return
          elseif (hour .lt. 4344) then
              mo = 6
              return
          elseif (hour .lt. 5088) then
              mo = 7
              return
          elseif (hour .lt. 5832) then
              mo = 8
              return
          elseif (hour .lt. 6552) then
              mo = 9
              return
          elseif (hour .lt. 7296) then
              mo = 10
              return
          elseif (hour .lt. 8016) then
              mo = 11
              return
          else
              mo = 12
              return
          endif
      endif

      end subroutine sia2_getmonth



! ======================================================================
! Subroutine: sia2_getmonth
! Purpose: what numerical (1-12) month are we in? mo is returned
! ======================================================================
      SUBROUTINE sia2_getmonth_dble(hour,leap_year,mo)
      double precision :: hour,mo
      logical :: leap_year

      if (hour .lt. 744.) then
          mo = 1. + hour/744.
          return
      elseif (leap_year) then
          if (hour .lt. 1440.) then
              mo = 2. + (hour-744.)/(1440.-744.)
              return
          elseif (hour .lt. 2184.) then
              mo = 3. + (hour-1440.)/(2184.-1440.)
              return
          elseif (hour .lt. 2904.) then
              mo = 4. + (hour-2184.)/(2904.-2184.)
              return
          elseif (hour .lt. 3648.) then
              mo = 5. + (hour-2904.)/(3648.-2904.)
              return
          elseif (hour .lt. 4368.) then
              mo = 6. + (hour-3648.)/(4368.-3648.)
              return
          elseif (hour .lt. 5112.) then
              mo = 7. + (hour-4368.)/(5112.-4368.)
              return
          elseif (hour .lt. 5856.) then
              mo = 8. + (hour-5112.)/(5856.-5112.)
              return
          elseif (hour .lt. 6576.) then
              mo = 9. + (hour-5856.)/(6576.-5856.)
              return
          elseif (hour .lt. 7320.) then
              mo = 10. + (hour-6576.)/(7320.-6576.)
              return
          elseif (hour .lt. 8040.) then
              mo = 11. + (hour-7320.)/(8040.-7320.)
              return
          else
              mo = 12. + (hour-8040.)/(8784.-8040.)
              return
          endif
      else
         if (hour .lt. 1416.) then
              mo = 2. + (hour-744.)/(1416.-744.)
              return
          elseif (hour .lt. 2160) then
              mo = 3. + (hour-1416.)/(2160.-1416.)
              return
          elseif (hour .lt. 2880) then
              mo = 4. + (hour-2160.)/(2880.-2160.)
              return
          elseif (hour .lt. 3624) then
              mo = 5. + (hour-2880.)/(3624.-2880.)
              return
          elseif (hour .lt. 4344) then
              mo = 6. + (hour-3624.)/(4344.-3624.)
              return
          elseif (hour .lt. 5088) then
              mo = 7. + (hour-4344.)/(5088.-4344.)
              return
          elseif (hour .lt. 5832) then
              mo = 8. + (hour-5088.)/(5832.-5088.)
              return
          elseif (hour .lt. 6552) then
              mo = 9. + (hour-5832.)/(6552.-5832.)
              return
          elseif (hour .lt. 7296) then
              mo = 10. + (hour-6552.)/(7296.-6552.)
              return
          elseif (hour .lt. 8016) then
              mo = 11. + (hour-7296.)/(8016.-7296.)
              return
          else
              mo = 12. + (hour-8016.)/(8760.-8016.)
              return
          endif
      endif

      end subroutine sia2_getmonth_dble



! ======================================================================
! Subroutine: sia2_isleapyear
! Purpose: is it a leap year?
! ======================================================================
      SUBROUTINE sia2_isleapyear(year,leap_year)
      integer :: year
      logical :: leap_year

      if ((year .eq. 1980) .or. (year .eq. 1984) .or. (year .eq. &
             1988) .or. (year .eq. 1992) .or. (year .eq. 1996) .or. (year &
             .eq. 2000) .or. (year .eq. 2004)) then

          leap_year = .true.
      else
          leap_year = .false.
      endif

      end subroutine sia2_isleapyear


! ======================================================================
! Subroutine: sia2_end_of_year
! Purpose: decides when it is the end of the year, for incrementing year
! ======================================================================
      SUBROUTINE sia2_end_of_year(year,hour,eoy,leap)

      integer :: year
      double precision :: hour
      logical :: eoy,leap

      if (leap) then
          if (hour/24. .ge. 366.) then
              eoy=.true.
          else
              eoy=.false.
          endif
      else
          if (hour/24. .ge. 365.) then
              eoy=.true.
          else
              eoy=.false.
          endif
      endif

      end subroutine sia2_end_of_year


! ======================================================================
! Subroutine: sia2_num_steps
! Purpose: calculate the number of total model time steps
! ======================================================================
      SUBROUTINE sia2_num_steps(y1,jd1,y2,jd2,dt,steps)

      integer :: y1,y2,jd1,jd2,steps,i
      double precision :: dt, temp_steps

      temp_steps=0

      if(y2-y1 .gt. 0) then
          do i=y1,y2-1
              temp_steps = temp_steps + 365*24/dt
              if ((i .eq. 1980) .or. (i .eq. 1984) .or. (i .eq. &
                      1988) .or. (i .eq. 1992) .or. (i .eq. 1996) .or. (i &
                      .eq. 2000) .or. (i .eq. 2004)) then
                  temp_steps = temp_steps + 24/dt
              endif
          enddo
      endif
      temp_steps = temp_steps + jd2*24/dt - jd1*24/dt
      steps = int(temp_steps)

      end subroutine sia2_num_steps


! ======================================================================
! Subroutine: sia2_increment_year
! Purpose: sets up timing variables for a new year
! ======================================================================
      SUBROUTINE sia2_increment_year(year, hour, leap)

      integer :: year
      double precision :: hour
      logical :: leap

      hour = hour - 365*24

      if (leap) then
          hour=hour-24
      endif
      year = year+1

      ! set leap year global for new year
      call sia2_isleapyear(year,leap)

      end subroutine sia2_increment_year


! ======================================================================
! Subroutine: sia2_days_1980
! Purpose: finds date80 give the year and hour of year (date80 = internal
! timing format measured in hours since 1980)
! ======================================================================
      SUBROUTINE sia2_days_1980(hour,year,date80)

      integer :: j_date,year,date80
      double precision :: hour

      j_date = int(hour/24)+1  ! possibly just use 'int' instead of ifix?

      !print *, 'timing:j_date: ',j_date

      date80 = (year-1980)*365+j_date

       if (date80 .gt. 60) then
	      date80=date80+1
          if (date80 .gt. 1460) then
	          date80=date80+1
              if (date80 .gt. 2921) then
	              date80=date80+1
                  if (date80 .gt. 4382) then
	                  date80=date80+1
                      if (date80 .gt. 5843) then
	                      date80=date80+1
                          if (date80 .gt. 7304) then
	                          date80=date80+1
                              if (date80 .gt. 8765) then
	                              date80=date80+1
                                  if (date80 .gt. 10226) then
	                                  date80=date80+1
	                              endif
	                          endif
	                      endif
	                  endif
	              endif
	          endif
	      endif
      endif

      if (j_date .eq. 366) then
          date80=date80-1
      endif

      end SUBROUTINE sia2_days_1980


! ======================================================================
! Subroutine: sia2_date_from_days_1980
! Purpose: converts date80 to a common format for easy creation of filenames
! ======================================================================
     SUBROUTINE sia2_date_from_days_1980(date80,yearmonthday)

      integer :: thisdate,thisyear,date80,internal_date80
      character(8) :: yearmonthday

      !print *, 'timing:date_80: ', date80

      internal_date80 = date80
      if (internal_date80 .gt. 8826) then
          internal_date80=internal_date80-1
      endif
      if (internal_date80 .gt. 7365) then
          internal_date80=internal_date80-1
      endif
      if (internal_date80 .gt. 5094) then
          internal_date80=internal_date80-1
      endif
      if (internal_date80 .gt. 4443) then
          internal_date80=internal_date80-1
      endif
      if (internal_date80 .gt. 2982) then
          internal_date80=internal_date80-1
      endif
      if (internal_date80 .gt. 1521) then
          internal_date80=internal_date80-1
      endif
      if (internal_date80 .gt. 60) then
          internal_date80=internal_date80-1
      endif

      thisyear = internal_date80/365+1980

      !print *, 'timing:thisyear:', thisyear
	  internal_date80 = internal_date80 - (thisyear-1980)*365
      !correct for case of julian day 365
      if (internal_date80 .eq. 0) then
          thisyear = thisyear-1
          thisdate=1231
      else

		  if ((thisyear .eq. 1980) .or. (thisyear .eq. 1984) .or. (thisyear .eq. &
		  1988) .or. (thisyear .eq. 1992) .or. (thisyear .eq. 1996) .or. (thisyear &
		  .eq. 2000) .or. (thisyear .eq. 2004)) then
			  !print *, 'timing:internal_date80 (leapyr): ', internal_date80
			  if(internal_date80 .lt. 32) then
				  thisdate = 100+internal_date80
			  else if (internal_date80 .lt. 61) then
					thisdate = 200+internal_date80-31
			  else if (internal_date80 .lt. 92) then
					thisdate = 300+internal_date80-60
			  else if (internal_date80 .lt. 122) then
					thisdate = 400+internal_date80-91
			  else if (internal_date80 .lt. 153) then
					thisdate = 500+internal_date80-121
			  else if (internal_date80 .lt. 183) then
					thisdate = 600+internal_date80-152
			  else if (internal_date80 .lt. 214) then
					thisdate = 700+internal_date80-182
			  else if (internal_date80 .lt. 245) then
					thisdate = 800+internal_date80-213
			  else if (internal_date80 .lt. 275) then
					thisdate = 900+internal_date80-244
			  else if (internal_date80 .lt. 306) then
					thisdate = 1000+internal_date80-274
			  else if (internal_date80 .lt. 336) then
					thisdate = 1100+internal_date80-305
			  else if (internal_date80 .lt. 367) then
					thisdate = 1200+internal_date80-335
			  else
				  print *, 'err in date_from_days_1980: Julian date greater than 366'
				  call exit(1003)
			  endif
		  else
			  !print *, 'internal_date80: ', internal_date80
			  if(internal_date80 .lt. 32) then
				  thisdate = 100+internal_date80
			  else if (internal_date80 .lt. 60) then
					thisdate = 200+internal_date80-31
			  else if (internal_date80 .lt. 91) then
					thisdate = 300+internal_date80-59
			  else if (internal_date80 .lt. 121) then
					thisdate = 400+internal_date80-90
			  else if (internal_date80 .lt. 152) then
					thisdate = 500+internal_date80-120
			  else if (internal_date80 .lt. 182) then
					thisdate = 600+internal_date80-151
			  else if (internal_date80 .lt. 213) then
					thisdate = 700+internal_date80-181
			  else if (internal_date80 .lt. 244) then
					thisdate = 800+internal_date80-212
			  else if (internal_date80 .lt. 274) then
					thisdate = 900+internal_date80-243
			  else if (internal_date80 .lt. 305) then
					thisdate = 1000+internal_date80-273
			  else if (internal_date80 .lt. 335) then
					thisdate = 1100+internal_date80-304
			  else if (internal_date80 .lt. 366) then
					thisdate = 1200+internal_date80-334
			  else
				  print *, 'err in date_from_days_1980: Julian date greater than 365'
				  call exit(1004)
			  endif
		  endif
      endif
      write(yearmonthday(1:4),'(I4)') thisyear
      write(yearmonthday(5:8),'(I4.4)') thisdate

      end SUBROUTINE sia2_date_from_days_1980

