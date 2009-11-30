 999     CALL SYSTEM_CLOCK(end_count)

         elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
c
c  Excluding assignmentgs
c
         if (elapsed_time .ne. 0.0) then

            write(*,'(/,a,f20.5,a,i3.3)') 'The elapsed time is',elapsed_time,
     &           ' seconds in configuration:',cfg

         endif

      end do                    !end cfg loop
