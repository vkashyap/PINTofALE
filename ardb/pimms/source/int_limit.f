*+INT_LIMIT
        subroutine INT_LIMIT( special, results, n_res )

        implicit none

        real frac
        parameter( frac = 0.02 )

        integer special, n_res
        real results( 0: n_res )

*       Description:
*         Given a point source count rate calculates the exposure time
*         for 5 sigma detection assuming a background rate.
*
*       Arguments:
*         special      (i)   : flag for instrument 
*                               801 IBIS/ISGRI
*                               802 IBIS/PICSIT
*                               811 SPI
*                              but we only expect 801 now
*         results      (i)   : count rate of the point source 
*
*       Origin:
*
*       Attempted an update in 2008 March, but this "special"
*         is currently turned off due to the complexity of ISGRI
*
*-INT_LIMIT

        real bgd, temp, threesigma, work, syst, threesig_b, res1
        integer k
        character*64 pw_strng

        real brate( 3 )
        character*14 bdesc( 3 )
        data brate / 617.5, 325.1, 292.4 /
        data bdesc / '18 keV-1 MeV', '18 keV-100 keV', '100 keV-1 MeV' /

        if( results( 0 ) .gt. 1.0e15 ) then
          call PWRITE( 'Predicted rate very high: mistake maybe?' )
          return
        else if( special .ne. 801 ) then
          call PWRITE( 'ERROR:: Instrument code unrecognized' )
          return
        else if( n_res .ne. 3 ) then
          call PWRITE( 'ERROR:: Mismatch in number of sub-bands' )
          return
        else
          do k = 1, n_res
            res1 = results( k )
            bgd = brate( k )
            pw_strng = '%%% in the ' // bdesc( k ) // ' band'
            call PWRITE( pw_strng )
            if( res1 .gt. 5000.0 .or. res1 .lt. 0.01 ) then
              write( pw_strng, 111 ) res1, bgd
 111          format( '  ', 1p, e9.2, ' source c/s +',
     &                              0p, f6.1, ' background c/s' )
            else
              write( pw_strng, 112 ) res1, bgd
 112          format( '  ', f9.3, ' source c/s +',
     &                              f6.1, ' background c/s' )
            end if
            call PWRITE( pw_strng )
            temp = res1 * res1
            if( temp .gt. 1e-30 ) then
              work = 9.0 * ( res1 + bgd )
              threesigma = work / temp
              pw_strng =
     &            '     3-sigma detection will be achieved in         s'
              if( threesigma .gt. 5000.0 ) then
                write( pw_strng( 43: 51 ), '(1p,e9.2)' ) threesigma
              else if( threesigma .gt. 0.01 ) then
                write( pw_strng( 43: 51 ), '(f9.3)' ) threesigma
              else
                write( pw_strng( 43: 51 ), '(1p,e9.2)' ) threesigma
              end if
              call PWRITE( pw_strng )
              syst = frac * bgd
              if( res1 .gt. 3.0 * syst ) then
                threesig_b = work / ( temp - 9.0 * syst * syst )
                pw_strng =
     &  '     (or in         s with 2% systematic uncertainties in bgd)'
                if( threesig_b .gt. 5000.0 ) then
                  write( pw_strng( 12: 20 ), '(1p,e9.2)' ) threesig_b
                else if( threesig_b .gt. 0.01 ) then
                  write( pw_strng( 12: 20 ), '(f9.3)' ) threesig_b
                else
                  write( pw_strng( 12: 20 ), '(1p,e9.2)' ) threesig_b
                end if
              else
                pw_strng =
     & '     (but undetectable with 2% systematic uncertainties in bgd)'
              end if
            else
              pw_strng
     &           = '     Count rate too low: 3-sigma detection unlikely'
            end if
            call PWRITE( pw_strng )
          end do
        end if

        end
