        subroutine SXS_GRADES( total )

        implicit none

        real psp_limit
        parameter( psp_limit = 630.0 )

        real total

*       Uses distribution of counts in the 36 pixels for a point source exactly
*       centered on the array, as provided by Takashi Okajima, integrates the
*       total number of Hi-res, mid-res primary etc. events using the pixel
*       by pixel formula as implimented in SXS_PERPIX.
*
*       KM, 2011 Apr

        real hi, mp, ms, lo
        real prate, phi, pmp, pms, plo
        integer k
        character*80 out_line

        real pix_sum, pix_fra( 36 )
        data pix_sum, pix_fra / 0.739681, 0.001659, 0.003808, 0.004764,
     &      0.004305, 0.003423, 0.001523, 0.003907, 0.012234, 0.022980,
     &      0.022729, 0.012713, 0.003543, 0.004543, 0.022722, 0.108968,
     &      0.106008, 0.023337, 0.004579, 0.004262, 0.022561, 0.112836,
     &      0.108747, 0.021570, 0.004784, 0.003923, 0.011928, 0.022291,
     &      0.021831, 0.012361, 0.004038, 0.001760, 0.003858, 0.004473,
     &                                    0.004916, 0.003964, 0.001834 /

        hi = 0.0
        mp = 0.0
        ms = 0.0
        lo = 0.0

        do k = 1, 36
          prate = total * pix_fra( k ) / pix_sum
          call SXS_PERPIX( prate, phi, pmp, pms, plo )
          hi = hi + phi * prate
          mp = mp + pmp * prate
          ms = ms + pms * prate
          lo = lo + plo * prate
        end do
        call PWRITE( ' ' )
        if( total .gt. psp_limit ) then
          out_line =
     &     '% WARNING: incoming count rate is higher than the PSP limit'
          call PWRITE( out_line )
        end if
        out_line = '% Count rates by event grade are estimated to be:'
        call PWRITE( out_line )
        out_line = '   high resolution   mid-res primary  ' //
     &                               'mid-res secondary  low resolution'
        call PWRITE( out_line )
        out_line = ' '
        if( hi .lt. 0.1 .or. hi .ge. 10000.0 ) then
          write( out_line( 7: 15 ), '(1p,e9.3)' ) hi
        else
          write( out_line( 7: 15 ), '(f9.4)' ) hi
        end if
        if( mp .lt. 0.1 .or. mp .ge. 10000.0 ) then
          write( out_line( 25: 33 ), '(1p,e9.3)' ) mp
        else
          write( out_line( 25: 33 ), '(f9.4)' ) mp
        end if
        if( ms .lt. 0.1 .or. ms .ge. 10000.0 ) then
          write( out_line( 43: 51 ), '(1p,e9.3)' ) ms
        else
          write( out_line( 43: 51 ), '(f9.4)' ) ms
        end if
        if( lo .lt. 0.1 .or. lo .ge. 10000.0 ) then
          write( out_line( 61: 69 ), '(1p,e9.3)' ) lo
        else
          write( out_line( 61: 69 ), '(f9.4)' ) lo
        end if
        call PWRITE( out_line )

        end




        subroutine SXS_PERPIX( total, hi, mp, ms, lo )

        implicit none

*       An implimentation by KM of the grade branching ratio formula
*       supplied by Caroline Kilbourne, 2011 Apr
*       Given an input count rate per pixel (total), returns the
*       rates of high-res, mid-res primary, mid-res secondary, and low-res
*       events.

        integer n, m, n4, m4
        real tconst
        parameter( n = 1024, m = 150, n4 = 256, m4 = 37 )
*       parameter( n = 2048, m = 300, n4 = 512, m4 = 75 ) for XRS
        parameter( tconst = 8.0e-5 )

	real total, hi, mp, ms, lo

        real tau

        tau = 2 * ( n - m ) * tconst
        hi = exp( -total * tau )
        tau = ( n - m + n4 - m4 ) * tconst
        mp = exp( -total * tau ) - hi
        tau = 2 * ( n4 - m4 ) * tconst
        ms = exp( -total * tau ) - hi - mp
        lo = 1.0 - hi - mp - ms

        end
