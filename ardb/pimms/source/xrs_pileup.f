        subroutine XRS_PILEUP( cps_given, special )

        implicit none

        real cps_given
	integer special

*       An implimentation, by KM, of the pile-up fraction formulae derived
*         from a series of xrssim simulation by Randall Smith.
*       "Special" is hardwired, must match the values in the *.special
*         files, expected to be 701 (OPEN/ND) or 702 (BE)
*         ++Note that the same formula is used regardless of source
*           spectrum or use of filter.  In reality, there should be
*           a dependence, due to the energy dependence of PSF++
*       Re-written for Astro-E2 on 2004 June 17 by KM
*       (Structure based on the ASTRO-E version written in 1999)

        real ndRate
        character*80 out_line

        call PWRITE( ' ' )
        out_line = 'Filter     Hi-Res  Hi+Med Res Med (S)  Low-Res'
        call PWRITE( out_line )
        if( special .eq. 701 ) then
          out_line = 'Open'
          call XRS_LINE( cps_given, out_line )
          call PWRITE( out_line )
          out_line = 'ND (10%)'
          ndRate = cps_given * 0.0953
          call XRS_LINE( ndRate, out_line )
          call PWRITE( out_line )
          out_line = '  (For a total of           ct/s'
     &                                     // ' with the 10% ND filter)'
          if( ndRate .lt. 0.1 .or. ndRate .ge. 10000.0 ) then
            write( out_line( 19: 27 ), '(1p,e9.3)' ) ndRate
          else
            write( out_line( 19: 27 ), '(f9.4)' ) ndRate
          end if
          call PWRITE( out_line )
        else
          out_line = 'Be'
          call XRS_LINE( cps_given, out_line )
          call PWRITE( out_line )
        end if

        end

        subroutine XRS_LINE( Rate, out_line )

        implicit none

        real Rate
        character*( * ) out_line

*       Wrapper to calculate the rate by filters and write it to a line
*         KM, 2004 June 17

        real hfrac, hmfrac, hmsfrac
        real hrate, hmrate, msrate, lrate

        real XRS_HFRAC, XRS_HMFRAC, XRS_HMSFRAC

        hfrac = XRS_HFRAC( Rate )
        hmfrac = XRS_HMFRAC( Rate )
        hmsfrac = XRS_HMSFRAC( Rate )

        hrate = Rate * hfrac
        hmrate = Rate * hmfrac
        msrate = Rate * ( hmsfrac - hmfrac )
        lrate = Rate * ( 1.0 - hmsfrac )

        if( hrate .lt. 0.1 .or. hrate .ge. 10000.0 ) then
          write( out_line( 10: 18 ), '(1p,e9.3)' ) hrate
        else
          write( out_line( 10: 18 ), '(f9.4)' ) hrate
        end if

        if( hmrate .lt. 0.1 .or. hmrate .ge. 10000.0 ) then
          write( out_line( 20: 28 ), '(1p,e9.3)' ) hmrate
        else
          write( out_line( 20: 28 ), '(f9.4)' ) hmrate
        end if

        if( msrate .lt. 0.1 .or. msrate .ge. 10000.0 ) then
          write( out_line( 30: 38 ), '(1p,e9.3)' ) msrate
        else
          write( out_line( 30: 38 ), '(f9.4)' ) msrate
        end if

        if( lrate .lt. 0.1 .or. lrate .ge. 10000.0 ) then
          write( out_line( 40: 48 ), '(1p,e9.3)' ) lrate
        else
          write( out_line( 40: 48 ), '(f9.4)' ) lrate
        end if

        if( Rate .ge. 160.0 ) then
          out_line( 50: 79 ) = 'Telemetry Saturation (160 c/s)'
        end if

        end



        real function XRS_HFRAC( c )

        implicit none

        real c

*       Fraction of Hi-res event, for a point source observed on-axis.
*         Tchebycheff polinomial fit to the ln of fraction as function
*         of count rate (* 0.001); linear extrapolation outside the
*         simulation range.
*       Written by KM, 2004 June 15

        double precision c3, temp
        integer k
        double precision coeffs( 0: 5 )
        data coeffs / 0.0025886842, -15.719536, 53.5910033,
     &                -140.072003, 181.983732, -88.5108818 /

        if( c .le. 1.3 ) then
          XRS_HFRAC = 1.0 - c / 1.3 * 0.017599761
        else if( c .le. 822.5821 ) then
          temp = 0.0d+00
          c3 = dble( c ) * 1.0d-03
          do k = 5, 0, -1
            temp = temp * c3 + coeffs( k )
          end do
          XRS_HFRAC = exp( temp )
        else
          c3 = ( c - 551.37 ) / 271.2121
c         log( 0.028907 ) = -3.54367
c         log( 0.0096235 ) - log( 0.028907 ) = -1.09988
          temp = -3.54367 - 1.09988 * c3
          XRS_HFRAC = exp( temp )
        end if

        end


       
        real function XRS_HMFRAC( c )

        implicit none

        real c

*       Fraction of Hi- and Med- res events, for an on-axis point source.
*         Tchebycheff polinomial fit to the ln of fraction as function
*         of count rate (* 0.001); linear extrapolation outside the
*         simulation range.
*       Written by KM, 2004 June 15

        double precision c3, temp
        integer k
        double precision coeffs( 0: 4 )
        data coeffs / 0.00103695749, -9.62000451, 17.9504393,
     &                -22.7044576, 10.8511646 /

        if( c .le. 1.0 ) then
          XRS_HMFRAC = 1.0 - c / 1.0 * 0.008528531
        else if( c .le. 822.5821 ) then
          temp = 0.0d+00
          c3 = dble( c ) * 1.0d-03
          do k = 4, 0, -1
            temp = temp * c3 + coeffs( k )
          end do
          XRS_HMFRAC = exp( temp )
        else
          c3 = ( c - 551.37 ) / 271.2121
c         log( 0.070854 ) = -2.64713
c         log( 0.032213 ) - log( 0.070854 ) = -0.78825
          temp = -2.64713 - 0.78825 * c3
          XRS_HMFRAC = exp( temp )
        end if

        end


       
        real function XRS_HMSFRAC( c )

        implicit none

        real c

*       Fraction of Hi-, Med-, and Mid (S) events, for an on-axis point source.
*         Tchebycheff polinomial fit to the ln of fraction as function
*         of count rate (* 0.001); linear extrapolation outside the
*         simulation range.
*       Written by KM, 2004 June 15

        double precision c3, temp
        integer k
        double precision coeffs( 0: 3 )
        data coeffs / 0.00306716614, -4.03737729, 3.42781658,
     &                 -1.64513958 /

        if( c .le. 5.0 ) then
          XRS_HMSFRAC = 1.0 - c / 5.0 * 0.016889989
        else if( c .le. 822.5821 ) then
          temp = 0.0d+00
          c3 = dble( c ) * 1.0d-03
          do k = 3, 0, -1
            temp = temp * c3 + coeffs( k )
          end do
          XRS_HMSFRAC = exp( temp )
        else
          c3 = ( c - 551.37 ) / 271.2121
c         log( 0.25317 ) = -1.37369
c         log( 0.14740 ) - log( 0.25317 ) = -0.54091
          temp = -1.37369 - 0.54091 * c3
          XRS_HMSFRAC = exp( temp )
        end if

        end
       
