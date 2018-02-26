*       Major re-write in 2015 July based on discussion with Michael Smith
*       to remove the obsolete mode, "Refreshed", for MOS and with post-launch
*       estimates of pile-up and dead-time fractions.

        subroutine EPIC_MOS( cps_given, n_res )

        implicit none

        real ln10, psffactor
        parameter( ln10 = 2.302585093 )
*         Using exp( ln10 * x ) instead of 10.0**x for fun and obfuscation
*        parameter( psffactor = 1.47 )
        parameter( psffactor = 1.0 )

        integer n_res
        real cps_given( 0: n_res )

*       A Fortran implimentation of the pile-up fraction formula given
*         by David Lumb
*       Last modified on 1999 Mar 07 to handle high count rate case
*         correctly
*       Modified back to psffactor = 1.0 coresponding to the new
*         effective area curve files for large extraction reginos,
*         as requested by XMM project at ESA.  2010 April 15 for Version 4.1.

        real work, pperc, factor, logr
        integer j, k
*        logical piledup
        character*80 p_strng

        character*9 NINE

        real coeffs( 0: 2, 3 )
        real dtf( 4 )
        character*10 desc( 4 )
        data coeffs / -1.332, 1.020, -0.185,
     &                -1.841, 1.191, -0.185,
     &                -2.452, 1.368, -0.185 /
        data dtf / 0.020, 0.007, 0.15, 0.0 /
        data desc / 'Full',
     &              'Large',
     &              'Small',
C     &              'Refreshed',               OBSOLETE - removed
     &              'Timing' /

        call PWRITE( ' ' )
        p_strng = '% Pile-up and dead-time corrected '
     &                                // 'count rates in 4 energy bands'
        call PWRITE( p_strng )
        call PWRITE( '  using various window options are:' )
        call PWRITE( ' ' )
        p_strng = ' Window    Pileup  Dead                ' //
     &                                      'Corrected Good Count Rates'
        call PWRITE( p_strng )
        p_strng = '   Option   frac.   Time    0.1-0.4   ' //
     &                           '0.4-1.0   1.0-2.5  2.5-10.0     Total'
        call PWRITE( p_strng )
        call PWRITE( ' ' )

        do j = 1, 3
          logr = log( cps_given( 0 ) * psffactor ) / ln10

*         When cps_given( 0 ) is for 15" extraction region, multiply by
*         psffactor=1.47 to get the total rate; now it's for a large
*         extraction region, so psffactor=1.00

          p_strng = desc( j )

*          piledup = .false.
          work = 0.0
          do k = 2, 0, -1
            work = work * logr + coeffs( k, j )
          end do
          pperc = exp( work * ln10 ) * 100.0
*         For humans, use % (that's why we multiyply by 100)
          work = pperc * 0.01
          if( pperc .lt. 0.001 ) then
            p_strng( 11: 17 ) = '<0.001%'
          else
            write( p_strng( 11: 17 ), '(F6.3,''%'')' ) pperc
          end if
          pperc = 100.0 * dtf( j )
          write( p_strng( 20: 24 ), '(F4.1,''%'')' ) pperc
          factor = ( 1.0 - work ) * ( 1.0 - dtf( j ) )
          p_strng( 27: 35 ) = NINE( factor * cps_given( 1 ) )
          p_strng( 37: 45 ) = NINE( factor * cps_given( 2 ) )
          p_strng( 47: 55 ) = NINE( factor * cps_given( 3 ) )
          p_strng( 57: 65 ) = NINE( factor * cps_given( 4 ) )
          p_strng( 69: 77 ) = NINE( factor * cps_given( 0 ) )
*          if( piledup ) p_strng( 79: 79 ) = '?'
          call PWRITE( p_strng )
        end do

*       Timing mode is special
        p_strng = desc( 4 )
        p_strng( 11: 17 ) = '-------'
        p_strng( 21: 24 ) = '0.0%'
        p_strng( 27: 35 ) = NINE( cps_given( 1 ) )
        p_strng( 37: 45 ) = NINE( cps_given( 2 ) )
        p_strng( 47: 55 ) = NINE( cps_given( 3 ) )
        p_strng( 57: 65 ) = NINE( cps_given( 4 ) )
        p_strng( 69: 77 ) = NINE( cps_given( 0 ) )
        call PWRITE( p_strng )
*
        call PWRITE( ' ' )

        end


        subroutine EPIC_PN( cps_given, n_res )

        implicit none

        real ln10, psffactor
        parameter( ln10 = 2.302585093 )
*         Using exp( ln10 * x ) instead of 10.0**x for fun and obfuscation
*        parameter( psffactor = 1.47 )
        parameter( psffactor = 1.00 )

        integer n_res
        real cps_given( 0: n_res )

*       A Fortran implimentation of the pile-up fraction formula given
*         by David Lumb
*       Modified on 1999 Mar 07 to handle high count rate case
*         correctly
*       Modified on 2006 April 14 --- changed dead time fraction for
*         burst mode from 99.8% (pre-launch estimate?) to 97%
*         a la XMM Users' handbook
*       Last modified in 2007 September, corresponding to the files that
*         have the PATTERN=0 only count rate mostly, with the last column
*         haveing total (including other PATTERNS) count rates.
*       Comments modified 2010 April 15 for Version 4.1

        real logc, work, pperc, factor
        integer j, k
        logical piledup
        character*80 p_strng

        character*9 NINE

        real coeffs( 0:2, 4 ), pu_lim( 4 )
        real dtf( 6 )
        character*9 desc( 6 )
        data coeffs / -1.957,  1.162, -0.134,
     &                -1.306,  1.001, -0.134,
     &                -2.179,  1.212, -0.134,
     &                -3.411,  1.459, -0.134 /
        data pu_lim / 75.0, 20.0, 150.0, 1000.0 /
        data dtf / 0.070, 0.070, 0.090, 0.290, 0.015, 0.998 /
        data desc / ' Full',
     &              ' FullExtd',
     &              ' Large',
     &              ' Small',
     &              ' Timing',
     &              ' Burst' /

        logc = log( cps_given( 0 ) * psffactor ) / ln10

*       cps_given( 5 ) is for 15" extraction region, all PATTERNs,
*       so multiply by psffactor=1.47 to get the total rate, even though
*       cps_given( 0 ) is back again to being for a large extraction region
*       Or simply cps_given( 0 ) with a psffactor of 1.0

        call PWRITE( ' ' )
        p_strng = '% Pile-up and dead-time corrected '
     &                                // 'count rates in 4 energy bands'
        call PWRITE( p_strng )
        call PWRITE( '  using various window options are:' )
        call PWRITE( ' ' )
        p_strng = ' Window    Pileup  Dead                ' //
     &                                      'Corrected Good Count Rates'
        call PWRITE( p_strng )
        p_strng = '   Option   frac.   Time    0.1-0.4   ' //
     &                           '0.4-1.0   1.0-2.5  2.5-10.0     Total'
        call PWRITE( p_strng )
        call PWRITE( ' ' )

        do j = 1, 4
          p_strng = desc( j )
          piledup = .false.
          if( cps_given( 0 ) .le. pu_lim( j ) ) then
            work = 0.0
            do k = 2, 0, -1
              work = work * logc + coeffs( k, j )
            end do
            pperc = 100.0 * exp( work * ln10 )
*           For humans, use % (that's why we multiyply by 100)
            work = pperc * 0.01
            if( pperc .gt. 10.0 ) then
              p_strng( 14: 17 ) = '>10%'
              if( pperc .ge. 100.0 ) work = 1.0
              piledup = .true.
            else if( pperc .lt. 0.001 ) then
              p_strng( 11: 17 ) = '<0.001%'
            else
              write( p_strng( 11: 17 ), '(F6.3,''%'')' ) pperc
            end if
          else
            piledup = .true.
            work = 0.1
            p_strng( 14: 17 ) = '>10%'
          end if
          pperc = 100.0 * dtf( j )
          write( p_strng( 20: 24 ), '(F4.1,''%'')' ) pperc
          factor = ( 1.0 - work ) * ( 1.0 - dtf( j ) )
          p_strng( 27: 35 ) = NINE( factor * cps_given( 1 ) )
          p_strng( 37: 45 ) = NINE( factor * cps_given( 2 ) )
          p_strng( 47: 55 ) = NINE( factor * cps_given( 3 ) )
          p_strng( 57: 65 ) = NINE( factor * cps_given( 4 ) )
          p_strng( 69: 77 ) = NINE( factor * cps_given( 0 ) )
          if( piledup ) p_strng( 79: 79 ) = '?'
          call PWRITE( p_strng )
        end do

*       Timing mode is special
        p_strng = desc( 5 )
        p_strng( 11: 17 ) = '-------'
        pperc = 100.0 * dtf( 5 )
        write( p_strng( 20: 24 ), '(F4.1,''%'')' ) pperc
        factor = 1.0 - dtf( 5 )
        p_strng( 27: 35 ) = NINE( factor * cps_given( 1 ) )
        p_strng( 37: 45 ) = NINE( factor * cps_given( 2 ) )
        p_strng( 47: 55 ) = NINE( factor * cps_given( 3 ) )
        p_strng( 57: 65 ) = NINE( factor * cps_given( 4 ) )
        p_strng( 69: 77 ) = NINE( factor * cps_given( 0 ) )
        call PWRITE( p_strng )

*       Burst mode is special
        p_strng = desc( 6 )
        p_strng( 11: 17 ) = '-------'
        pperc = 100.0 * dtf( 6 )
        write( p_strng( 20: 24 ), '(F4.1,''%'')' ) pperc
        factor = 1.0 - dtf( 6 )
        p_strng( 27: 35 ) = NINE( factor * cps_given( 1 ) )
        p_strng( 37: 45 ) = NINE( factor * cps_given( 2 ) )
        p_strng( 47: 55 ) = NINE( factor * cps_given( 3 ) )
        p_strng( 57: 65 ) = NINE( factor * cps_given( 4 ) )
        p_strng( 69: 77 ) = NINE( factor * cps_given( 0 ) )
        call PWRITE( p_strng )
*
        call PWRITE( ' ' )
        call PWRITE(
     &       '% Any pile-up predictions over 10% are highly uncertain' )
        call PWRITE( ' ' )

        end



        subroutine EPIC_MOS_OPEN( cps_given )

        implicit none

        real ln10
        parameter( ln10 = 2.302585093 )
*         Using exp( ln10 * x ) instead of 10.0**x for fun and obfuscation

        real cps_given

*       A Fortran implimentation of the pile-up fraction formula given
*         by David Lumb
*       This assumes that "events per second DETECTED" means
*         "events per second DETECTED, if it weren't for the pile-up"

        real logc, work, pperc
        integer j, k
        character*80 p_strng

        real coeffs( 0:2, 3 )
        character*24 desc( 3 )
        data coeffs / -1.8,  0.92, -0.215,
     &                -2.8,  1.31, -0.215,
     &                -2.08, 1.04, -0.215 /
        data desc / ' MOS Full',
     &              ' MOS Window 2 (=100x100)',
     &              ' MOS Window 3 (=400x400)' /

        logc = log10( cps_given )

        do j = 1, 3
          p_strng = '* For ' // desc( j )
          call PWRITE( p_strng )
          work = 0.0
          do k = 2, 0, -1
            work = work * logc + coeffs( k, j )
          end do
          pperc = 100.0 * exp( work * ln10 )
*         For humans, use % (that's why we multiyply by 100)
          if( pperc .ge. 100.0 ) then
            p_strng = '  >=100% pile-up calculated, check the numbers?'
          else if( pperc .lt. 0.001 ) then
            p_strng = '  Pile-up is predicted to be <0.001%'
          else
            write( p_strng,
     &        '(''  Pile-up is predicted to be '', f6.3, ''%'')' ) pperc
          end if
          call PWRITE( p_strng )
        end do

        end


        subroutine EPIC_PN_OPEN( cps_given, special, restr )

        implicit none

        real cps_given
        integer special
*                    611-614 are XMM PN Open, Thin, Med, and Thick
        logical restr

*       This routine lists the count rates in different window mode
*         and calls PN_PILEUP when appropriate.

        real scaleit, cps_new
        integer k
        character*80 p_strng

        real open_scales( 2: 5 ), other_scales( 2: 5 )
        character*6 winmode( 2: 5 )
        data open_scales / 0.947368, 0.736842, 1.03684, 0.00178947 /
        data other_scales / 0.967742, 0.752688, 1.05914, 0.00182796 /
        data winmode / 'LARGE', 'SMALL', 'TIMING', 'BURST' /

        call PWRITE( '  in the FULL window mode.' )
        if( .not. restr ) then
          call PN_OPEN_PILEUP( cps_given, 1 )
        end if

        do k = 2, 5
          if( special .eq. 611 ) then
            scaleit = open_scales( k )
          else
            scaleit = other_scales( k )
          end if
          cps_new = cps_given * scaleit
          call PWRITE( ' ' )
          p_strng = '% In '  // winmode( k ) // ' window mode'
          call PWRITE( p_strng )
          p_strng = '  The rate is FULL mode rate * 0.00183 = '
          write( p_strng( 32: 38 ), '(f7.5)' ) scaleit
          if( cps_new .gt. 1000.0 ) then
            write( p_strng( 42: 51 ), '(1p,e10.3)' ) cps_new
          else if( cps_new .ge. 0.001 ) then
            write( p_strng( 42: 51 ), '(f10.5)' ) cps_new
          else
            write( p_strng( 42: 51 ), '(1p,e10.3)' ) cps_new
          end if
          call PWRITE( p_strng )
          if( .not. restr ) then
            call PN_OPEN_PILEUP( cps_new, k )
          end if
        end do

        end


        subroutine PN_OPEN_PILEUP( cps_given, window )

        implicit none

        real ln10
        parameter( ln10 = 2.302585093 )
*         Using exp( ln10 * x ) instead of 10.0**x for fun and obfuscation

        real cps_given
        integer window

*       A Fortran implimentation of the pile-up fraction formula given
*         by David Lumb

        real logc, work, pperc
        integer k
        character*80 p_strng

        real coeffs( 0:2, 1:3 )
*                         1 for Full, 2 for Large, and 3 for Small
        data coeffs / -3.348, 1.55, -0.174,
     &                -3.963, 1.68, -0.174,
     &                -5.48,  1.97, -0.174 /

        if( window .ge. 1 .and. window .le. 3 ) then
*         Simply return if window is 4 (Timing) or 5 (Burst)

          logc = log10( cps_given )

          work = 0.0
          do k = 2, 0, -1
            work = work * logc + coeffs( k, window )
          end do
          pperc = 100.0 * exp( work * ln10 )
*         For humans, use % (that's why we multiyply by 100)
          if( pperc .ge. 100.0 ) then
            p_strng = '* >=100% pile-up calculated, check the numbers?'
          else if( pperc .lt. 0.001 ) then
            p_strng = '* Pile-up is predicted to be <0.001%'
          else
            write( p_strng,
     &        '(''* Pile-up is predicted to be '', f6.3, ''%'')' ) pperc
          end if
          call PWRITE( p_strng )

        end if

        end



        character*9 Function NINE( num )

        implicit none

        real num

*       Writes a real number into a string using Fortran format
*       	1PE9.4 if it is less than 0.01
*               F9.4 if it is between 0.01 and 999.49999
*               1PE9.4 again if it is greater than or equal to 999.5
*       Might be ugly if the num is negative

        character*9 buffer

        if( num .lt. 0.01 ) then
          write( buffer, '(1P,E9.2)' ) num
        else if( num .lt. 999.5 ) then
          write( buffer, '(F9.4)' ) num
        else
          write( buffer, '(1P,E9.2)' ) num
        end if

        NINE = buffer

        end
