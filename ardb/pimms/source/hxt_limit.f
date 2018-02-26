*+HXT_LIMIT
        subroutine HXT_LIMIT( results, n_res, special )

        implicit none

        integer n_res
        real results( 0: n_res )
        integer special

*       Description:
*         Outputs XTE HEXTE specific information to the screen via PWRITE
*
*       Arguments:
*         results (i)  : Predicted count rates
*         n_res   (i)  : Number of standard bands
*
*       Dependencies:
*         PWRITE
*
*       Origin:
*         Created by KM
*
*       Author
*         Koji Mukai, 1994 Dec
*         Koji Mukai, 1996 May --- Modified to handle different LLDs
*         Koji Mukai, 1997 Jun --- Using the more rigorous bgd systematics
*-HXT_LIMIT

        character*48 xw_strng
        character*70 bgd_strng
        real temp, thisbgd, thissys, fivesigma
        integer j, i

        real bgd_rate( 0: 4, 0: 5 ), bgd_syst( 0: 4, 0: 5 )
        character*20 bounds( 4, 0: 5 )
        data bgd_rate / 73.08, 11.86, 17.93, 21.94, 21.35,
     &                  74.06, 12.84, 17.93, 21.94, 21.35,
     &                  71.57, 10.34, 17.93, 21.94, 21.35,
     &                  69.33,  8.10, 17.93, 21.94, 21.35,
     &                  67.54,  6.31, 17.93, 21.94, 21.35,
     &                  61.22, 17.93, 21.94, 21.35,  0.00 /
        data bgd_syst / 2.39e-2, 1.59e-2, 5.88e-3, 1.93e-3, 2.31e-4,
     &                  2.64e-2, 1.84e-2, 5.88e-3, 1.93e-3, 2.31e-4,
     &                  1.99e-2, 1.19e-2, 5.88e-3, 1.93e-3, 2.31e-4,
     &                  1.46e-2, 6.56e-3, 5.88e-3, 1.93e-3, 2.31e-4,
     &                  1.08e-2, 2.76e-3, 5.88e-3, 1.93e-3, 2.31e-4,
     &                  8.04e-3, 5.88e-3, 1.93e-3, 2.31e-4, 0.00e-0 /
        data bounds / '    5- 29   12-  30', '   30- 61   30-  62',
     &                '   62-125   62- 126', '  126-250  126- 250',
     &                '    5- 29   10-  30', '   30- 61   30-  62',
     &                '   62-125   62- 126', '  126-250  126- 250',
     &                '    5- 29   15-  30', '   30- 61   30-  62',
     &                '   62-125   62- 126', '  126-250  126- 250',
     &                '    5- 29   20-  30', '   30- 61   30-  62',
     &                '   62-125   62- 126', '  126-250  126- 250',
     &                '    5- 29   25-  30', '   30- 61   30-  62',
     &                '   62-125   62- 126', '  126-250  126- 250',
     &                '   30- 61   30-  62', '   62-125   62- 126',
     &                '  126-250  126- 250', '                   ' /
C        data bgd_rate / 34.4, 4.52, 5.66, 7.77, 16.4 /
C        data bounds / '   15- 29   15-  30', '   30- 61   30-  61',
C     &                '   62-125   61- 125', '  126-250  125- 250' /

        if( special .eq. 202 ) then
          i = 0
        else
          i = special - 210
        end if
        if( i .eq. 5 ) then
          if( n_res .ne. 3 ) then
            call PWRITE(
     &                 'SEVERE ERROR:: PIMMS is confused in HXT_LIMIT' )
            return
          end if
        else
          if( n_res .ne. 4 ) then
            call PWRITE(
     &                 'SEVERE ERROR:: PIMMS is confused in HXT_LIMIT' )
            return
          end if
        end if

        thisbgd = bgd_rate( 0, i )
        thissys = 25.0 * bgd_syst( 0, i ) * bgd_syst( 0, i )
        bgd_strng = '   (Source-only count rate in 1 cluster; ' //
     &                                   'BGD rate is 34.4 per cluster)'
        write( bgd_strng( 53: 57 ), '(f5.1)' ) thisbgd
        call PWRITE( bgd_strng )
        if( results( 0 ) .gt. 1e15 ) then
          xw_strng = 'Predicted count rate very high: mistake maybe?'
        else
          temp = results( 0 ) * results( 0 )
          if( temp .gt. thissys ) then
            fivesigma = 25.0 * ( results( 0 ) + 2.33 * thisbgd )
     &                                              / ( temp - thissys )
            xw_strng = '5-sigma detection will be achieved in         s'
            if( fivesigma .gt. 500000.0 ) then
              write( xw_strng( 38: 46 ), '(1p,e9.2)' ) fivesigma
            else if( fivesigma .gt. 1.0 ) then
              write( xw_strng( 38: 46 ), '(f9.1)' ) fivesigma
            else
              write( xw_strng( 38: 46 ), '(1p,e9.2)' ) fivesigma
            end if
          else
            xw_strng =
     &                'Count rate too low: 5-sigma detection impossible'
          end if
        end if
        call PWRITE( xw_strng )
        call PWRITE( ' ' )
        call PWRITE( 'Results in the 4 canonical XTE HEXTE bands are:' )
        call PWRITE( '             (per HEXTE cluster)' )
        call PWRITE( ' ' )
        call PWRITE( '  Channels Nominal   Source   BGD   5-sigma' )
        call PWRITE(
     &         '           E (keV)   (cps)    (cps) detection (s)' )
        if( i .eq. 5 ) then
          call PWRITE( '    0- 29   30-  30    0.0     0.0  ********' )
        end if
        do j = 1, n_res
          xw_strng = bounds( j, i )
          if( results( j ) .gt. 500.0 ) then
            write( xw_strng( 21: 29 ), '(1p,e9.2)' ) results( j )
          else if( results( j ) .gt. 1.0 ) then
            write( xw_strng( 23: 27 ), '(f5.1)' ) results( j )
          else
            write( xw_strng( 21: 29 ), '(1p,e9.2)' ) results( j )
          end if
          thisbgd = bgd_rate( j, i )
          thissys = 25.0 * bgd_syst( j, i ) * bgd_syst( j, i )
          write( xw_strng( 30: 35 ), '(f6.2)' ) thisbgd
          if( results( j ) .gt. 1e15 ) then
            xw_strng( 36: 44 ) = ' ********'
          else
            temp = results( j ) * results( j )
            if( temp .gt. thissys ) then
              fivesigma = 25.0 * ( results( j ) + 2.33 * thisbgd )
     &                                              / ( temp - thissys )
              if( fivesigma .gt. 500000.0 ) then
                write( xw_strng( 36: 44 ), '(1p,e9.2)' ) fivesigma
              else if( fivesigma .gt. 1.0 ) then
                write( xw_strng( 36: 44 ), '(f9.1)' ) fivesigma
              else
                write( xw_strng( 36: 44 ), '(1p,e9.2)' ) fivesigma
              end if
            else
              xw_strng( 36: 44 ) = ' ********'
            end if
          end if
          call PWRITE( xw_strng )
        end do
        call PWRITE(
     & ' (The default 16-s rocking cycle is assumed for detection time)'
     &                )

        end
