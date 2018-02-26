*+SIS_LIMIT
        subroutine SIS_LIMIT( result )

        implicit none

        integer m_type
        parameter( m_type = 2 )

        real result

*       Description:
*         Outputs SIS specific information to the screen via PWRITE
*
*       Arguments:
*         result (i)  : Predicted count rate
*
*       Dependencies:
*         SIS_OBSTM, PWRITE
*
*       Origin:
*         Created by KM
*
*       Author
*         Koji Mukai, 1993 Mar, first official version
*         Koji Mukai, 1993 Apr 6 --- patch to avoid crashes at high count rates.
*         Koji Mukai, 1995 Aug: A major update to make it AO-4 compatible
*-SIS_LIMIT

        character*80 pw_strng
        real fm, fh, bm, bh

        real time, r_opt, s_opt, b_opt
        integer i

        real f_rate( m_type )
        character*9 f_note( m_type )
        data f_rate / 5.334, 6.697 /
        data f_note / '(nominal)', '(warm)' /

        call PWRITE( ' ' )
        if( result .gt. 1.0e+04 ) then
          call PWRITE( '* The source is apparently too bright for SIS' )
          call PWRITE( '  If this is unexpected, '
     &                            // 'check the unit of input rate' )
        else
          call SIS_OBSTM( 2, result, 2.12e-07, time, 5.0,
     &                                             r_opt, s_opt, b_opt )
*         SIS_OBSTM calculates r_opt in unbinned pixels, times 1.59/60
*         to put it in arcmin
          r_opt = r_opt * 0.0265
          if( time .gt. 1.0e+07 ) then
            pw_strng = '* It appears to take more than 10 million s'
     &                                     // ' for a 5-sigma detection'
            call PWRITE( pw_strng( : 68 ) )
          else
            write( pw_strng, 105 ) time
105         format( '* An exposure of ', f10.2, 's is required for a',
     &                              ' 5-sigma detection using optimal' )
            call PWRITE( pw_strng( : 79 ) )
            write( pw_strng, 106 ) r_opt, s_opt, b_opt
106         format( '  extraction radius of ', f4.2, ' arcmin (', 1p,
     &                            e9.2, ' src and ', e9.2, ' bgd cps)' )
            call PWRITE( pw_strng( : 72 ) )
          end if
          call PWRITE( ' ' )
          call PWRITE( '* Extrapolated to 1998 Dec for S1C3 in ' //
     &                             '1-CCD mode at 2 SIS temperatures:' )
          call PWRITE( ' ' )
          do i = 1, m_type
            fm = ( result + f_rate( i ) ) / 8.0 * 100
            fh = ( result + f_rate( i ) ) / 64.0 * 100
            bm = ( result + f_rate( i ) ) / 32.0 * 100
            bh = ( result + f_rate( i ) ) / 256.0 * 100
            if( fm .lt. 99.95 ) then
              write( pw_strng, 110 ) fm, fh, bm, bh, f_note( i )
110           format( '  Telemetry is ', f4.1, '% full, ', f4.1,
     &           '% full, ', f4.1, '% full, and ', f4.1, '% full ', a9 )
            else if( bm .lt. 99.95 ) then
              write( pw_strng, 111 ) fh, bm, bh, f_note( i )
111           format( '  Telemetry is  saturated, ', f4.1,
     &           '% full, ', f4.1, '% full, and ', f4.1, '% full ', a9 )
            else if( fh .lt. 99.95 ) then
              write( pw_strng, 112 ) fh, bh, f_note( i )
112           format( '  Telemetry is  saturated, ', f4.1,
     &                 '% full,  saturated, and ', f4.1, '% full ', a9 )
            else if( bh .lt. 99.95 ) then
              write( pw_strng, 113 ) bh, f_note( i )
113           format( '  Telemetry is  saturated,  saturated,',
     &                        '  saturated, and ', f4.1, '% full ', a9 )
            else
              pw_strng = '  Telemetry is  saturated,  saturated,' //
     &                     '  saturated, and  saturated ' // f_note( i )
            end if
            call PWRITE( pw_strng( : 75 ) )
          end do
          pw_Strng = '         in    Faint  (M), Faint  (H), ' //
     &                                'Bright (M), and Bright (H) modes'
          call PWRITE( pw_strng )
          pw_Strng = '    (S0C1 has lower flickering pixel rate; ' //
     &                                     'decide modes based on S1C3)'
          call PWRITE( pw_Strng )
        end if

        end
