*+SGD_LIMIT
        subroutine SGD_LIMIT( results, n_res )

        implicit none

        integer n_res_fixed
        parameter( n_res_fixed = 2 )

        integer n_res
        real results( 0: n_res )

*       Description:
*         Outputs ASTRO-H SGD specific information to the screen via PWRITE
*
*       Arguments:
*         results (i)  : Predicted count rates
*         n_res   (i)  : Number of standard bands
*
*       Dependencies:
*         PCA_LIMIT_DO
*         PWRITE
*
*       Origin:
*         Created by KM
*
*       Author
*         Koji Mukai, 2012 May
*-SGD_LIMIT

        real bgd_rate( n_res_fixed ), frac, work, temp, dttime
        integer k
        character*11 bounds( n_res_fixed )
        character*78 xw_strng

        data frac / 0.03 /
        data bgd_rate / 0.057151978, 0.027951659 /
        data bounds / ' 40-150 keV', '150-600 keV' /

        if( n_res .ne. 2 ) then
          call PWRITE( 'SEVERE ERROR:: PIMMS is confused in SGD_LIMIT' )
          return
        end if
        if( results( 0 ) .gt. 1.0e15 ) then
          call PWRITE( '% Predicted count rate very high - error?' )
        else
          call PWRITE( '% Source & bgd count rates in 2 bands '
     &                                // 'and exposure time necessary' )
          xw_strng = '  for 3 and 5 sigma detections (assuming  % '
     &                             // 'bgd systematic uncertainty) are:'
          k = nint( frac * 100.0 )
          write( xw_strng( 42: 42 ), '(i1)' ), k
          call PWRITE( xw_strng )
          call PWRITE( '    Energy     SRC       BGD     '
     &                                          // '3-sigma   5-sigma' )
          do k = 1, 2
            xw_strng = ' '
            xw_strng( 2: 12 ) = bounds( k )
            if( results( k ) .gt. 99.0 ) then
              write( xw_strng( 14: 22 ), '(1p,e9.2)' ) results( k )
            else if( results( k ) .gt. 0.001 ) then
              write( xw_strng( 14: 22 ), '(f9.6)' ) results( k )
            else
              write( xw_strng( 14: 22 ), '(1p,e9.2)' ) results( k )
            end if
            write( xw_strng( 24: 31 ), '(f8.6)' ) bgd_rate( k )
            work = bgd_rate( k ) * frac
            temp = results( k ) * results( k ) / 9.0 - work * work
            if( temp .le. 0.0 ) then
              xw_strng( 33: 41 ) = 'Non-Detct'
            else
              dttime = ( results( k ) + bgd_rate( k ) ) / temp
              if( dttime .ge. 1.0e+6 ) then
                write( xw_strng( 33: 41 ), '(1p,e9.2)' ) dttime
              else
                write( xw_strng( 33: 41 ), '(f9.2)' ) dttime
              end if
            end if
            temp = results( k ) * results( k ) / 25.0 - work * work
            if( temp .le. 0.0 ) then
              xw_strng( 43: 51 ) = 'Non-Detct'
            else
              dttime = ( results( k ) + bgd_rate( k ) ) / temp
              if( dttime .ge. 1.0e+6 ) then
                write( xw_strng( 43: 51 ), '(1p,e9.2)' ) dttime
              else
                write( xw_strng( 43: 51 ), '(f9.2)' ) dttime
              end if
            end if
            call PWRITE( xw_strng )
          end do
        end if

        end
