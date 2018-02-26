*+HA4_LEVIN
        subroutine HA4_LEVIN( results, n_res )

        implicit none

        integer n_res
        real results( 0: n_res )

*       Description:
*         Outputs HEAO-1 A4 specific information to the screen via PWRITE
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
*         Koji Mukai, 1995 Mar
*-HA4_LEVIN

        character*48 xw_strng
        integer j

        character*16 bounds( 4 )
        data bounds / '     A:   15- 25', '     B:   25- 40',
     &                '     C:   40- 80', '     D:   80-180' /

        if( n_res .ne. 4 ) then
          call PWRITE( 'SEVERE ERROR:: PIMMS is confused in HA4_LEVIN' )
          return
        end if
        call PWRITE( ' ' )
        call PWRITE( 'Results in the 4 Levine et al bands are:' )
        call PWRITE( ' ' )
        call PWRITE( ' Levine  Nominal   Source' )
        call PWRITE( '  Band   E (keV)   (cps)' )
        do j = 1, 4
          xw_strng = bounds( j )
          if( results( j ) .gt. 500000.0 ) then
            write( xw_strng( 18: 26 ), '(1p,e9.2)' ) results( j )
          else if( results( j ) .gt. 1.0 ) then
            write( xw_strng( 18: 26 ), '(f9.1)' ) results( j )
          else
            write( xw_strng( 18: 26 ), '(1p,e9.2)' ) results( j )
          end if
          call PWRITE( xw_strng )
        end do

        end
