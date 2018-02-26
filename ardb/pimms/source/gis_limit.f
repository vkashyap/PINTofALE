*+GIS_LIMIT
        subroutine GIS_LIMIT( result )

        implicit none

        real result

*       Description:
*         Outputs GIS specific information to the screen via PWRITE
*
*       Arguments:
*         result (i)  : Predicted count rate
*
*       Dependencies:
*         PWRITE
*
*       Origin:
*         Created by KM
*
*       Author
*         Koji Mukai, 1993 Mar, first official version
*         Koji Mukai, 1995 Aug, radically updated recommendation
*-GIS_LIMIT

        character*78 pw_strng
        real sandb, rh, rm, rl
        integer irh, irm, irl

        real bgd
        data bgd / 1.0 /

*       Detection limit/time for detction info not available

*       Now the telemtry limit (128, 16 and 4 cps/GIS limits are hardwired)

        sandb = result + bgd
        rh = ( sandb / 128.0 ) * 100.0
        irh = nint( rh )
        rm = ( sandb / 16.0 ) * 100.0
        irm = nint( rm )
        rl = ( sandb / 4.0 ) * 100.0
        irl = nint( rl )
        pw_strng = '  Telemetry will be saturated, saturated, ' //
     &                            'and saturated in H, M and L bit rate'
        if( sandb .lt. 0.4 ) then
          write( pw_strng( 21: 29 ), '(f3.1,''% full'')' ) rh
          write( pw_strng( 32: 40 ), '(f3.1,''% full'')' ) rm
          write( pw_strng( 47: 55 ), '(f3.1,''% full'')' ) rl
        else if( sandb .lt. 1.6 ) then
          write( pw_strng( 21: 29 ), '(f3.1,''% full'')' ) rh
          write( pw_strng( 32: 40 ), '(f3.1,''% full'')' ) rm
          write( pw_strng( 47: 55 ), '(i3,''% full'')' ) irl
        else if( sandb .lt. 4.0 ) then
          write( pw_strng( 21: 29 ), '(f3.1,''% full'')' ) rh
          write( pw_strng( 32: 40 ), '(i3,''% full'')' ) irm
          write( pw_strng( 47: 55 ), '(i3,''% full'')' ) irl
        else if( sandb .lt. 12.8 ) then
          write( pw_strng( 21: 29 ), '(f3.1,''% full'')' ) rh
          write( pw_strng( 32: 40 ), '(i3,''% full'')' ) irm
        else if( sandb .lt. 16.0 ) then
          write( pw_strng( 21: 29 ), '(i3,''% full'')' ) irh
          write( pw_strng( 32: 40 ), '(i3,''% full'')' ) irm
        else if( sandb .lt. 128.0 ) then
          write( pw_strng( 21: 29 ), '(i3,''% full'')' ) irh
        end if
        calL PWRITE( pw_strng )
        call PWRITE( '  using PH mode (MPC mode is not recommended) ' //
     &                             'allowing for ~1 cps for X-ray and' )
        call PWRITE( '  particle background and the on-board ' //
     &                                           'calibration source.' )

        end
