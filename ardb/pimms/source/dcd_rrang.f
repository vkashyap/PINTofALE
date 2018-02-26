*+DCD_RRANG
        subroutine DCD_RRANG( string, lo, hi, flag )

        implicit none

        character*( * ) string
        real lo, hi
        integer flag

*       Description:
*         Fortran 77 routine to decode string of the form r1-r2
*
*       Arguments:
*         string      (i) : input string
*         lo, hi      (o) : results
*         flag        (o) : 0 if no errors
*
*       Dependencies:
*         LENTRIM
*
*       Origin:
*         Created by KM
*
*       Author:
*         Koji Mukai, 1993 Feb 05, integer decode version
*         Koji Mukai, 1993 March, real version
*-DCD_RRANG

        integer total, minus

        integer LENTRIM

*       Begin

*       First, find the length of string
        flag = 0
        total = LENTRIM( string )
        if( total .lt. 0 ) then
          flag = -1
          return
        end if

*       Next, find where the "-" is in the string
        minus = INDEX( string, '-' )
        if( minus .le. 0 ) then
          flag = -1
          return
        else if( minus .gt. 10 ) then
          flag = -2
          return
        end if

        call RD_REAL( string( : minus - 1 ), lo, flag )
        if( flag .ge. 0 ) then
          call RD_REAL( string( minus + 1: total ), hi, flag )
        end if

        if( flag .gt. 0 ) then
          flag = 0
        end if

        end
