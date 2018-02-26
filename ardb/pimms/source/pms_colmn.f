*+PMS_COLMN

        subroutine PMS_COLMN( lun, line, column )

        implicit none

        integer lun
        character*( * ) line
        integer column

*       Description:
*         Count number of columns in an ASCII file
*
*       Arguments:
*         lun         (i) : Input logical unit number for the file
*         line        (o) : A line in the file
*         column      (o) : The number of columns in that line
*
*       Dependencies:
*         LENTRIM
*
*       Origin:
*         Lost in time
*
*       Author:
*         Koji Mukai, 1994 Nov, original version
*-PMS_COLMN

	integer len, l, m
	logical on
	integer LENTRIM

        column = 0
        line = ' '
        read( lun, '(a)', end = 900 ) line
        len = LENTRIM( line )
        on = .false.
        do l = 1, len
          m = ichar( line( l: l ) )
          if( ( m .le. 32 .or. m .eq. 44 ) .and. on ) then
            on = .false.
          else if( ( m .gt. 32 .and. m .ne. 44 ) .and. .not. on ) then
            on = .true.
            column = column + 1
          end if
        end do
        return

900     continue
        column = -1

        end
