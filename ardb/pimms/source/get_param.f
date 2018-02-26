	Subroutine INQ_PARAM( num_C, num_F, num_I )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer num_C, num_F, num_I
C
	num_C = np_C
	num_F = np_F
	num_I = np_I
C
	End
C
C
C
	Subroutine PARAM_C( I, Param, p_String )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer I
	Character*( * ) Param, p_String
C
	Integer num_C, num_F, num_I
C
	Call INQ_PARAM( num_C, num_F, num_I )
	If( I .ge. 1 .and. I .le. num_C ) Then
	  Call GET_PARAM_C( I, Param, *100 )
	Else If( I .lt. 1 ) Then
          call PWRITE( 'ERROR:: non-positive number not allowed' )
	Else
          call GET_PSTR( p_string, param )
	End If
        return

 100	continue
*       Alternate return --- error in GET_PARAM_C
C
	End
C
C
C
	Subroutine GET_PARAM_C( I, Param, * )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer I
	Character*( * ) Param
C
	Integer p_Length
C
	If( I .gt. np_C .or. I .le. 0 ) Then
          call PWRITE( 'ERROR in GET_PARAM_C::' //
     &			' Trying to access a non-existent parameter' )
	  Param = ' '
	  Return 1
	Else
	  p_Length = Len( Param )
	  If( p_length .le. pc_F( I ) - pc_B( I ) ) Then
            call PWRITE( 'ERROR in GET_PARAM_C:: ' //
     &				'parameter is too long...truncating' )
	    Param = p_Char( pc_B( I ): pc_B( I ) + p_Length - 1 )
	  Else
	    Param = p_Char( pc_B( I ): pc_F( I ) )
	  End If
	End If
C
	End
C
C
C
	Subroutine PARAM_F( I, Param, Default, p_String )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer I
	Real Param, Default
	Character*( * ) p_String
C
	Integer num_C, num_F, num_I
	Character*80 dummy_Line
C
	Call INQ_PARAM( num_C, num_F, num_I )
	If( I .ge. 1 .and. I .le. num_F ) Then
	  Call GET_PARAM_F( I, Param, *100 )
	Else If( I .lt. 1 ) Then
          call PWRITE( 'ERROR:: non-positive number not allowed' )
	Else
          call GET_PSTR( p_string, dummy_line )
	  If( dummy_Line .eq. ' ' ) Then
	    Param = Default
	  Else
	    Read( dummy_Line, * ) Param
	  End If
	End If
        return

 100	continue
*       alternate return --- error in GET_PARAM_F
C
	End
C
C
C
	Subroutine GET_PARAM_F( I, Param, * )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer I
	Real Param
C
	If( I .gt. np_F .or. I .le. 0 ) Then
          call PWRITE( 'ERROR in GET_PARAM_F:: ' //
     &			'Trying to access a non-existent parameter' )
	  Param = -99.0
	  Return 1
	Else
	  Param = p_Float( I )
	End If
C
	End
C
C
C
	Subroutine PARAM_I( I, Param, Default, p_String )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer I
	Integer Param, Default
	Character*( * ) p_String
C
	Integer num_C, num_F, num_I
	Character*80 dummy_Line
C
	Call INQ_PARAM( num_C, num_F, num_I )
	If( I .ge. 1 .and. I .le. num_I ) Then
	  Call GET_PARAM_I( I, Param, *100 )
	Else If( I .lt. 1 ) Then
          call PWRITE( 'ERROR:: non-positive number not allowed' )
	Else
          call GET_PSTR( p_string, dummy_line )
	  If( dummy_Line .eq. ' ' ) Then
	    Param = Default
	  Else
	    Read( dummy_Line, * ) Param
	  End If
	End If
        return

 100	continue
*       alternate return --- error in GET_PARAM_I
C
	End
C
C
C
	Subroutine GET_PARAM_I( I, Param, * )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer I
	Integer Param
C
	If( I .gt. np_I .or. I .le. 0 ) Then
          call PWRITE( 'ERROR in GET_PARAM_I:: ' //
     &			'Trying to access a non-existent parameter' )
	  Param = -999
	  Return 1
	Else
	  Param = p_Int( I )
	End If
C
	End
C
C
C
	Subroutine PARAM_N( I, Param, Default, p_String )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer I
	Real Param, Default
	Character*( * ) p_String
C
	Integer num_C, num_F, num_I
	Character*80 dummy_Line
C
	Call INQ_PARAM( num_C, num_F, num_I )
	If( I .ge. 1 .and. I .le. num_F + num_I ) Then
	  Call GET_PARAM_N( I, Param, *100 )
	Else If( I .lt. 1 ) Then
          call PWRITE( 'ERROR:: non-positive number not allowed' )
	Else
          call GET_PSTR( p_string, dummy_line )
	  If( dummy_Line .eq. ' ' ) Then
	    Param = Default
	  Else
	    Read( dummy_Line, * ) Param
	  End If
	End If
        return

 100	continue
*	alternate return --- error in GET_PARAM_N
C
	End
C
C
C
	Subroutine GET_PARAM_N( I, Param, * )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer I
	Real Param
C
	If( I .gt. np_N .or. I .le. 0 ) Then
          call PWRITE( 'ERROR in GET_PARAM_N:: ' //
     & 			'Trying to access a non-existent parameter' )
	  Param = -99.0
	  Return 1
	Else If( n_List( I ) .gt. 100 ) Then
	  Param = p_Int( n_List( I ) - 100 )
	Else
	  Param = p_Float( n_List( I ) )
	End If
C
	End



        subroutine PARAM_SUB( nc, np_num, np_beg )

        implicit none

        include 'xcomm.inc'

        integer nc
        integer np_beg
        integer np_num

*       PIMMS will take <model name> <model par> <model par>
*                      [<another model name> <model par> <model par>...]
*       so needs to know how many numerical parameter follows a
*       particular character string parameter (NB range is taken
*       as a character string by this FCI routine).
*
*       nc        (i) : String parameter number
*       np_num    (o) : Number of numerical parameters belonging to this string
*                       -99 for error, 0 for none
*       np_beg    (o) : First numerical parameter (pass to PARAM_N)

        np_beg = po_c( nc ) + 1
        if( nc .le. 0 .or. nc .gt. np_c ) then
*         Invalid number
          np_num = -99
        else if( nc .eq. np_c ) then
          np_num = np_total - np_beg + 1
        else
          np_num = po_c( nc + 1 ) - np_beg
        end if

        end



        subroutine GT2_PARAM_N( i, param )

        implicit none

        include 'xcomm.inc'

        integer i
        real param

        integer ii
        do ii = 1, np_n
          if( po_n( ii ) .eq. i ) then
            if( n_list( ii ) .gt. 100 ) then
              param = p_int( n_list( ii ) - 100 )
            else
              param = p_float( n_list( ii ) )
            end if
            return
          end if
        end do
        call PWRITE( 'ERROR in GT2_PARAM_N:: ' //
     & 			'Trying to access a non-existent parameter' )
        param = -99.0

        end



        subroutine GET_PSTR( query, result )

        implicit none

        include 'xcomm.inc'

        character*( * ) query, result

        integer q_len, r_len
        character*256 w_strng
        integer LENTRIM

        q_len = LENTRIM( query )
        if( from_file ) then
          w_strng = '! ' // query( : q_len )
          call WRITEN( w_strng )
          read( file_unit, '(a)' ) result
          r_len = LENTRIM( result )
          write( *, '(a)' ) result( : r_len )
        else
          call WRITEN( query( : q_len ) )
          read( *, '(a)' ) result
          r_len = LENTRIM( result )
        end if
        w_strng = '! ' // query( : q_len ) // ' ' // result( : r_len )
        call PWRITEL( w_strng )

        end
