	Subroutine FORCE_COMM( force_Line )
C
	Implicit None

        include 'xcomm.inc'
C
	Character*( * ) force_Line
C
	Integer line_Length
C
	Integer K
C
	line_Length = Len( force_Line )
	If( line_Length .gt. l_Buffer ) Then
	  Stop 'FCI FATAL ERROR:: Command line too long'
	End If
	If( .not. Empty ) Then
	  n_Buffer = n_Buffer + 1
	  If( n_Buffer .gt. m_Buffer ) Then
	    Stop 'FCI FATAL ERROR:: Command buffer overload'
	  End If
	  If( n_Buffer .eq. 1 ) Then
	    buff_Line( 1 ) = current_Line
	    current_Line = force_Line
	  Else
	    K = n_Buffer - 1
	    buff_Line( n_Buffer ) = buff_Line( K )
	    buff_Line( K ) = force_Line
	  End If
	Else
	  n_Buffer = n_Buffer + 1
	  Do K = n_Buffer, 2, -1
	    buff_Line( K ) = buff_Line( K - 1 )
	  End Do
	  buff_Line( 1 ) = force_Line
	End If
C
	End
