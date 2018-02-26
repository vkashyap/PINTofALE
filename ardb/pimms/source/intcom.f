        block data INTCOM_DATA
*
*       Initialize a common block variable
*                            JCM noticed its necessity 1997 Sep 27
*
        include 'xcomm.inc'

        data empty, from_file / .true., .false. /

        end


        subroutine INTCOM( list_in, n_comm_in, prompt_in )

        implicit none

        include 'xcomm.inc'

        integer n_comm_in
        character*( * ) list_in( n_comm_in )
        character*( * ) prompt_in

        character*( l_comm ) temp1, temp2
        integer i, j, k, l, labor1, labor2

        integer FCI_LEN

        character*( l_comm ) list_reserved( n_res )
C        data empty / .true. /          Moved to block data
	Data list_Reserved / 	'STOP', 'HELP' /
*                          also <EOF>    also ?
C
	Length = FCI_LEN( Prompt_in )
	If( n_Comm_in .gt. m_Comm ) Then
	  Stop 'FCI FATAL ERROR: Too many commands, cannot initialize'
	Else If( Length .gt. l_Comm ) Then
	  Stop 'FCI FATAL ERROR: Prompt string too long to initialize'
	End if
	Prompt = Prompt_in
	n_Comm = n_Comm_in
C
C		Sort the commands in alphabetical order
	List( 1 ) = List_in( 1 )
	Call FCI_UPPER_CASE( List( 1 ), l_Comm )
	Do L = 1, n_Res
	  If( List( 1 ) .eq. list_Reserved( L ) ) Then
	    fci_Reserved( L ) = 1
	  Else
	    fci_Reserved( L ) = 0
	  End If
	End Do
	Order( 1 ) = 1
	Do I = 2, n_Comm_in
	  Temp1 = List_in( I )
	  Labor1 = I
	  Call FCI_UPPER_CASE( Temp1, l_Comm )
	  Do L = 1, n_Res
	    If( Temp1 .eq. list_Reserved( L ) ) Then
	      fci_Reserved( L ) = I
	    End If
	  End Do
	  J = 1
	  Do While( J .lt. I .and. Temp1 .gt. List( J ) )
	    J = J + 1
	  End Do
	  If( J .lt. I ) Then
	    Do K = J, I - 1
	      Temp2 = List( K )
	      Labor2 = Order( K )
	      List( K ) = Temp1
	      Order( K ) = Labor1
	      Temp1 = Temp2
	      Labor1 = Labor2
	    End Do
	  End If
	  K = I
	  List( K ) = Temp1
	  Order( K ) = Labor1
	End Do
C
C		Now see if the reserved commands are in the list; if not, add
	I = n_Comm_in
	Do L = 1, n_Res
	  If( fci_Reserved( L ) .eq. 0 ) Then
	    I = I + 1
	    fci_Reserved( L ) = I
	    If( I .gt. m_Comm ) Then
	      Stop
     &           'FCI FATAL ERROR: Too many commands, cannot initialize'
	    End If
	    Temp1 = list_Reserved( L )
	    Labor1 = I
	    J = 1
	    Do While( J .lt. I .and. Temp1 .gt. List( J ) )
	      J = J + 1
	    End Do
	    If( J .lt. I ) Then
	      Do K = J, I - 1
		Temp2 = List( K )
		Labor2 = Order( K )
		List( K ) = Temp1
		Order( K ) = Labor1
		Temp1 = Temp2
		Labor1 = Labor2
	      End Do
	    End If
	    K = I
	    List( K ) = Temp1
	    Order( K ) = Labor1
	  End If
	End Do
	n_Comm = I
C
C       Now check if the software was invoked with a command line argument
C       and, if so, copy it as the current line.  (KM, 2000 Mar 29)
C
        call ARKGCL( current_line )
        if( current_line .ne. ' ' ) empty = .false.
C
	end
