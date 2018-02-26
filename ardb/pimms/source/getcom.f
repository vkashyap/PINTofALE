	Subroutine GETCOM( choice_Comm )
C
C	New version of Koji's Fortran Command Interpreter
C	This subroutine actually returns the command number
C
c
c       Note on 2000 March 29 --- the command stack (buff_line)
c            never seems to get filled.
c            Added @file capability
c       Note on 2010 Feb 4 --- added echo of the command line to
c            the log file when the input is from the terminal
c            (that is, if the log file is open)
c
	Implicit None

        include 'xcomm.inc'
C
	Integer choice_Comm
C
	Character*( l_Comm ) Short
	Character one_Char
	Integer I, J, K, flag
	Integer beg_P, fin_P, l_Typed
	Logical param_On, quote_On, look_for_Parameters, separator_Found
        character*( l_buffer ) pw_strng
C
1000	Continue
c                			! Restart after a SPAWN
C
	choice_Comm = 0
c               			! Clear
	Call CLR_PARAM( )
C					  Find a command line
	If( Empty ) Then
c                       		! Nothing in the current command line
	  If( n_Buffer .gt. 0 ) Then
c                               	! But something in the Command Stack
	    current_Line = buff_Line( n_Buffer )
	    n_Buffer = n_Buffer - 1
	  Else
c       				! Nothing there, either..read input
            if( from_file ) then
c             If reading from a file...
              read( file_unit, '(A)', end = 390 ) current_line
c             Echo to STDOUT
              pw_strng = '! ' // prompt( 1: length ) // ' '
     &                                                   // current_line
              call PWRITE( pw_strng )
            else
c             read from terminal
              call WRITEN( prompt( 1: length ) )
              Read( *, '(A256)', End = 400 ) current_Line
              pw_strng = '! ' // prompt( 1: length ) // ' '
     &                                                   // current_line
              call PWRITEL( pw_strng )
            end if
	  End If
	End If
C			Find the first non-empty character of the input command
	I = 1
	Call FCI_NOSPC( current_Line, I )
	If ( I .gt. l_Buffer ) then
          choice_comm=-1
          Return
c                ! No commands - return choice=-1
        end if
	If( current_Line( I: I ) .eq. ';' ) Then
          call PWRITE( 'FCI ERROR:: Illegal command separator' )
	  Empty = .True.
	  n_Buffer = 0
          if( from_file ) then
            close( file_unit )
            file_unit = -99
            from_file = .false.
          end if
	  Return
	End If
C
	If( current_Line( I: I ) .eq. '$' ) Then
c                                                       ! Spawn rest of line
          call SPAWN( current_line( i + 1: ) )
	  Empty = .True.
	  Goto 1000
c                					! Go back
	Else If( current_Line( I: I ) .eq. '@' ) Then
c                                                	! @ - input from a file
	  empty = .True.
          if( from_file ) then
c           Already reading from a file
            call PWRITE( 'FCI ERROR:: already reading from a file' )
            goto 1000
          end if
          call ARKOPN( file_unit, ' ', current_line( i + 1: ),
     &                 'xco', 'OLD', 'READONLY', 'FORMATTED',
     &                 'SEQUENTIAL', 1, flag )
          if( flag .lt. 0 ) then
c           couldn't open the command file
            call PWRITE( 'FCI ERROR:: Failed to open command file' )
          else
            from_file = .true.
          end if
          goto 1000
c
	Else If( current_Line( I: I ) .eq. '?' ) Then
c                                                	! ?=HELP
	  Empty = .True.
c                              			! Assume rest of the line and
	  n_Buffer = 0
c                        				! the buffer is garbage
	  choice_Comm = fci_Reserved( 2 )
          if( from_file ) then
            close( file_unit )
            file_unit = -99
            from_file = .false.
          end if
	  Return
	End If
	Call FCI_UPPER_CASE( current_Line, l_Buffer )
C		Find the end of the command
	K = 2
	Short( 1: 1 ) = current_Line( I: I )
	I = I + 1
	separator_Found = .False.
	Do While( I .le. l_Buffer .and. .not. separator_Found )
	  one_Char = current_Line( I: I )
	  If( one_Char .eq. ';' ) Then
	    separator_Found = .True.
	    look_for_Parameters = .False.
	  Else If( one_Char .eq. ' ' ) Then
	    separator_Found = .True.
	    look_for_Parameters = .True.
	  Else
	    Short( K: K ) = one_Char
	    K = K + 1
	    If( K .gt. l_Comm ) Then
	      I = I + 1
	      If( I .le. l_Buffer ) Then
		one_Char = current_Line( I: I )
		If( one_Char .eq. ';' ) Then
		  separator_Found = .True.
		  look_for_Parameters = .False.
		Else If( one_Char .eq. ' ' ) Then
		  separator_Found = .True.
		  look_for_Parameters = .True.
		Else
                  write( pw_strng, 110 ) l_comm
110               format( 'FCI WARNING:: Command too long...',
     &                      'truncated to ', i2, ' letters' )
                  call PWRITE( pw_strng )
		  Do While( I .lt. l_Buffer
     &                            .and. .not. separator_Found )
		    I = I + 1
		    one_Char = current_Line( I: I )
		    If( one_Char .eq. ';' ) Then
		      separator_Found = .True.
		      look_for_Parameters = .False.
		    Else If( one_Char .eq. ' ' ) Then
		      separator_Found = .True.
		      look_for_Parameters = .True.
		    End If
		  End Do
		  If( .not. separator_Found ) Then
		    separator_Found = .True.
		    look_for_Parameters = .False.
		  End If
		End If
	      Else
		separator_Found = .True.
		look_for_Parameters = .False.
	      End If
	    End If
	  End If
	  I = I + 1
	End Do
	l_Typed = K - 1
                                                          
C		Now look for parameters, if necessary
	If( look_for_Parameters ) Then
	  param_On = .False.
	  quote_On = .False.
	  J = I
	  separator_Found = .False.
	  Do While( J .le. l_Buffer .and. .not. separator_Found )
	    one_Char = current_Line( J: J )
	    If( quote_On ) Then
	      If( one_Char .eq. '''' ) Then
		J = J + 1
		If( J .le. l_Buffer ) Then
		  one_Char = current_Line( J: J )
		  If( one_Char .ne. '''' ) Then
		    quote_On = .False.
		    fin_P = J - 1
		    Call SET_PARAM( current_Line, beg_P, fin_P )
		    J = J - 1
		    param_On = .False.
		  End If
	        Else
		  fin_P = l_Buffer
		  Call SET_PARAM( current_Line, beg_P, fin_P )
		  Empty = .True.
		End If
	      End If
	    Else
	      If( one_Char .eq. ';' ) Then
		If( param_On ) Then
		  fin_P = J - 1
		  Call SET_PARAM( current_Line, beg_P, fin_P )
		End If
		Do K = J + 1, l_Buffer
		  current_Line( K - J: K - J ) = current_Line( K: K )
		End Do
		Do K = l_Buffer - J + 1, l_Buffer
		  current_Line( K: K ) = ' '
		End Do
		K = 1
		Call FCI_NOSPC( current_Line, K )
		If( K .gt. l_Buffer ) Then
		  Empty = .True.
		Else
		  Empty = .False.
		End If
		separator_Found = .True.
	      Else If( one_Char .eq. ' ' ) Then
		If( param_On ) Then
		  fin_P = J - 1
		  Call SET_PARAM( current_Line, beg_P, fin_P )
		  param_On = .False.
		End If
	      Else
		If( .not. param_On ) Then
		  If( one_Char .eq. '''' ) Then
		    quote_On = .True.
		  End If
		  beg_P = J
	          param_On = .True.
		End If
	      End If
	    End If
	    J = J + 1
	  End Do
	  If( J .gt. l_Buffer ) Then
	    Empty = .True.
	  End If
	Else
	  If( I .gt. l_Buffer ) Then
	    Empty = .True.
	  Else
	    J = I - 1
	    Do K = J + 1, l_Buffer
	      current_Line( K - J: K - J ) = current_Line( K: K )
	    End Do
	    Do K = l_Buffer - J + 1, l_Buffer
	      current_Line( K: K ) = ' '
	    End Do
	    K = 1
	    Call FCI_NOSPC( current_Line, K )
	    If( K .gt. l_Buffer ) Then
	      Empty = .True.
	    Else
	      Empty = .False.
	    End If
	  End If
	End If
       
C		Now find out which command it is
	Do I = 1, n_Comm
	  If( Short( 1: l_Typed ) .eq. List( I )( 1: l_Typed ) ) Then
	    If( I .eq. n_Comm ) Then
	      choice_Comm = Order( I )
	      If( choice_Comm .eq. fci_Reserved( 2 ) ) Then
c                                                        	! HELP
		Empty = .True.
		n_Buffer = 0
                if( from_file ) then
                  close( file_unit )
                  file_unit = -99
                  from_file = .false.
                end if
	      End If
	      Return
	    Else If( Short( 1: l_Typed )
     &                        .ne. List( I + 1 )( 1: l_Typed ) ) Then
	      choice_Comm = Order( I )
	      If( choice_Comm .eq. fci_Reserved( 2 ) ) Then
c                                                              	! HELP
		Empty = .True.
                if( from_file ) then
                  close( file_unit )
                  file_unit = -99
                  from_file = .false.
                end if
		n_Buffer = 0
	      End If
	      Return
	    Else
	      Write( pw_strng, 310 ) Short( 1: l_Typed )
310           Format( 'FCI ERROR:: Ambiguous command: "', A, '"' )
              call PWRITE( pw_strng )
	      Empty = .True.
	      n_Buffer = 0
              if( from_file ) then
                close( file_unit )
                file_unit = -99
                from_file = .false.
              end if
	      Return
	    End If
	  End If
	End Do
	Write( pw_strng, 320 ) Short( 1: l_Typed )
320     Format( 'FCI ERROR:: Unknown command :"', A, '"' )
        call PWRITE( pw_strng )
	Empty = .True.
	n_Buffer = 0
        if( from_file ) then
          close( file_unit )
          file_unit = -99
          from_file = .false.
        end if
	Return
c
390     continue
c       Come here if you've been reading from a file and detected End-of-file
        close( file_unit )
        file_unit = -99
        from_file = .false.
        empty = .true.
        goto 1000

C
400	Continue
c       Come here if you've been reading from STDIN and detected End-of-file
c                        		! ^Z=STOP
	choice_Comm = fci_Reserved( 1 )
C
	End
