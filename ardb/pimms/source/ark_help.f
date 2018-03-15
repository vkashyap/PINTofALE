	Subroutine ARK_HELP( dir_name, lib_Name, String )
C
C	A brand new ARK routine to simulate VAX/VMS on-line help
C	on UNIX systems.  Requires an ark help library (pre-processed
C	by AHMAKE.FOR) and takes a string as the starting position.
C
C       Updated using Judy Chen (CfA)'s modifications re initializations

	Implicit None
C
	Integer m_Level, m_Item
	Parameter( m_Level = 16 )
	Parameter( m_Item = 256 )
C
	Character*( * ) dir_name, lib_Name, String
C
	Character*256 dummy_Line
	Character*80 Msg, one_Line
	Character*30 to_Find( m_Level ), c_Item( m_Item ), Temp
	Character*30 s_Title( m_Level ), blank_Line, Work
	Integer l_Find( m_Level ), b_Up( m_Level ), l_Title( m_Level )
	Integer l_Item( m_Item ), p_Item( m_Item ), r_Found( m_Item )
	Integer Head, Body, Level, n_Req, Req, K, Pointer, Item, L
	Integer n_Item, n_Line, Found, Ln, l_Page, F, l_Msg, Upto
        integer ierr, ii
	Logical not_Found, Ambiguous, gone_Up
C
	blank_Line = '                                '

	do ii = 1, m_level
          s_title( ii ) = blank_line
	end do

        call ARKOPN( k, dir_name, lib_name, 'AHL', 'OLD',
     &               'READONLY', 'UNFORMATTED', 'DIRECT', 80, ierr )
        if( ierr .lt. 0 ) then
          call PWRITE( 'SEVERE ERROR:: Failed to open HELP file' )
          stop
        end if

	Level = 1
	Head = 1
	dummy_Line = String
	Do While( .True. )
          Found = 0

	  Call AH_PROC_REQ( dummy_Line, to_Find, l_Find,
     1						m_Level, n_Req )
	  not_Found = .False.
	  Ambiguous = .False.
	  Call AH_GET_HEADER( K, Head, n_Item, n_Line, b_Up( Level ),
     1				Body, c_Item, l_Item, p_Item, m_Item )
	  Do Req = 1, n_Req
	    If( n_Item .eq. 0 ) Then
	      Head = b_Up( Level )
	      Level = Level - 1
	      Call AH_GET_HEADER( K, Head, n_Item, n_Line,
     1					b_Up( Level ), Body,
     2				c_Item, l_Item, p_Item, m_Item )
	      not_Found = .True.
	      Goto 110
	    End If
	    Found = 0
	    Do Item = 1, n_Item
	      Temp = c_Item( Item )
	      Call FCI_UPPER_CASE( Temp, l_Item( Item ) )
	      Ln = l_Find( Req )
	      Work = to_Find( Req )( 1: Ln )
	      Call FCI_UPPER_CASE( Work, Ln )
	      If( Work( 1: Ln ) .eq. Temp( 1: Ln ) ) Then
		Found = Found + 1
		r_Found( Found ) = Item
	      End If
	    End Do
	    If( Req .lt. n_Req .and. Found .gt. 1 ) Then
	      Ambiguous = .True.
	      Goto 110
	    Else If( Found .eq. 0 ) Then
	      not_Found = .True.
	      Goto 110
	    End If
	    Level = Level + 1
	    s_Title( Level ) = c_Item( r_Found( 1 ) )
	    l_Title( Level ) = l_Item( r_Found( 1 ) )
	    Head = p_Item( r_Found( 1 ) )
	    Call AH_GET_HEADER( K, Head, n_Item, n_Line, b_Up( Level ),
     1				Body, c_Item, l_Item, p_Item, m_Item )
	  End Do
110	  Continue
	  If( Ambiguous ) Then
            call PWRITE( to_find( req )( : ln ) // ' is ambiguous' )
	    l_Page = 1
	  Else If( not_Found ) Then
	    Ln = l_Find( Req )
            call PWRITE( to_find( req )( : ln ) // ' unknown' )
	    l_Page = 1
	  Else
C		Now display the correct info
	    l_Page = 0
	    Do F = 1, Found
	      Do L = 1, Level
		If( L .eq. 1 ) Then
c                 Space added by TTN
                  call PWRITE( s_title( l ) )
		Else
                  call PWRITE( blank_line( : l - 1 ) // s_title( l ) )
		End If
		Call AH_MORE( l_Page, dummy_Line, *200, *900 )
	      End Do
              call PWRITE( ' ' )
	      Call AH_MORE( l_Page, dummy_Line, *200, *900 )
	      Pointer = Body
	      Do Pointer = Body, Body + n_Line - 1
		Read( K, Rec = Pointer ) one_Line
c               Space added by TTN
                call PWRITE( one_line )
	      End Do
	      If( F .lt. Found ) Then
c		The folloing two lines inserted to correct a minor bug
		s_Title( Level ) = c_Item( r_Found( F + 1 ) )
		l_Title( Level ) = l_Item( r_Found( F + 1 ) )
		Head = p_Item( r_Found( F + 1 ) )
		Call AH_GET_HEADER( K, Head, n_Item, n_Line,
     1					b_Up( Level ), Body,
     2				c_Item, l_Item, p_Item, m_Item )
                call PWRITE( ' ' )
		Call AH_MORE( l_Page, dummy_Line, *200, *900 )
	      End If
	    End Do
	    If( Found .gt. 1 .or. n_Item .eq. 0 ) Then
	      Head = b_Up( Level )
	      Level = Level - 1
	      Call AH_GET_HEADER( K, Head, n_Item, n_Line,
     1					b_Up( Level ), Body,
     2				c_Item, l_Item, p_Item, m_Item )
	      gone_Up = .True.
	    Else
	      gone_Up = .False.
	    End If
	  End If
          call PWRITE( ' ' )
	  Call AH_MORE( l_Page, dummy_Line, *200, *900 )
c	  If( gone_Up ) Then
c	    Write( *, 151 )
c151	    Format( 'Information available on:' )
c	  Else
	  If( .not. gone_Up ) Then
            call PWRITE( 'More information available on:' )
c	  End If
	    Call AH_MORE( l_Page, dummy_Line, *200, *900 )
            call PWRITE( ' ' )
	    Call AH_MORE( l_Page, dummy_Line, *200, *900 )
	    Msg = ' '
	    l_Msg = 0
	    Do Item = 1, n_Item
	      Upto = ( l_Item( Item ) / 8 + 1 ) * 8
	      If( l_Msg + Upto .gt. 80 ) Then
c               Space added by TTN.
                call PWRITE( msg )
		Call AH_MORE( l_Page, dummy_Line, *200, *900 )
		Msg = ' '
		l_Msg = 0
	      End If
	      Msg( l_Msg + 1: l_Msg + Upto ) = c_Item( Item )
	      l_Msg = l_Msg + Upto
	    End Do
	    If( l_Msg .gt. 0 ) Then
c             Space added by TTN
              call PWRITE( msg )
	      Call AH_MORE( l_Page, dummy_Line, *200, *900 )
	    End If
            call PWRITE( ' ' )
	   Call AH_MORE( l_Page, dummy_Line, *200, *900 )
	  End If
	  Msg = ' '
	  l_Msg = 0
C	  BUG FIX by Koji, 1992 Mar 19 --- extended the special treatment of
C			Level=1 case, started DO loop from L=2
	  If( Level .eq. 1 ) Then
	    Msg = 'Topic ?'
	    l_Msg = 7
	  Else
	    Do L = 2, Level
	      If( l_Msg .eq. 0 ) Then
		l_Msg = l_Title( L )
		Msg = s_Title( L )( 1: l_Title( L ) )
	      Else
		If( l_Msg + l_Title( L ) + 1 .gt. 80 ) Then
c                 Space added by TTN
                  call PWRITE( msg( : l_msg ) )
		  Msg = ' '
		  Msg = s_Title( L )( 1: l_Title( L ) )
		  l_Msg = l_TItle( L )
		Else
		  Msg = Msg( 1: l_Msg ) // ' '
     1				// s_Title( L )( 1: l_Title( L ) )
		  l_Msg = l_Msg + 1 + l_Title( L )
		End If
	      End If
	    End Do
	    If( l_Msg + 12 .gt. 80 ) Then
c             Space added by TTN
              call PWRITE( msg( : l_msg ) )
	      Msg = 'Sub-Topic ?'
	      l_Msg = 11
	    Else
	      Msg = Msg( 1: l_Msg ) // ' Sub-Topic ?'
	      l_Msg = l_Msg + 12
	    End If
	  End If
	  Call WRITEN( Msg( 1: l_Msg ) )
	  Read( *, '(A)', End = 900 ) dummy_Line
200	  Continue
	  Do While( dummy_Line .eq. ' ' )
	    Head = b_Up( Level )
	    Level = Level - 1
	    If( Head .eq. -1 .or. Level .eq. 0 ) Goto 910
	    If( Level .gt. 1 ) Then
	      l_Msg = l_Title( Level )
	      Msg = s_Title( Level )( 1: l_Msg ) // ' Sub-Topic ?'
	      l_Msg = l_Msg + 12
	      Call WRITEN( Msg( 1: l_Msg ) )
	    Else
	      Call WRITEN( 'Topic ?' )
	    End If
	    Read( *, '(A)', End = 900 ) dummy_Line
	  End Do
	  If( dummy_Line .eq. '?' ) Then
c	    If( .not. gone_Up ) Then
c	      Head = b_Up( Level )
c	      Level = Level - 1
c	      Call AH_GET_HEADER( K, Head, n_Item, n_Line,
c     1					b_Up( Level ), Body,
c     2				c_Item, l_Item, p_Item, m_Item )
c	    End If
	    dummy_Line = ' '
C	    BUGFIX attempt by Koji...seems to work
C	    (without this, crashed at ? if last Found > 1)
	    Found = 1
	  End If
	End Do
C
900	Continue
        call PWRITE( ' ' )
910	Continue
	Close( K )
C
	End
C
C
C
	Subroutine AH_PROC_REQ
     1			( dummy_Line, to_Find, l_Find, m_Level, n_Req )
C
C	Decomposes the user input into distinctive words, places them
C	in to_Find array (length of each word in l_Find )
C
	Implicit None
C
	Character*( * ) dummy_Line
	Integer m_Level
	Character*30 to_Find( m_Level )
	Integer l_Find( m_Level ), n_Req
C
	Integer l_Line, Beg, Fin, One
C
	l_Line = Len( dummy_Line )
	Beg = 1
	Call FCI_NOSPC( dummy_Line, Beg )
	n_Req = 0
	Do While( Beg .le. l_Line )
	  n_Req = n_Req + 1
	  Fin = Beg
	  One = Ichar( dummy_Line( Fin: Fin ) )
	  Do While( One .gt. 32 .and. Fin .lt. l_Line )
	    Fin = Fin + 1
	    One = Ichar( dummy_Line( Fin: Fin ) )
	  End Do
	  If( One .le. 32 ) Fin = Fin - 1
	  l_Find( n_Req ) = Fin - Beg + 1
	  If( l_Find( n_Req ) .gt. 30 ) Then
c           Changed by TTN to free format.  Really long keywords overran
c           the statement.
            call PWRITE( 'ARKHELP ERROR:: keyword too long' )
            call PWRITE( '...truncating to 30 characters' )
	    to_Find( n_Req ) = dummy_Line( Beg: Beg + 31 )
	    l_Find( n_Req ) = 30
	  Else
	    to_Find( n_Req ) = dummy_Line( Beg: Fin )
	  End If
	  Beg = Fin + 1
	  If( Beg .le. l_Line ) Call FCI_NOSPC( dummy_Line, Beg )
	End Do
C
	End
C
C
C
	Subroutine AH_GET_HEADER( Unit, Head, n_Item, n_Line, Up, Body,
     1					c_Item, l_Item, p_Item, m_Item )
C
C	Reads the Header information from a ARK Help Library Entry
C
	Implicit None
C
	Integer Unit, Head, n_Item, n_Line, Up, m_Item
	Character*30 c_Item( m_Item )
	Integer l_Item( m_Item ), p_Item( m_Item )
C
	Character*80 header_Line
	Integer Body, Item
c       Body was declared twice, cut to once by TTN
C
	Read( Unit, Rec = Head ) header_Line
	Read( header_Line, * ) n_Item, n_Line, Up
	Body = Head
	If( n_Item .gt. m_Item )
     1		Stop 'EMERGENCY STOP:: Buffer overflow in ARK_HELP'
	Do Item = 1, n_Item, 2
	  Body = Body + 1
	  Read( Unit, Rec = Body ) header_Line
	  c_Item( Item ) = header_Line( 1: 30 )
c         Format changed to I5 by TTN.
	  Read( header_Line( 31: 35 ), '(I5)' ) l_Item( Item )
	  Read( header_Line( 36: 40 ), '(I5)' ) p_Item( Item )
	  If( Item .lt. n_Item ) Then
	    c_Item( Item + 1 ) = header_Line( 41: 70 )
c           Format changed to I5 by TTN.
	    Read( header_Line( 71: 75 ), '(I5)' ) l_Item( Item + 1 )
	    Read( header_Line( 76: 80 ), '(I5)' ) p_Item( Item + 1 )
	  End If
	End Do
	Body = Body + 1
C
	End
C
C
C
	Subroutine AH_MORE( l_Page, dummy_Line, *, * )
C
	Implicit None
C
	Integer l_Page
	Character*( * ) dummy_Line
C
	l_Page = l_Page + 1
	If( l_Page .eq. 20 ) Then
	  Call WRITEN( 'Type <CR> for more...' )
	  Read( *, '(A)', End = 900 ) dummy_Line
          call PWRITE( ' ' )
	  If( dummy_Line .ne. ' ' ) Return 1
	  l_Page = 0
	End If
	Return
C
900	Continue
	Return 2
C
	End
