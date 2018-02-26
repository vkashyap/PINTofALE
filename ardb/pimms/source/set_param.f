	Subroutine SET_PARAM( Line, Beg, Fin )
C
	Implicit None

        include 'xcomm.inc'
C
	Character*( * ) Line
	Integer Beg, Fin
C
	Integer I, asci_Code, Last, K, L
	Logical Dot
C					! First check for integers
	Dot = .False.
	I = Fin
	asci_Code = Ichar( Line( I: I ) )
	If( asci_Code .eq. 46 ) Goto 100
	If( asci_Code .lt. 48 .or. asci_Code .gt. 57 ) Goto 300
	I = I - 1
	Do While( I .gt. Beg )
	  asci_Code = Ichar( Line( I: I ) )
	  If( asci_Code .eq. 46 ) Then
	   Goto 100
	  Else If( asci_Code .eq. 45 .or. asci_Code .eq. 43 ) Then
	    I = I - 1
	    asci_Code = Ichar( Line( I: I ) )
	    If( asci_Code .eq. 69 .or. asci_Code .eq. 101 ) Then
	      Goto 200
	    Else
	      Goto 300
	    End If
	  Else If( asci_Code .eq. 69 .or. asci_Code .eq. 101 ) Then
	    Goto 200
	  Else If( asci_Code .lt. 48 .or. asci_Code .gt. 57 ) Then
	    Goto 300
	  End If
	  I = I - 1
	End Do
	If( I .lt. Beg ) Then
	  np_I = np_I + 1
	  If( np_I .gt. m_Par ) Then
	    Print *,
     &      'FCI WARNING:: Too many integer parameters...ignoring one'
	    np_I = np_I - 1
	  Else
	    Read( Line( Beg: Fin ), * ) p_Int( np_I )
	    np_N = np_N + 1
	    n_List( np_N ) = np_I + 100
            np_total = np_total + 1
            p_kind( np_total ) = 3
            po_i( np_i ) = np_total
            po_n( np_n ) = np_total
	  End If
	  Return
	End If
	asci_Code = Ichar( Line( I: I ) )
	If( asci_Code .eq. 46 ) Then
	  np_F = np_F + 1
	  If( np_F .gt. m_Par ) Then
	    Print *,
     &'FCI WARNING:: Too many floating-number parameters...ignoring one'
	    np_F = np_F - 1
	  Else
	    Read( Line( Beg: Fin ), * ) p_Float( np_F )
	    np_N = np_N + 1
	    n_List( np_N ) = np_F
            np_total = np_total + 1
            p_kind( np_total ) = 2
            po_f( np_i ) = np_total
            po_n( np_n ) = np_total
	  End If
	  Return
	Else If( asci_Code .eq. 45 .or. asci_Code .eq. 43
     1	  .or. ( asci_Code .ge. 48 .and. asci_Code .le. 57 ) ) Then
	  np_I = np_I + 1
	  If( np_I .gt. m_Par ) Then
	    Print *,
     &    'FCI WARNING:: Too many integer parameters...ignoring one'
	    np_I = np_I - 1
	  Else
	    Read( Line( Beg: Fin ), * ) p_Int( np_I )
	    np_N = np_N + 1
	    n_List( np_N ) = np_I + 100
            np_total = np_total + 1
            p_kind( np_total ) = 3
            po_i( np_i ) = np_total
            po_n( np_n ) = np_total
	  End If
	  Return
	Else
	  Goto 300
	End If
C
100	Continue		! Check for F format
	I = I - 1
	Do While( I .gt. Beg )
	  asci_Code = Ichar( Line( I: I ) )
	  If( asci_Code .lt. 48 .or. asci_Code .gt. 57 ) Goto 300
	  I = I - 1
	End Do
	asci_Code = Ichar( Line( I: I ) )
	If( asci_Code .eq. 45 .or. asci_Code .eq. 43
     1	  .or. ( asci_Code .ge. 48 .and. asci_Code .le. 57 ) ) Then
	  np_F = np_F + 1
	  If( np_F .gt. m_Par ) Then
	    Print *,
     &'FCI WARNING:: Too many floating-number parameters...ignoring one'
	    np_F = np_F - 1
	  Else
	    Read( Line( Beg: Fin ), * ) p_Float( np_F )
	    np_N = np_N + 1
	    n_List( np_N ) = np_F
            np_total = np_total + 1
            p_kind( np_total ) = 2
            po_f( np_i ) = np_total
            po_n( np_n ) = np_total
	  End If
	  Return
	Else
	  Goto 300
	End If
C
200	Continue			! Check for E format
	I = I - 1
	Do While( I .gt. Beg )
	  asci_Code = Ichar( Line( I: I ) )
	  If( asci_Code .eq. 46 ) Then
	    If( Dot ) Then
	      Goto 300
	    Else
	      Dot = .True.
	    End If
	  Else If( asci_Code .lt. 48 .or. asci_Code .gt. 57 ) Then
	    Goto 300
	  End If
	  I = I - 1
	End Do
	asci_Code = Ichar( Line( I: I ) )
	If( asci_Code .eq. 46 ) Then
	  If( Dot ) Goto 300
	  np_F = np_F + 1
	  If( np_F .gt. m_Par ) Then
	    Print *,
     &'FCI WARNING:: Too many floating-number parameters...ignoring one'
	    np_F = np_F - 1
	  Else
	    Read( Line( Beg: Fin ), * ) p_Float( np_F )
	    np_N = np_N + 1
	    n_List( np_N ) = np_F
            np_total = np_total + 1
            p_kind( np_total ) = 2
            po_f( np_i ) = np_total
            po_n( np_n ) = np_total
	  End If
	  Return
	Else If( asci_Code .eq. 45 .or. asci_Code .eq. 43
     1	  .or. ( asci_Code .ge. 48 .and. asci_Code .le. 57 ) ) Then
	  np_F = np_F + 1
	  If( np_F .gt. m_Par ) Then
	    Print *,
     &'FCI WARNING:: Too many floating-number parameters...ignoring one'
	    np_F = np_F - 1
	  Else
	    Read( Line( Beg: Fin ), * ) p_Float( np_F )
	    np_N = np_N + 1
	    n_List( np_N ) = np_F
            np_total = np_total + 1
            p_kind( np_total ) = 2
            po_f( np_i ) = np_total
            po_n( np_n ) = np_total
	  End If
	  Return
	Else
	  Goto 300
	End If
C
300	Continue			! Read it as a character string
	np_C = np_C + 1
	If( np_C .gt. m_Par ) Then
	  Print *, ' FCI WARNING:: '
     &         // 'Too many character string parameters...ignoring one'
	  np_C = np_C - 1
	Else
	  If( np_C .eq. 1 ) Then
	    Last = 0
	  Else
	    Last = pc_F( np_C - 1 )
	  End If
	  pc_B( np_C ) = Last + 1
	  K = Beg
	  L = Last + 1
	  If( Line( K: K ) .eq. ''''
     &                .and. Line( K: K + 1 ) .ne. '''''' ) Then
	    K = K + 1
	    Do While( K .lt. Fin .and. L .le. l_Buffer )
	      p_Char( L: L ) = Line( K: K )
	      If( p_Char( L: L ) .eq. '''' ) K = K + 1
	      K = K + 1
	      L = L + 1
	    End Do
	    If( L .gt. l_Buffer .and. K .ne. Fin ) Then
	      Print *,
     &        'FCI WARNING:: Character string too long...truncating'
	      pc_F( np_C ) = l_Buffer
	    Else
	      pc_F( np_C ) = L - 1
	    End If
	  Else
	    If( Fin - Beg .ge. l_Buffer - Last ) Then
	      Print *,
     &        'FCI WARNING:: Character string is too long...truncating'
	      pc_F( np_C ) = l_Buffer
	      p_Char( pc_B( np_C ): pc_F( np_C ) )
     1				= Line( Beg: Beg + l_Buffer - Last - 1 )
	    Else
	      pc_B( np_C ) = Last + 1
	      pc_F( np_C ) = Fin - Beg + Last + 1
	      p_Char( pc_B( np_C ): pc_F( np_C ) ) = Line( Beg: Fin )
	    End If
	  End If
          np_total = np_total + 1
          p_kind( np_total ) = 1
          po_c( np_c ) = np_total
	End If
C
	End
C
C
C
	Subroutine CLR_PARAM( )
C
	Implicit None

        include 'xcomm.inc'
C
	Integer I
C
	p_Char = ' '
	Do I = 1, m_Par
	  pc_B( I ) = 0
	  pc_F( I ) = 0
	  p_Float( I ) = 0.0
	  p_Int( I ) = 0
          po_c( i ) = 0
          po_f( i ) = 0
          po_i( i ) = 0
          po_n( i ) = 0
	End Do
	Do I = 1, m_Par * 2
	  n_List( I ) = 0
	End Do
        do i = 1, m_par * 3
          p_kind( i ) = 0
        end do
	np_C = 0
	np_F = 0
	np_I = 0
	np_N = 0
        np_total = 0
C
	End
