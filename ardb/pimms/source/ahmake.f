	Program AHMAKE
C
C	This is stage 1 of VMS help emulator.
C	AHMAKE takes a text file, the kind that can also be passed to
C	$ LIBRARY/HELP/CREATE=KEYSIZE:31
C	and creates an ARK help library file (equivalent of .HLB file).
C
C	This is a version with hardwired file names, for PIMMS
C
C       Added initialization statements for "where" and other arrays
C       2011 Sep
C
	Implicit None
C
	Integer m_Total, m_Level
	Parameter( m_Total = 1024 )
	Parameter( m_Level = 16 )
C
	Integer Where( m_Total ), l_Parent( m_Level ), Up( m_Total )
	Integer l_Zero( m_Total ), n_Item( m_Total ), Index( m_Total )
	Integer l_Title( m_Total ), n_Line( m_Total ), p_Title( m_Total )
	Integer l_Up( m_Level )
	Integer Level, new_Level, Item, Line, One, Beg, Fin
	Integer Parent, Zero, Pointer, I, L, N
	Character*30 Title( m_Total )
	Character*80 one_Line, out_Line, blank_Line
        integer lun1, lun7, flag

        data where / m_total * 0 /
        data l_parent  / m_level * 0 /
        data up / m_total * 0 /
        data l_zero / m_total * 0 /
        data n_item / m_total * 0 /
        data index / m_total * 0 /
        data l_title / m_total * 0 /
        data n_line / m_total * 0 /
        data p_title / m_total * 0 /
        data l_up / m_level * 0 /
        
C
	blank_Line( 1: 40 ) = '                                        '
	blank_Line( 41: 80 ) = blank_Line( 1: 40 )

        call ARKOPN( lun1, ' ', 'pimms.hlp', 'hlp',
     &               'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &               1, flag )

        call ARKOPN( lun7, ' ', 'pimms.ahl', 'ahl',
     &               'NEW', 'OVERWRITE', 'UNFORMATTED', 'DIRECT',
     &                80, flag )


	Level = 0
	Item = 1
	Index( 1 ) = 0
	Line = 0
	l_Up( 1 ) = 1
	Up( 1 ) = -1
	p_Title( 1 ) = 1
	l_Zero( 1 ) = 0
	l_Parent( 1 ) = 1
C		FIRST PASS
	Do While( .True. )
	  Read( lun1, '(A)', End = 500 ) one_Line
	  Line = Line + 1
	  One = Ichar( one_Line( 1: 1 ) )
	  If( One .ge. 49 .and. One .le. 57 ) Then
c							! Found a new HELP item
	    Fin = 1
	    Do While( One .ge. 48 .and. One .le. 57 )
	      Fin = Fin + 1
	      One = Ichar( one_Line( Fin: Fin ) )
	    End Do
	    Read( one_Line( 1: Fin - 1 ), * ) new_Level
c		! New item is at level = new_level
	    Beg = Fin
c           "ichar" inserted by TTN.
	    One = ichar(one_Line( Beg: Beg ))
	    Do While( One .le. 32 )
	      Beg = Beg + 1
c             "ichar" inserted by TTN.
	      One = ichar(one_Line( Beg: Beg ))
	    End Do
	    Fin = Beg
	    Do While( One .gt. 32 )
	      Fin = Fin + 1
c             "ichar" inserted by TTN.
	      One = ichar(one_Line( Fin: Fin ))
	    End Do
	    Fin = Fin - 1
c					! Update INDEX
	    Item = Item + 1
	    Where( Item ) = Line
	    p_Title( Item ) = p_Title( Item - 1 ) + n_Line( Item - 1 )
     1				+ ( n_Item( Item - 1 ) + 1 ) / 2 + 1
c		It's the Item'th item overall, at line Line, the header
c		line should eventual end up on line p_Title( Item )
	    If( new_Level .le. Level ) Then
c		Okay, the previous item was the end of a search branch
	      Parent = l_Parent( new_Level )
c		The new item is a branch of the Parent's item
	      Zero = l_Zero( Parent ) + n_Item( Parent ) + 1
	      Do I = Item, Zero + 1, -1
		Index( I ) = Index( I - 1 )
	      End Do
	      Index( Zero ) = Item
c		All branches of the parent a sorted in order
	      Do I = 1, Item
		If( l_Zero( I ) .ge. Zero - 1 ) Then
		  l_Zero( I ) = l_Zero( I ) + 1
		End If
	      End Do
	    Else If( new_Level .eq. Level + 1 ) Then
c		This item is the branch of the item immediately preceeding
	      Parent = Item - 1
	      Index( Parent ) = Item
	    Else
	      Stop 'EMERGENCY STOP:: Illegal jump in item level'
	    End If
	    l_Parent( new_Level + 1 ) = Item
	    Up( Item ) = l_Up( new_Level )
	    Level = new_Level
	    n_Item( Parent ) = n_Item( Parent ) + 1
	    If( Mod( n_Item( Parent ), 2 ) .eq. 1 ) Then
	      Do I = Parent + 1, Item
		p_Title( I ) = p_Title( I ) + 1
	      End Do
	      Do I = 1, Item
		If( Up( I ) .gt. p_Title( Parent ) ) Then
		  Up( I ) = Up( I ) + 1
		End If
	      End Do
	    End If
	    l_Up( new_Level + 1 ) = p_Title( Item )
	    Title( Item ) = one_Line( Beg: Fin )
	    l_Title( Item ) = Fin - Beg + 1
	    l_Zero( Item ) = Item - 1
	    n_Item( Item ) = 0
	    n_Line( Item ) = 0
	  Else
c							! HELP text
	    n_Line( Item ) = n_Line( Item ) + 1
	  End If
	End Do
C
500	Continue
	Rewind( lun1 )
C		SECOND PASS
	Item = 1
	Pointer = 1
	Do L = 0, Line
	  If( L .ge. 1 ) Read( lun1, '(A)' ) one_Line
	  If( L .eq. 0 .or. L .eq. Where( Item ) ) Then
c		Found a header item
	    Write( out_Line, * )
     1			n_Item( Item ), n_Line( Item ), Up( Item )
	    Write( lun7, Rec = Pointer ) out_Line
	    Pointer = Pointer + 1
	    Zero = l_Zero( Item )
	    Do N = 1, n_Item( Item ), 2
	      I = Index( Zero + N )
	      out_Line( 1: 30 ) = Title( I )
c             Format changed to I5 by TTN.
	      Write( out_Line( 31: 40 ), '(2I5)' )
     1					l_Title( I ), p_Title( I )
	      If( N .lt. n_Item( Item ) ) Then
		I = Index( Zero + N + 1 )
		out_Line( 41: 70 ) = Title( I )
c               Format changed to I5 by TTN.
		Write( out_Line( 71: 80 ), '(2I5)' )
     1					l_Title( I ), p_Title( I )
	      Else
		out_Line( 41: 80 ) = ' '
	      End If
	      Write( lun7, Rec = Pointer ) out_Line
	      Pointer = Pointer + 1
	    End Do
	    Item = Item + 1
	  Else
c		Write out with a suitable number of leading space
	    If( Level .eq. 0 ) Then
	      out_Line = one_Line
	    Else
	      out_Line = blank_Line( 1: Level )
     1					// one_Line( 1: 80 - Level )
	    End If
	    Write( lun7, Rec = Pointer ) out_Line
	    Pointer = Pointer + 1
	  End If
	End Do
C
	End
