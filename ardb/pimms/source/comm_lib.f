	Subroutine CLRCOM( )
C
	Implicit None

        include 'xcomm.inc'
C
	Empty = .True.
	n_Buffer = 0
C
	End
C
C
C
	Integer Function FCI_LEN( String )
C
C	Find outs the length of a string
C
	Implicit None
C
	Character*( * ) String
C
	Integer Length, K, L
C
	Length = Len( String )
	Do K = 1, Length
	  If( String( K: K ) .gt. ' ' ) L = K
	End Do
	FCI_LEN = L
C
	End
C
C
C
	Subroutine FCI_UPPER_CASE( String, Length )
C
	Implicit none
C
	Integer Length
	Character String( Length )
C
	Integer I, Element
	Logical quote_On
C
	I = 1
	quote_On = .False.
	Do While( I .le. Length )
	  Element = Ichar( String( I ) )
	  If( quote_On ) Then			! Within a quote
	    If( Element .eq. 39 ) Then
	      I = I + 1
	      If( I .le. Length ) Then
		Element = Ichar( String( I ) )
		If( Element .ne. 39 ) Then	! Wasn't ''
		  quote_On = .False.
		  I = I - 1
		End If
	      End If
	    End If
	  Else
	    If( Element .eq. 39 ) Then		! Found a single quote
	      quote_On = .True.
	    Else If( Element .ge. 97 .and. Element .le. 122 ) Then
	      String( I ) = Char( Element - 32 )
	    Else If( Element .lt. 32 ) Then
	      String( I ) = ' '
	    End If
	  End If
	  I = I + 1
	End Do
C
	End
C
C
C
	Subroutine FCI_NOSPC( Line, I )
C
	Implicit None
C
	Character*( * ) Line
	Integer I, Length
C
	Length = Len( Line )
	Do While( I .le. Length )
	  If( Line( I: I ) .ne. ' ' ) Return
	  I = I + 1
	End Do
C
	End
