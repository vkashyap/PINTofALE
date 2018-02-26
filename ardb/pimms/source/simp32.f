*+SIMP32

        real function SIMP32( FX, lower, upper, error )

        implicit none

        real FX
        external FX
        real lower, upper
        real error

*       Description:
*         Integrates external function FX from lower to uper
*         with more divisions of this range until it stops improving
*         and there are at least 32 divisions.
*
*       This routine is deprecated, here for historical reasons only.
*
*       Arguments:
*         FX        (i) : External reference to a one-argument function
*         lower     (i) : lower boundary for integration
*         upper     (i) : upper boundary for integration
*         error     (i) : fractional improvement to aim for
*         <SIMP32>  (r) : result of integration
*
*       Dependencies:
*         None
*
*       Origin:
*         A program by Graziella Branduardi=Raymont
*
*       Author:
*         Koji Mukai, 1993 March, first official PIMMS version
*-SIMP32

        real h, x, even, odd, edge, sum, save, diff, comp
        integer i, m, n

        SIMP32 = 0.0
        if( upper .le. lower ) then
          call PWRITE( 'ERROR in SIMP32:: reverse boundary' )
          return
        end if
        edge = FX( lower ) + FX( upper )
        n = 1
        even = 0.0
        h = upper - lower
        diff = abs( edge )
        comp = diff * error - 0.1
        sum = edge
        do while( diff .gt. comp .or. n .le. 16 )
          n = n * 2
          h = h * 0.5
          save = sum
          m = n - 1
          odd = 0.0
          do i = 1, m, 2
            x = lower + i * h
            odd = odd + FX( x )
          end do
          sum = h * ( edge + even * 2.0 + odd * 4.0 ) * 0.333333333
          even = even + odd
          diff = abs( sum - save )
          comp = abs( sum * error )
        end do

        SIMP32 = sum

        end


*+SIMP2K

        real function SIMP2K( FX, lower, upper, error )

        implicit none

        real FX
        external FX
        real lower, upper
        real error

*       Description:
*         Integrates external function FX from lower to uper
*         with more divisions of this range until it stops improving
*         and there are at least 2048 divisions.  Internal calculations
*         are mostly done in double precision in this version.
*
*       Arguments:
*         FX        (i) : External reference to a one-argument function
*         lower     (i) : lower boundary for integration
*         upper     (i) : upper boundary for integration
*         error     (i) : fractional improvement to aim for
*         <SIMP2K>  (r) : result of integration
*
*       Dependencies:
*         None
*
*       Origin:
*         A program by Graziella Branduardi-Raymont
*
*       Author:
*         Koji Mukai, 2000 April, copied SIMP32
*-SIMP2K

        double precision h, even, odd, edge, sum, save, diff, comp
        real x
        integer i, m, n

        SIMP2K = 0.0
        if( upper .le. lower ) then
          call PWRITE( 'ERROR in SIMP2K:: reverse boundary' )
          return
        end if
        edge = FX( lower ) + FX( upper )
        n = 1
        even = 0.0
        h = upper - lower
        diff = abs( edge )
        comp = diff * error - 0.1
        sum = edge
        do while( diff .gt. comp .or. n .le. 1024 )
          n = n * 2
          if( n .ge. 262144 ) then
            call PWRITE( 'SIMP2K Warning:: Taking too long to '
     &                              // 'converge, cutting at n=131072' )
            call PWRITE( '                 the result may be '
     &                                       // 'somewhat inaccurate.' )
            SIMP2K = sum
            return
          end if
          h = h * 0.5
          save = sum
          m = n - 1
          odd = 0.0
          do i = 1, m, 2
            x = lower + i * h
            odd = odd + FX( x )
          end do
          sum = h * ( edge + even * 2.0 + odd * 4.0 ) * 0.333333333
          even = even + odd
          diff = abs( sum - save )
          comp = abs( sum * error )
        end do

        SIMP2K = sum

        end
