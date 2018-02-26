*+BAT_LIMIT
        subroutine BAT_LIMIT( results, n_res )

        implicit none

        real alpha, Ndet, snrat
        parameter( alpha = 0.28 )
        parameter( Ndet = 32768.0 )
        parameter( snrat = 5.0 )

        integer n_res
        real results( 0: n_res )

*       Description:
*         Outputs Swift BAT specific information to the screen via PWRITE
*
*       Arguments:
*         results (i)  : Predicted count rates
*         n_res   (i)  : Number of standard bands
*
*       Dependencies:
*         PCA_LIMIT_DO
*         PWRITE
*
*       Origin:
*         Created by KM based on algorithm by Craig Markwardt
*
*       Author
*         Koji Mukai, 2006 Apr
*-BAT_LIMIT

        integer k
        real work, temp, fivesigma
        character*78 pw_strng
        real brate( 0: 4 )
        character*20 bounds( 4 )
c        data brate / 7297.27143, 2133.71118, 2325.52686, 1800.08606,
c     &                                                      1037.94629 /
        data brate / 0.262011021, 0.0766116455, 0.0834988579,
     &                                      0.0646327361, 0.0372678265 /
        data bounds / ' 3-7 /  14-24 keV   ', ' 8-19/  24-48.9 keV ',
     &                '20-41/48.9-98.8 keV ', '42-78/98.8-194.9 keV' /

        if( n_res .ne. 4 ) then
          call PWRITE( 'SEVERE ERROR:: PIMMS is confused in BAT_LIMIT' )
          return
        end if
        if( results( 0 ) .gt. 1.0e3 ) then
          pw_strng = 'Predicted count rate very high: mistake maybe?'
          call PWRITE( pw_strng )
        else
          work = snrat * snrat * ( 2.0 * brate( 0 ) + results( 0 ) )
          temp = results( 0 ) * results( 0 ) * alpha * Ndet 
          fivesigma = work / temp
          pw_strng = '  With nnnn.n background cps, ' //
     &                    'it can be detected at S/N=abc in n.nnnE+mm s'
          write( pw_strng( 8: 13 ), '(f6.4)' ) brate( 0 )
          write( pw_strng( 57: 59 ), '(f3.1)' ) snrat
          if( fivesigma .gt. 5000.0 ) then
            write( pw_strng( 64: 72 ), '(1p,e9.3)' ) fivesigma
          else if( results( 0 ) .gt. 0.01 ) then
            write( pw_strng( 64: 72 ), '(f9.3)' ) fivesigma
          else
            write( pw_strng( 64: 72 ), '(1p,e9.3)' ) fivesigma
          end if
          call PWRITE( pw_strng )
          call PWRITE( ' ' )
          pw_strng =
     &           'Results in the 4 standard Swift/BAT survey bands are:'
          call PWRITE( pw_strng )
          call PWRITE( ' ' )
          pw_strng =
     &            '    Channel / Energy     Source     BGD    Detection'
          call PWRITE( pw_strng )
          do k = 1, 4
            pw_strng =
     &          '  1-3              keV  n.nnnE+mm  ffff.f  n.nnnE+mm s'
            pw_strng( 3: 22 ) = bounds( k )
            if( results( k ) .gt. 5000.0 ) then
              write( pw_strng( 25: 33 ), '(1p,e9.3)' ) results( k )
            else if( results( k ) .gt. 0.01 ) then
              write( pw_strng( 25: 33 ), '(f9.3)' ) results( k )
            else
              write( pw_strng( 25: 33 ), '(1p,e9.3)' ) results( k )
            end if
            write( pw_strng( 36: 41 ), '(f6.4)' ) brate( k )
            work = snrat * snrat * ( 2.0 * brate( k ) + results( k ) )
            temp = results( k ) * results( k ) * alpha * Ndet 
            fivesigma = work / temp
            if( fivesigma .gt. 5000.0 ) then
              write( pw_strng( 44: 52 ), '(1p,e9.3)' ) fivesigma
            else if( fivesigma .gt. 0.01 ) then
              write( pw_strng( 44: 52 ), '(f9.3)' ) fivesigma
            else
              write( pw_strng( 44: 52 ), '(1p,e9.3)' ) fivesigma
            end if
            call PWRITE( pw_strng )
          end do
        end if

        end
