*+PIN_LIMIT
        subroutine PIN_LIMIT( results, n_res )

        implicit none

        integer n_res
        real results( 0: n_res )

*       Description:
*         Outputs Suzaku HXD/PIN specific information to the screen via PWRITE
*
*       Arguments:
*         results (i)  : Predicted count rates
*         n_res   (i)  : Number of standard bands
*
*       Dependencies:
*         PWRITE
*
*       Origin:
*         Created by KM
*
*       Author
*         Koji Mukai, 2008 Oct
*         Modified 2008 Nov: now accounts for the statistical errors
*            of the background model as well as source (with 10x over
*            sampling)
*         Updated 2010 Oct using new AO-6 background numbers
*-PIN_LIMIT

        character*76 xw_strng
        real temp, labour, work, threesigma, fivesigma, bgd11
        real frac, bgd_rate, cxb_rate

        data frac, bgd_rate, cxb_rate / 0.013, 0.2783, 0.0018 /

        xw_strng = '% In the 15-40 keV band, PIMMS predicts'
        call PWRITE( xw_strng )
        xw_strng = '  n.nnnE+mm source cps, n.nnnn CXB cps, and '
     &                                          // 'n.nnnn NXB cps'
        if( results( 1 ) .gt. 5000.0 ) then
          write( xw_strng( 3: 11 ), '(1p,e9.3)' ) results( 1 )
        else if( results( 1 ) .gt. 0.01 ) then
          write( xw_strng( 3: 11 ), '(f9.3)' ) results( 1 )
        else
          write( xw_strng( 3: 11 ), '(1p,e9.3)' ) results( 1 )
        end if
        write( xw_strng( 25: 30 ), '(f6.4)' ) cxb_rate
        write( xw_strng( 45: 50 ), '(f6.4)' ) bgd_rate
        call PWRITE( xw_strng )
        temp = results( 1 ) * results( 1 )
        labour = frac * bgd_rate
        bgd11 = cxb_rate + bgd_rate * 1.1
        if( temp .gt. 1.0e-30 ) then
          work = temp / 9.0 - labour * labour
          if( work .lt. 0.0 ) then
            xw_strng = '  The source is undetectable at 3-sigma'
          else
            threesigma = ( results( 1 ) + bgd11 ) / work
            xw_strng = '  The source is detectable at 3-sigma '
     &                                             // 'in n.nnnE+nn s'
            if( threesigma .gt. 5000.0 ) then
              write( xw_strng( 42: 50 ), '(1p,e9.3)' ) threesigma
            else if( threesigma .gt. 0.01 ) then
              write( xw_strng( 42: 50 ), '(f9.3)' ) threesigma
            else
              write( xw_strng( 42: 50 ), '(1p,e9.3)' ) threesigma
            end if
          end if
          call PWRITE( xw_strng )
          work = temp / 25.0 - labour * labour
          if( work .lt. 0.0 ) then
            xw_strng = '  The source is undetectable at 5-sigma'
          else
            fivesigma = ( results( 1 ) + bgd11 ) / work
            xw_strng = '  The source is detectable at 5-sigma '
     &                                             // 'in n.nnnE+nn s'
            if( fivesigma .gt. 5000.0 ) then
              write( xw_strng( 42: 50 ), '(1p,e9.3)' ) fivesigma
            else if( fivesigma .gt. 0.01 ) then
              write( xw_strng( 42: 50 ), '(f9.3)' ) fivesigma
            else
              write( xw_strng( 42: 50 ), '(1p,e9.3)' ) fivesigma
            end if
          end if
        else
          xw_strng = '  The source is undetectable at 3-sigma'
        end if
        call PWRITE( xw_strng )
        xw_strng = '  (considering the x.x% systematic '
     &                     // 'uncertainty in the NXB estimation)'
        write( xw_strng( 20: 22 ), '(f3.1)' ) frac * 100.0
        call PWRITE( xw_strng )

        end


*+GSO_LIMIT
        subroutine GSO_LIMIT( results, n_res )

        implicit none

        integer n_res
        real results( 0: n_res )

*       Description:
*         Outputs Suzaku HXD/GSO specific information to the screen via PWRITE
*
*       Arguments:
*         results (i)  : Predicted count rates
*         n_res   (i)  : Number of standard bands
*
*       Dependencies:
*         PWRITE
*
*       Origin:
*         Created by KM
*
*       Author
*         Koji Mukai, 2008 Oct
*         Modified 2008 Nov: now accounts for the statistical errors
*            of the background model as well as source
*         Updated 2010 Oct using new AO-6 background numbers
*-GSO_LIMIT

        character*76 xw_strng
        real temp, labour, work, threesigma, fivesigma
        real frac( 2 ), bgd_rate( 2 )
        character*7 range( 2 )
        integer k

        data frac, bgd_rate / 0.0064, 0.0059, 8.638, 12.03 /
        data range / ' 50-100', '100-200' /

        do k = 1, 2
          xw_strng = '% In the nnn-nnn keV band, PIMMS predicts'
          xw_strng( 10: 16 ) = range( k )
          call PWRITE( xw_strng )
          xw_strng = '  n.nnnE+mm source cps and nn.nnn NXB cps'
          if( results( k ) .gt. 5000.0 ) then
            write( xw_strng( 3: 11 ), '(1p,e9.3)' ) results( k )
          else if( results( k ) .gt. 0.01 ) then
            write( xw_strng( 3: 11 ), '(f9.3)' ) results( k )
          else
            write( xw_strng( 3: 11 ), '(1p,e9.3)' ) results( k )
          end if
          write( xw_strng( 28: 33 ), '(f6.3)' ) bgd_rate( k )
          call PWRITE( xw_strng )
          temp = results( k ) * results( k )
          labour = frac( k ) * bgd_rate( k )
          if( temp .gt. 1.0e-30 ) then
            work = temp / 9.0 - labour * labour
            if( work .lt. 0.0 ) then
              xw_strng = '  The source is undetectable at 3-sigma'
            else
              threesigma = ( results( k ) + bgd_rate( k ) * 2.0 ) / work
              xw_strng = '  The source is detectable at 3-sigma '
     &                                  // 'in n.nnnE+nn s'
              if( threesigma .gt. 5000.0 ) then
                write( xw_strng( 42: 50 ), '(1p,e9.3)' ) threesigma
              else if( threesigma .gt. 0.01 ) then
                write( xw_strng( 42: 50 ), '(f9.3)' ) threesigma
              else
                write( xw_strng( 42: 50 ), '(1p,e9.3)' ) threesigma
              end if
            end if
            call PWRITE( xw_strng )
            work = temp / 25.0 - labour * labour
            if( work .lt. 0.0 ) then
              xw_strng = '  The source is undetectable at 5-sigma'
            else
              fivesigma = ( results( k ) + bgd_rate( k ) * 2.0 ) / work
              xw_strng = '  The source is detectable at 5-sigma '
     &                                             // 'in n.nnnE+nn s'
            if( threesigma .gt. 5000.0 ) then
              write( xw_strng( 42: 50 ), '(1p,e9.3)' ) fivesigma
            else if( threesigma .gt. 0.01 ) then
              write( xw_strng( 42: 50 ), '(f9.3)' ) fivesigma
            else
              write( xw_strng( 42: 50 ), '(1p,e9.3)' ) fivesigma
            end if
          end if
          else
            xw_strng = '  The source is undetectable at 3-sigma'
          end if
          call PWRITE( xw_strng )
          xw_strng = '  (considering the x.xx% systematic '
     &                     // 'uncertainty in the NXB estimation)'
          write( xw_strng( 20: 23 ), '(f4.2)' ) frac( k ) * 100.0
          call PWRITE( xw_strng )
        end do

        end
