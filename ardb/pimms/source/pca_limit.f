*+PCA_LIMIT
        subroutine PCA_LIMIT( results, n_res )

        implicit none

        real frac
        parameter( frac = 0.01 )

        integer n_res
        real results( 0: n_res )

*       Description:
*         Outputs XTE PCA specific information to the screen via PWRITE
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
*         Created by KM
*
*       Author
*         Koji Mukai, 1994 Dec
*-PCA_LIMIT

        character*64 xw_strng

        if( n_res .ne. 6 ) then
          call PWRITE( 'SEVERE ERROR:: PIMMS is confused in PCA_LIMIT' )
          return
        end if
        if( results( 0 ) .gt. 1.0e15 ) then
          xw_strng = 'Predicted count rate very high: mistake maybe?'
        else
          call PWRITE( ' ' )
          call PWRITE( '%%%        With 3 PCUs operational:' )
          call PWRITE( '           (Use these numbers in RPS)' )
          call PWRITE( ' ' )
          call PCA_LIMIT_DO( results, 3.0, frac )
          call PWRITE( ' ' )
          call PWRITE( '%%% ...and with 2 PCUs operational:' )
          call PWRITE( ' ' )
          call PCA_LIMIT_DO( results, 2.0, frac )
        end if

        end



*+PCA_LIMIT_DO
        subroutine PCA_LIMIT_DO( results, frac_PCU, frac_SYS )

        implicit none

        integer n_res_fixed
        parameter( n_res_fixed = 6 )

        real results( 0: n_res_fixed )
        real frac_PCU
        real frac_SYS

*       Description:
*         Sub-subroutine of PCA_LIMIT, which outputs XTE PCA
*         specific information to the screen via PWRITE
*
*       Arguments:
*         results (i)  : Predicted count rates
*         frac_PCU (i) : Fraction of PCUs operating
*         frac_SYS (i) : Systematic error of the background estimates
*
*       Dependencies:
*         PWRITE
*
*       Origin:
*         Created by KM
*
*       Author
*         Koji Mukai, 1999 June
*-PCA_LIMIT_DO

        character*74 xw_strng
        real temp, fivesigma, work, syst, fivesig_b
        real fr_res, fr_bgd
        integer j

        real bgd_rate( 0: n_res_fixed )
        character*20 bounds( n_res_fixed )
C The Epoch 4 Bgd rate per PCU
        data bgd_rate / 30.46, 3.52, 1.154, 1.52, 2.622, 2.966, 18.684 /
C The Epoch 4 Bgd rate for 5 PCUs
C        data bgd_rate / 152.3, 17.60, 5.77, 7.60, 13.11, 14.83, 93.42 /
C The Epoch 3 Bgd rate
C        data bgd_rate / 95.9, 10.5, 4.3, 5.3, 8.2, 8.7, 59.2 /
C        data bgd_rate / 109.5, 7.0, 3.5, 5.0, 6.5, 7.5, 80.0 /
C        data bgd_rate / 90.17, 5.84, 1.80, 2.76, 5.28, 5.29, 69.2 /
C These are the Epoch 4 boundaries
        data bounds / '  0- 13  0.00- 6.14', ' 14- 17  6.14- 7.90',
     &                ' 18- 23  7.90- 10.5', ' 24- 35  10.5- 15.8',
     &                ' 36- 49  15.8- 22.1', ' 50-249  22.1-116.0' /
C These are the Epoch 3 boundaries
C        data bounds / '  0- 13  0.00- 5.01', ' 14- 17  5.01- 6.47',
C     &                ' 18- 23  6.47- 8.68', ' 24- 35  8.68-13.12',
C     &                ' 36- 49 13.12-18.34', ' 50-249 18.34-98.52' /
C These are the Epoch 1 boundaries
C        data bounds / '  0- 13  0.00- 3.25', ' 14- 17  3.25- 4.25',
C     &                ' 18- 23  4.25- 5.75', ' 24- 35  5.75- 8.75',
C     &                ' 36- 49  8.75-12.25', ' 50-249 12.25-62.25' /

        xw_strng = 'PIMMS predicts n.nnnE+mm cps from the source '
     &                                // 'plus n.nnnE+mm background cps'
        fr_res = results( 0 ) * frac_PCU
        if( fr_res .gt. 5000.0 ) then
          write( xw_strng( 16: 24 ), '(1p,e9.2)' ) fr_res
        else if( fr_res .gt. 0.01 ) then
          write( xw_strng( 16: 24 ), '(f9.3)' ) fr_res
        else
          write( xw_strng( 16: 24 ), '(1p,e9.2)' ) fr_res
        end if
        fr_bgd = bgd_rate( 0 ) * frac_PCU
        if( fr_bgd .gt. 5000.0 ) then
          write( xw_strng( 51: 59 ), '(1p,e9.2)' ) fr_bgd
        else if( fr_bgd .gt. 0.01 ) then
          write( xw_strng( 51: 59 ), '(f9.3)' ) fr_bgd
        else
          write( xw_strng( 51: 59 ), '(1p,e9.2)' ) fr_bgd
        end if
        call PWRITE( xw_strng )
        temp = fr_res * fr_res
        if( temp .gt. 1e-30 ) then
          work = 25.0 * ( fr_res + fr_bgd )
          fivesigma = work / temp
          xw_strng = '5-sigma detection will be achieved in         s'
          if( fivesigma .gt. 5000.0 ) then
            write( xw_strng( 38: 46 ), '(1p,e9.2)' ) fivesigma
          else if( fivesigma .gt. 0.01 ) then
            write( xw_strng( 38: 46 ), '(f9.3)' ) fivesigma
          else
            write( xw_strng( 38: 46 ), '(1p,e9.2)' ) fivesigma
          end if
          call PWRITE( xw_strng )
          syst = frac_SYS * fr_bgd
          if( fr_res .gt. 5.0 * syst ) then
            fivesig_b = work / ( temp - 25.0 * syst * syst )
            xw_strng =
     &       '(or in         s with 1% systematic uncertainties in bgd)'
            if( fivesigma .gt. 5000.0 ) then
              write( xw_strng( 7: 15 ), '(1p,e9.2)' ) fivesig_b
            else if( fivesigma .gt. 0.01 ) then
              write( xw_strng( 7: 15 ), '(f9.3)' ) fivesig_b
            else
              write( xw_strng( 7: 15 ), '(1p,e9.2)' ) fivesig_b
            end if
          else
            xw_strng =
     &      '(but undetectable with 1% systematic uncertainties in bgd)'
          end if
        else
          xw_strng = 'Count rate too low: 5-sigma detection unlikely'
        end if
        call PWRITE( xw_strng )
        call PWRITE( ' ' )
        call PWRITE( 'Results in the 6 canonical XTE PCA bands are:' )
        call PWRITE( ' ' )
        call PWRITE(
     &          'Channels  Nominal    Source   BGD   5-sigma    (+1%)' )
        call PWRITE(
     &             '          E (keV)    (cps)    (cps) detection (s)' )
        do j = 1, n_res_fixed
          fr_res = results( j ) * frac_PCU
          fr_bgd = bgd_rate( j ) * frac_PCU
          xw_strng = bounds( j )
          xw_strng( 46: 56 ) = '(         )'
          if( fr_res .gt. 5000.0 ) then
            write( xw_strng( 21: 29 ), '(1p,e9.2)' ) fr_res
          else if( fr_res .gt. 0.01 ) then
            write( xw_strng( 21: 29 ), '(f9.3)' ) fr_res
          else
            write( xw_strng( 21: 29 ), '(1p,e9.2)' ) fr_res
          end if
          write( xw_strng( 30: 35 ), '(f6.2)' ) fr_bgd
          if( fr_res .gt. 1.0e15 ) then
            xw_strng( 36: 44 ) = ' ********'
            xw_strng( 47: 55 ) = '*********'
          else
            temp = fr_res * fr_res
            if( temp .gt. 1e-30 ) then
              work = 25.0 * ( fr_res + fr_bgd )
              fivesigma = work / temp
              if( fivesigma .gt. 5000.0 ) then
                write( xw_strng( 36: 44 ), '(1p,e9.2)' ) fivesigma
              else if( fivesigma .gt. 0.01 ) then
                write( xw_strng( 36: 44 ), '(f9.3)' ) fivesigma
              else
                write( xw_strng( 36: 44 ), '(1p,e9.2)' ) fivesigma
              end if
              syst = frac_SYS * fr_bgd
              if( fr_res .gt. 5.0 * syst ) then
                fivesig_b = work / ( temp - 25.0 * syst * syst )
                if( fivesigma .gt. 5000.0 ) then
                  write( xw_strng( 47: 55 ), '(1p,e9.2)' ) fivesig_b
                else if( fivesigma .gt. 0.01 ) then
                  write( xw_strng( 47: 55 ), '(f9.3)' ) fivesig_b
                else
                  write( xw_strng( 47: 55 ), '(1p,e9.2)' ) fivesig_b
                end if
              else
                xw_strng( 47: 55 ) = '*********'
              end if
            else
              xw_strng( 36: 44 ) = ' ********'
              xw_strng( 47: 55 ) = '*********'
            end if
          end if
          call PWRITE( xw_strng )
        end do

        end
