*+PMS_INTLZ
        subroutine PMS_INTLZ( )

        implicit none

        include 'pms_index.inc'

*       Description:
*         Initializes the mission list common block for PIMMS
*
*       Arguments:
*         None
*
*       Dependencies:
*         Index subroutines in PMS_INDEX.FOR; LENTRIM
*
*       Origin:
*         Cooked up by KM to handle the PIMMS-specific database.
*
*       Author:
*         Koji Mukai, 1993 Feb 23, Original version
*
*         Koji Mukai, 2010 Jan 15, bug fix for mission with no choice
*                     of detectors or filters
*-PMS_INTLZ

        integer msn_unit
        integer flag

        integer n_detect, n_filter, last_index, next_len
        integer i, j, j1, k, k1, status
        character*( max_len ) next_sat, next_det, next_fil

        integer LENTRIM, PDX_GETDT, PDX_GETFL

        include 'sitespec.inc'

        call ARKOPN( msn_unit, ddir_name, 'pms_mssn.lst', 'lst',
     &               'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &               1, flag )
        if( flag .ne. 0 )  then
          call PWRITE(
     &             'FATAL ERROR:: failed to open mission list file' )
          call PWRITE( 'PIMMS was looking for ''pms_mssn.lst'' in' )
          call PWRITE( ddir_name( : LENTRIM( ddir_name ) ) )
          stop
        end if

        n_mssn = 0
        last_index = 0
        do while( .true. )
          read( msn_unit, *, end = 800 ) next_sat, n_detect
          call LOCASE( next_sat )
          next_len = LENTRIM( next_sat )
          if( next_len .gt. max_len ) then
            flag = -2
            goto 900
          end if
          do i = 1, n_mssn
            if( next_sat .eq. nm_array( ind_mssn( i ) ) ) then
              flag = -11
              goto 900
            end if
          end do
          n_mssn = n_mssn + 1
          if( n_mssn .gt. max_mssn ) then
            flag = -12
            goto 900
          end if
          last_index = last_index + 1
          if( last_index .gt. max_index ) then
            flag = -3
            goto 900
          end if
          ind_mssn( n_mssn ) = last_index
          nm_array( last_index ) = next_sat
          ln_array( last_index ) = next_len
          if( n_detect .eq. 0 ) then
*           No subdivision into detectors etc. here
            st_array( last_index ) = -1
            call PDX_CKALL( i, 0, 0, status )
            st_array( last_index ) = -status
          else
            st_array( last_index ) = n_detect
            do j = 1, n_detect
              read( msn_unit, * ) next_det, n_filter
              call LOCASE( next_det )
              next_len = LENTRIM( next_det )
              if( next_len .gt. max_len ) then
                flag = -2
                goto 900
              end if
              do j1 = 1, j - 1
                if( next_det .eq. nm_array( PDX_GETDT( i, j1 ) ) ) then
                  flag = -11
                  goto 900
                end if
              end do
              last_index = last_index + 1
              if( last_index .gt. max_index ) then
                flag = -3
                goto 900
              end if
              call PDX_SETDT( i, j, last_index )
              nm_array( last_index ) = next_det
              ln_array( last_index ) = next_len
              if( n_filter .eq. 0 ) then
*               No subdivision into filters here
                call PDX_CKALL( i, j, 0, status )
                st_array( last_index ) = -status
              else
                st_array( last_index ) = n_filter
                do k = 1, n_filter
                  read( msn_unit, * ) next_fil
                  call LOCASE( next_fil )
                  next_len = LENTRIM( next_fil )
                  if( next_len .gt. max_len ) then
                    flag = -2
                    goto 900
                  end if
                  do k1 = 1, k - 1
                    if( next_fil
     &                    .eq. nm_array( PDX_GETFL( i, j, k1 ) ) ) then
                      flag = -11
                      goto 900
                    end if
                  end do
                  last_index = last_index + 1
                  if( last_index .gt. max_index ) then
                    flag = -3
                    goto 900
                  end if
                  call PDX_SETFL( i, j, k, last_index )
                  nm_array( last_index ) = next_fil
                  ln_array( last_index ) = next_len
                  call PDX_CKALL( i, j, k, status )
                  st_array( last_index ) = status
                end do
              end if
            end do
          end if
        end do

800     continue
*       Normal end
        close( msn_unit )
        return

900     continue
*       Error condition occured

        end
