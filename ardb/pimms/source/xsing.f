        character*128 ifname, ofname, drivel
        integer lun1, lun7, last, LENTRIM

        call WRITEN( 'Input file name >' )
        read '(a)', ifname
        call WRITEN( 'Output file name >' )
        read '(a)', ofname
        call ARKOPN( lun1, ' ', ifname, 'dat',
     &               'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &               1, last )
        read( lun1, '(a)' ) drivel
        read( lun1, '(a)' ) drivel
        read( lun1, '(a)' ) drivel
        call ARKOPN( lun7, ' ', ofname, 'mdl',
     &               'NEW', 'OVERWRITE', 'FORMATTED', 'SEQUENTIAL',
     &               1, lsat )
        do while( .true. )
          read( lun1, '(a)', end = 100 ) drivel
          last = LENTRIM( drivel )
          read( drivel, * ) x, dx, y
          do while( drivel( last: last ) .eq. '-' )
            read( lun1, '(a)', end = 100 ) drivel
            last = LENTRIM( drivel )
          end do
          write( lun7, * ) x, y
        end do
100     continue
        close( lun1 )
        close( lun7 )
        end
