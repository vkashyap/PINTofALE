*+PMS_RAREA
        subroutine PMS_RAREA( mission, detector, filter, m_bin, m_band,
     &                          e, a, n_bin, n_band, in_lo, in_hi )

        implicit none

        integer mission, detector, filter
        integer m_bin, m_band, n_bin, n_band
        real e( m_bin ), a( m_bin, 0: m_band )
        real in_lo( 0: m_band ), in_hi( 0: m_band )

*       Description:
*         Reads in the effective area file for the given combination of
*         instruments into the given array
*
*       Arguments:
*         mission      (i) : input mission number
*         detector     (i) : input detector number
*         filter       (i) : input filter number
*         m_bin        (i) : declared size of the arrays below
*         m_band       (i) : declared size of the arrays below
*         e, a         (o) : arrays containing energy and area
*         n_bin        (o) : actual size, from the file
*         n_band       (o) : actual number of PHA sets.
*         in_lo, in_hi (o) : energy range, from the file
*
*       Dependencies:
*         XANLIB routines, PIMMS Index routines
*
*       Origin:
*         Created by KM for PIMMS
*
*       Author:
*         Koji Mukai, 1993 March, first official PIMMS version
*         Modified by KM, 1994 Dec, for the XTE small PIMMS
*         Modified by KM, 2002 Jan:
*             Fixed the bug so that the limit is m_bin, not m_bin - 1
*-PMS_RAREA

        real drivel
        integer lun, j, cols, k, flag
        character*256 af_name, line

        include 'sitespec.inc'

        call PDX_MKNAM
     &              ( mission, detector, filter, 'AREA', af_name, flag )
        if( flag .lt. 0 ) then
          call PWRITE(
     &          'FATAL ERROR in PMS_RAREA:: file name cannot be found' )
          stop
        end if
        call ARKOPN( lun, ddir_name, af_name, 'area', 'OLD', 'READONLY',
     &                  'FORMATTED', 'SEQUENTIAL', 1, flag )
        if( flag .ne. 0 ) then
          call PWRITE(
     &         'FATAL ERROR in PMS_RAREA:: file cannot be opened' )
          stop
        end if
        call PMS_COLMN( lun, line, cols )
        if( cols .lt. 2 ) then
          call PWRITE(
     &        'FATAL ERROR in PMS_RAREA:: file in unrecognized format' )
          stop
        end if
        n_band = cols - 2
        read( line, * ) e( 1 ), ( a( 1, k ), k = 0, n_band )
        j = 2
        do while( j .le. m_bin )
          read( lun, *, end = 100 ) e( j ), ( a( j, k ), k = 0, n_band )
          j = j + 1
        end do
        read( lun, *, end = 100 ) drivel
*       if program control comes here, we've overfilled the array!
        call PWRITE( 'FATAL ERROR in PMS_RAREA:: too many points' )
        stop

100     continue
        close( lun )
        n_bin = j - 1
        do k = 0, n_band
          in_lo( k ) = e( n_bin )
          in_hi( k ) = e( 1 )
          do j = 1, n_bin
            if( a( j, k ) .gt. 0.0 ) then
              in_lo( k ) = min( in_lo( k ), e( j ) )
              in_hi( k ) = max( in_hi( k ), e( j ) )
            end if
          end do
        end do

        end
