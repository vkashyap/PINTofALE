*+PMS_OUTPT
        subroutine PMS_OUTPT( flag )

        implicit none

        integer flag

        include 'pimms.inc'

*       Description:
*         Writes out the photon spectrum to a file in the specified
*         energy range and a specified step.  For multi-component
*         models, the format is E: total: comp1: comp2: comp3 ...
*
*       Author:
*         Koji Mukai 2000 April, original version
*-PMS_OUTPT

        real cpnts( m_mdls )
        real lo_e, hi_e, delta, e_now
        integer lun_out, j, n_bins, m
        character*64 of_name

        real SPEC

        call PARAM_C( 1, of_name, 'Enter output file name >' )
        call LOCASE( of_name )
        flag = 0
        call ARKOPN( lun_out, ' ', of_name, 'mdl', 'NEW', 'OVERWRITE',
     &                            'FORMATTED', 'SEQUENTIAL', 1, flag )
        if( flag .lt. 0 ) then
          call PWRITE( 'Error opening output file ' )
          flag = -349
          return
        end if
        call PARAM_N( 1, lo_e, 2.0,
     &                       'Enter lowest energy of interest (keV) >' )
        call PARAM_N( 2, hi_e, 10.0,
     &                      'Enter highest energy of interest (keV) >' )
        if( hi_e .le. lo_e .or. lo_e .le. 0.0 ) then
          call PWRITE( 'ERROR in energy range' )
          flag = -350
          return
        end if
        call PARAM_N( 3, delta, 0.002, 'Enter energy step (keV) >' )
        if( delta .le. 0.0 ) then
          call PWRITE( 'ERROR in energy step' )
          flag = -351
          return
        end if
        n_bins = ( hi_e - lo_e ) / delta
        if( n_comp .ge. 2 ) then
          do j = 0, n_bins
            e_now = lo_e + j * delta
            do m = 1, n_comp
              m_crrnt = m
              cpnts( m ) = SPEC( e_now )
            end do
            m_crrnt = 0
            write( lun_out, 100 ) e_now, SPEC( e_now ),
     &                                     ( cpnts( m ), m = 1, n_comp )
 100        format( f9.4, ' ', 1p, 9( e10.3, ' ' ) )
          end do
        else
          do j = 0, n_bins
            e_now = lo_e + j * delta
            write( lun_out, * ) e_now, SPEC( e_now )
          end do
        end if
        close( lun_out )

        end
