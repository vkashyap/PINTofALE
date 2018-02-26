*+SIS_OBSTM
        subroutine SIS_OBSTM
     &  ( flag, rate, bg_rate, exp_time, sn_ratio, r_opt, s_opt, b_opt )

        implicit none

        integer m_psf
        parameter( m_psf = 512 )

        integer flag
        real rate, bg_rate
        real exp_time, sn_ratio
        real r_opt, s_opt, b_opt

*       Description:
*         Given a point source count rate and a background rate, SIS_OBSTM
*         calculates the S/N for a given exposure time or vice versa.
*
*       Arguments:
*         flag         (i)   : 1 for exposure time to S/N, 2 for reverse
*         rate         (i)   : SIS count rate of the point source 
*         bg_rate      (i)   : Guesstimated background count rate
*                              per pixel per second.  Use 2.12E-07 (nominal)
*                              unless the user wants something different!
*         exp_time   (i or o): Exposure time in question
*         sn_ratio   (o or i): S/N ratio in question
*         r_opt        (o)   : Radius of the optimal extraction cell
*         s_opt        (o)   : Source rate in the above
*         b_opt        (o)   : Background rate in the above
*
*       Dependences:
*         The PSF file must be read in by this subroutine when it is run for
*         the first time.  In this version, the file I/O is badly hard-wired.
*         In *THIS* version, ARKio is used
*
*       Origin:
*         The source detection uses a circular aperture, background is assumed
*         to be well determined (maybe from the data itself, outside the source
*         aperture, maybe from some model).  Using the (assumed to be uniform)
*         background level and the source count rate, the radius of the
*         aperture is optimized, then the S/N from a 1 sec exposure is
*         calculated.  The desired results follow from there.
*
*       Author:
*         Koji Mukai,  1992 Dec 17, Original version
*         Koji Mukai,  1995 Aug, added extra outputs
*
*-SIS_OBSTIM

        real r( m_psf ), eef( m_psf )
        integer n_psf
*                        Above are the psf (eef) data and its actual size
        real temp, b_now, sn1, max_sn1
        integer j, lun, ierr

        logical first

        include 'sitespec.inc'

        data first / .true. /

        if( first ) then
          call ARKOPN( lun, ddir_name, 'asca_sis_psf.dat', 'dat',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, ierr )
          if( ierr .lt. 0 ) then
            call PWRITE( 'SEVERE ERROR:: Datafile not found' )
            stop
          end if
          j = 1
          do while( .true. )
            read( lun, *, end = 100 ) r( j ), eef( j )
            j = j + 1
          end do
100       continue
          close( lun )
          n_psf = j - 1
          first = .false.
        end if

        r_opt = r( 1 )
        s_opt = rate * eef( 1 )
        b_opt = bg_rate * 6.2831853 * r( 1 ) * r( 1 )
        max_sn1 = s_opt / sqrt( s_opt + b_opt )
*       SN1 = Encircled source counts / Sq-root (Enc. src. counts + bgd)
*            bgd = background rate per pixel times area
*       This is the S/N for a 1 sec observation
        do j = 2, n_psf
          temp = rate * eef( j )
          b_now = bg_rate * 6.2831853 * r( j ) * r( j )
          sn1 = temp / sqrt( temp + b_now )
          if( sn1 .gt. max_sn1 ) then
            r_opt = r( j )
            s_opt = temp
            b_opt = b_now
            max_sn1 = sn1
          end if
        end do

        if( flag .eq. 1 ) then
*         do exposure time to S/N
          sn_ratio = max_sn1 * sqrt( exp_time )
        else if( flag .eq. 2 ) then
*         do S/N to exposure time
          temp = sn_ratio / max_sn1
          exp_time = temp * temp
        end if

        end
