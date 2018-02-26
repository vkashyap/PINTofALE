*+ESPEC
        real function ESPEC( e )

        implicit none

        real e

*       Returns X-ray spectrum in keV/cm/cm/s/keV

        real SPEC

        ESPEC = E * SPEC( E )

        end


*+SPEC
        real function SPEC( E )

        implicit none

        include 'pimms.inc'

        real e

*       Description:
*         Returns photons/cm/cm/s/keV X-ray spectrum at energy e,
*         for several internal models or an external model red from file
*
*       Arguments:
*         e      (i) : energy (keV) of interst
*         <SPEC> (r) : photon spectrum
*
*       Dependencies:
*         None
*
*       Origin:
*         The three standard models from MSSL
*
*       Author and revision history:
*         Koji Mukai (1992) PIMMS version with external model capability
*         Koji Mukai (1993 Mar) first official PIMMS version
*         Koji Mukai (1996 Jun) Protection against overflow in brems
*         Koji Mukai (2000 Apr) Complete re-write to enable multi-component
*                           model.
*-SPEC

        real ez, gal_abs, int_abs, temp
        integer m, m_lo, m_hi

        real TRANMM
        real PMS_BBODY, PMS_BREMS, PMS_POWLW, PMS_EFILE, PMS_GAUSS

        gal_abs = TRANMM( e, nh_g )
	ez = e * z1
*       Morrison & McCannon interstellar absorption

        if( m_crrnt .eq. 0 ) then
          m_lo = 1
          m_hi = n_comp
        else if( m_crrnt .ge. 1 .and. m_crrnt .le. n_comp ) then
          m_lo = m_crrnt
          m_hi = m_crrnt
        else
          call PWRITE( 'ERROR in SPEC:: Model number incorrect' )
          stop
        end if

        temp = 0.0
        do m = m_lo, m_hi
          int_abs = TRANMM( ez, nh( m ) )
          if( model_( m ) .eq. 1 ) then
*           Blackbody
            temp = temp + norm( m ) * int_abs
     &                                  * PMS_BBODY( ez, param( 1, m ) )

          else if( model_( m ) .eq. 2 ) then
*           Bremsstrahlung
            temp = temp + norm( m ) * int_abs
     &                                  * PMS_BREMS( ez, param( 1, m ) )

          else if( model_( m ) .eq. 3 ) then
*           Power Law with optional cut-off
            temp = temp + norm( m ) * int_abs
     &    * PMS_POWLW( ez, param( 1, m ), param( 2, m ), param( 3, m ) )

          else if( model_( m ) .eq. 5 ) then
*           Gaussian
            temp = temp + norm( m ) * int_abs
     &                   * PMS_GAUSS( ez, param( 1, m ), param( 2, m ) )

          else
*           External model, read in from a file.
            temp = temp + norm( m ) * int_abs * PMS_EFILE( ez, m )
          end if

        end do

        SPEC = temp * gal_abs
*       1s in rest frame is stretched to z1 s in observer frame
*       1keV in rest frame is scrunched to 1/z1 keV in observer frame
*       Those cancel out, while photons are conserved and /cm/cm are the same
*       This affects model normalization.

        end



        real function PMS_BBODY( e, kT )

*       implicit none

        real e, kT

        if( e .le. kT * 64.0 ) then
          PMS_BBODY = e * e / ( exp( e / kT ) - 1.0 )
        else
          PMS_BBODY = 0.0
        end if

        end



        real function PMS_BREMS( e, kT )

*       implicit none

        real const
        parameter( const = 0.5513288954 )
*                  const = sqrt( 3.0 ) / 3.14159265359

        real e, kT

        real t, t1_5, nu_9, g

        if( e .gt. kT * 70.0 ) then
          PMS_BREMS = 0.0
        else
          t = kT / 8.61735e-08
*         Temperature (keV) / k
          t1_5 = t * sqrt( t )
          nu_9 = e / 4.135701e-09
*         e / h * 1e-9
          if( t .le. 3.0e+05 ) then
            if( nu_9 .le. t1_5 ) then
              g = const * ( 17.7 + log( t1_5 / nu_9 + 1.0e-09 ) )
            else
              g = 1.0
            end if
          else
            if( e .le. kT ) then
              if( e .ge. kT / 2.79 ) then
                g = 1.0
              else
                g = const * log( 2.2 * kT / e )
              end if
            else
              g = ( kT / e ) ** 0.4
            end if
          end if
          PMS_BREMS = g / e / exp( e / kT )
        end if

        end



        real function PMS_POWLW( e, alpha, highe, ecut )

*       implicit none

        real e, alpha, highe, ecut

        if( highe .eq. 0.0 ) then
*         Power Law
          PMS_POWLW = e ** alpha
        else if( ecut .eq. 0.0 ) then
*         cutoffpl
          PMS_POWLW = e ** alpha * exp( -e / highe )
        else
*         highecut powerlaw
          if( e .le. highe ) then
            PMS_POWLW = e ** alpha
          else
            PMS_POWLW = e ** alpha * exp( ( highe - e ) / ecut )
          end if
        end if

        end



        real function PMS_GAUSS( e, e_cen, sigma )

*       implicit none

        real sqr2pi
        parameter( sqr2pi = 2.506628275 )
*                  Sqrt( 2 * pi )
*       Returns a Gaussian --- normalized such that the integrated area
*                 of the function is 1.0

        real e, e_cen, sigma

        real temp

        if( sigma .gt. 0.0 ) then
*         Actual Gaussian
          temp = ( e - e_cen ) / sigma
          temp = temp * temp
          if( temp .lt. 64.0 ) then
            PMS_GAUSS = exp( -0.5 * temp ) / ( sqr2pi * sigma )
          else
            PMS_GAUSS = 0.0
          end if
        else
*         A delta function
          if( e .eq. e_cen ) then
            PMS_GAUSS = 1.0
          else
            PMS_GAUSS = 0.0
          end if
        end if

        end


        real function PMS_EFILE( e, m )

*       implicit none
*
*       Bug fixed on 2001 May 03 (.lt. rather than .le. in do while)
*           for PIMMS v3.1c

        include 'pimms.inc'

        real e
        integer m

        real e_lo, e_hi, frac
        integer j

*       Do it from a file

        if( e .lt. e_in( beg_m( m ) )
     &                             .or. e .gt. e_in( fin_m( m ) ) ) then
          outob( m ) = .true.
        end if
        if( e .gt. e_old( m ) ) then
          j = j_old( m ) + 1
          e_lo = e_old( m )
          e_hi = e_in( j )
        else
          j_old( m ) = beg_m( m ) - 1
          e_lo = 0.0
          j = beg_m( m )
          e_hi = e_in( j )
        end if
        do while( e .gt. e_hi .and. j .lt. fin_m( m ) )
          j_old( m ) = j
          e_lo = e_hi
          j = j + 1
          e_hi = e_in( j )
        end do
        frac = ( e - e_lo ) / ( e_hi - e_lo )
        e_old( m ) = e_lo
        PMS_EFILE = f_in( j_old( m ) ) * ( 1.0 - frac )
     &                                                + f_in( j ) * frac

        end
