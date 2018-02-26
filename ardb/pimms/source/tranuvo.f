        real function TRANUVO( E, Nh )

        implicit none

        real E, Nh

        double precision const, kev2mic, magconst
        parameter( const = 8.5416667d-22 )
c       Eb-v = Nh / K, where K = 4.8d+21
c                      from Line 6 of Table 2 of Bohlin et al. (1978)
c       A_b = ( 1.0 + R_v ) * Eb-v, where R_v = 3.1
        parameter( kev2mic = 1.239854d-03 )
c       keV to micron conversion.
        parameter( magconst = -0.921034037d+00 )
c       -ln( 10.0 ) / 2.5, so exp( magconst * mag ) gives the correct anser.

C       UV/Optical/IR extinction routine, based on input from
C       Stephen Holland.
C       This version by Koji Mukai, 2003 September

        double precision lambda, a_b, extinct2

        double precision EXTINCT

        if( E .le. 0.0 ) then
          TRANUVO = 1.0
          return
        end if
        lambda = kev2mic / E
        a_b = Nh * const
        extinct2 = EXTINCT( lambda, 1, a_b )
        if( extinct2 .le. 0.0d+00 ) then
          TRANUVO = 1.0
        else
          TRANUVO = exp( magconst * extinct2 )
        end if

        end
c
c
c
        block data SET_EXTIN_COEFFS

        integer maxterms, maxlaws
        parameter( maxterms = 6 )
        parameter( maxlaws = 3 )

        double precision a_i( maxterms, maxlaws )
        double precision lambda_i( maxterms, maxlaws )
        double precision b_i( maxterms, maxlaws )
        double precision n_i( maxterms, maxlaws )

        common / EXTIN_COEFFS / a_i, lambda_i, b_i, n_i

        data a_i / 165.0d+00, 14.0d+00, 0.045d+00,
     &                      0.002d+00, 0.002d+00, 0.012d+00,
     &             175.0d+00, 19.0d+00, 0.023d+00,
     &                      0.005d+00, 0.006d+00, 0.020d+00,
     &             185.0d+00, 27.0d+00, 0.005d+00,
     &                      0.010d+00, 0.012d+00, 0.030d+00 /

        data lambda_i / 0.047d+00, 0.08d+00, 0.22d+00,
     &                      9.7d+00, 18.0d+00, 25.0d+00,
     &                  0.046d+00, 0.08d+00, 0.22d+00,
     &                      9.7d+00, 18.0d+00, 25.0d+00,
     &                  0.042d+00, 0.08d+00, 0.22d+00,
     &                      9.7d+00, 18.0d+00, 25.0d+00 /

        data b_i / 90.0d+00, 4.0d+00, -1.95d+00,
     &                      -1.95d+00, -1.80d+00, 0.0d+00,
     &             90.0d+00, 5.5d+00, -1.95d+00,
     &                      -1.95d+00, -1.80d+00, 0.0d+00,
     &             90.0d+00, 5.5d+00, -1.95d+00,
     &                      -1.95d+00, -1.80d+00, 0.0d+00 /

        data n_i / 2.0d+00, 6.5d+00, 2.0d+00,
     &                      2.0d+00, 2.0d+00, 2.0d+00,
     &             2.0d+00, 4.5d+00, 2.0d+00,
     &                      2.0d+00, 2.0d+00, 2.0d+00,
     &             2.0d+00, 4.5d+00, 2.0d+00,
     &                      2.0d+00, 2.0d+00, 2.0d+00 /

        end

C LAW=0 is no extinction, LAW=1 is Milkyway, LAW=2 is LMC, and LAW=3 is SMC
C TERM=1 is BKG, TERM=2 is FUV, TERM=3 is 2175AA, TERM=4 is 9.7 um,
C                TERM=5 is 18 um, and TERM=6 is FIR
C Data from Pei (1992) Table 4.

c*****
      DOUBLE PRECISION FUNCTION EXTINCT(LAMBDA,LAW,A_B)
c
c/===================================================================
c/ Computes the extinction based on Pei, Y. C., 1992, ApJ, 395, 130
c/===================================================================
c
      IMPLICIT NONE


        integer maxterms, maxlaws
        parameter( maxterms = 6 )
        parameter( maxlaws = 3 )

        double precision a_i( maxterms, maxlaws )
        double precision lambda_i( maxterms, maxlaws )
        double precision b_i( maxterms, maxlaws )
        double precision n_i( maxterms, maxlaws )

        common / EXTIN_COEFFS / a_i, lambda_i, b_i, n_i
c
        DOUBLE PRECISION LAMBDA
c                                     ! Wavelength in microns
        INTEGER LAW
c                                     ! Extinction law to use
        DOUBLE PRECISION A_B
c                                     ! B-band extinction
c
        DOUBLE PRECISION XI
c                                     ! Defined to be A_lambda / A_B
        DOUBLE PRECISION TERM(MAXTERMS)
c
        INTEGER I, J
c

c--------------------------------------------------------------------
c

        IF ( LAW .EQ. 0 ) THEN
           EXTINCT = 0.0d0
        ELSE IF ( ( LAW .GE. 1 ) .AND. ( LAW .LE. MAXLAWS ) ) THEN
           XI = 0.0d0
           DO I = 1, MAXTERMS, 1
              TERM(I) = A_I(I,LAW)
     &           / ( ( LAMBDA / LAMBDA_I(I,LAW) )**N_I(I,LAW)
     &           + ( LAMBDA_I(I,LAW) / LAMBDA)**N_I(I,LAW)
     &           + B_I(I,LAW) )
              XI = XI + TERM(I)
           END DO
           EXTINCT = XI * A_B
        ELSE
c         WRITE( UNIT=*, FMT=1000 ) LAW
c 1000    FORMAT(' extinct: invalid extinction law (LAW=',1I1,')')
           EXTINCT = -1.0d0
c         IERR = 1
        END IF
c
c      WRITE(*,*) 'extinct: xi = ',xi
c
        RETURN
        END
