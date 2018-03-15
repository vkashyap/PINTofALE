*       Updated 2014 January with inputs from Antonella Fruscione

*       Note --- Identical parameters appear in several subroutine
*                in this file.  If you change one, change all.

        subroutine ACIS_CHOICE( cps_given )

        implicit none

        real frame_def, frame_fast, accept_pu
        parameter( frame_def = 3.2 )
*                  Default frame time is 3.2 s
        parameter( frame_fast = 0.2 )
*                  Fastest single CCD readout is 0.2 s
        parameter( accept_pu = 0.1 )
*                  Proposers' Guide Sec 5.9 contains the guideline
*                  "pileup fraction >10% means greatly affected"

        real cps_given

*       This subroutine is used for ACIS non-grating observations
*       Check for the second numerical parameter on the command line
*       and if present, try to interpret as the optional frame time
*       argument, and estimate the pile-up fraciton accordingly.
*       If not, suggest a sensible frame time to minimize pile-up.

        real frame_time, cps_piled, pileup_frac, rev_ret
        integer num_c, num_f, num_i, pu_flag
        character*80 p_strng

        p_strng = ' '
        call INQ_PARAM( num_c, num_f, num_i )
        if( num_f + num_i .eq. 2 ) then
*         The user has typed an extra number, try to interpret that as
*            the (non-standard) frame time and calculate the pileup
          call GET_PARAM_N( 2, frame_time )
          if( frame_time .gt. frame_def )then
            call PWRITE( '  The input frame time is too long.' )
          else if( frame_time .lt. frame_fast ) then
            call PWRITE( '  The input frame time is too short.' )
          else
            call DO_PILEUP
     &                 ( cps_given, frame_time, cps_piled, pileup_frac )
            p_strng( : 31 ) = '  For a frame time of       s, '
            write( p_strng( 23: 27 ), '(f5.3)' ) frame_time
            p_strng( 32: 63 ) = 'the count rate after pile-up is '
            if( cps_piled .ge. 100.0 ) then
              write( p_strng( 64: 73 ), '(1p,e10.3)' ) cps_piled
            else if( cps_piled .ge. 0.001 ) then
              write( p_strng( 64: 73 ), '(f10.5)' ) cps_piled
            else
              write( p_strng( 64: 73 ), '(1p,e10.3)' ) cps_piled
            end if
            call PWRITE( p_strng )
            p_strng = '  (or a pile-up fraction of        %)'
            write( p_strng( 28: 34 ), '(f7.3)' ) pileup_frac * 100.0
            call PWRITE( p_strng )
          end if

        else
*         Otherwise, run PILEUP for the default frame time, and if
*            pileup looks important, try to guess what frame time
*            might be necessary to reduce it to a manageable level.
          call PILEUP_REV( cps_given, pu_flag, rev_ret )
          if( pu_flag .eq. -1 ) then
*           No (significant) pile-up at normal frame time
            p_strng( : 30 ) = '  Pile-up is not significant ('
            write( p_strng( 31: 36 ), '(f6.3)' ) rev_ret * 100.0
            p_strng( 37: 67 ) = ' %) at normal frame-time (3.2s)'
            write( p_strng( 63: 65), '(f3.1)' ) frame_def
            call PWRITE( p_strng )
          else if( pu_flag .eq. 0 ) then
*           Frame time given such that pile-up is at threshold value.
            p_strng( : 41 ) =
     &                       '  Pile-up is generally tolerable (     %)'
            write( p_strng( 35: 38 ), '(f4.1)' ) accept_pu * 100.0
            p_strng( 42: 69 ) = ' at a frame-time of        s'
            write( p_strng( 62: 67 ), '(f6.3)' ) rev_ret
            call PWRITE( p_strng )
          else
*           Pile-up severe even at the fastest frame time
            p_strng( : 38 ) = '  Pile-up is too high (     %) at the '
            write( p_strng( 24: 27 ), '(f4.1)' ) rev_ret * 100.0
            p_strng( 39: 76 ) = 'fastest single-chip frame time (0.2 s)'
            write( p_strng( 71: 73 ), '(f3.1)' ) frame_fast
            call PWRITE( p_strng )
            p_strng = '  Consult the Chandra POG for mitigation methods'
            call PWRITE( p_strng )
          end if

        end if

        end


        subroutine PILEUP_REV( cps_given, pu_flag, rev_ret )

        implicit none

*       Given the naively calculated count rate, returns a status flag
*       pu_flag and rev_ret.

        real frame_def, frame_fast, accept_pu
        integer max_iter
        parameter( frame_def = 3.2 )
*                  Default frame time is 3.2 s
        parameter( frame_fast = 0.2 )
*                  Fastest single CCD readout is 0.2 s
        parameter( accept_pu = 0.1 )
*                  Proposers' Guide Sec 5.9 contains the guideline
*                  "pileup fraction >10% means greatly affected"
        parameter( max_iter = 10 )
*                  How many iterations in bisect method?

*
*       Input
        real cps_given
*                      : Naively calculated count rate
*       Output
        integer pu_flag
*                       : -1 if no pileup at normal (frame_def) frame time
*                          0 if in between
*                          1 if significant pileup at fastest frame time
*                               without going to continuous readout (frame_fast)
*             where significant pileup is defined by the
*             parameter 'accept_pu' currently set to 10%
        real rev_ret
*          : if pu_flag is non-0, the pileup fraction at def or fast frame time
*            if pu_flag is 0, frame time such that pileup fraction = accept_pu

        real cps_piled, pileup_frac
        real frame_lo, frame_hi, frame_test
        integer j

        call DO_PILEUP( cps_given, frame_def, cps_piled, pileup_frac )
        if( pileup_frac .le. accept_pu ) then
          pu_flag = -1
          rev_ret = pileup_frac
        else
          call DO_PILEUP
     &                 ( cps_given, frame_fast, cps_piled, pileup_frac )
          if( pileup_frac .gt. accept_pu ) then
            pu_flag = +1
            rev_ret = pileup_frac
          else
            pu_flag = 0
            frame_lo = frame_fast
            frame_hi = frame_def
            do j = 1, max_iter
              frame_test = ( frame_lo + frame_hi ) * 0.5
              call DO_PILEUP
     &                 ( cps_given, frame_test, cps_piled, pileup_frac )
              if( pileup_frac .eq. accept_pu ) then
                goto 100
              else if( pileup_frac .lt. accept_pu ) then
                frame_lo = frame_test
              else
                frame_hi = frame_test
              end if
            end do
            frame_test = ( frame_lo + frame_hi ) * 0.5
 100        continue
            rev_ret = frame_test
          end if
        end if

        end



        subroutine ACIS_NOCHOICE( cps_given )

        implicit none

        real frame_def
        parameter( frame_def = 3.2 )
*                  Default frame time is 3.2 s

        real cps_given
*                      : Input count rate
*       This routine is called for 0th order calculations
*         --- just run PILEUP, print out the results.

        real cps_piled, pileup_frac
        character*80 p_strng

        call DO_PILEUP( cps_given, frame_def, cps_piled, pileup_frac )

        p_strng = ' '
        p_strng( : 34 ) = '  The piled-up count rate is '
        if( cps_piled .ge. 100.0 ) then
          write( p_strng( 35: 44 ), '(1p,e10.3)' ) cps_piled
        else if( cps_piled .ge. 0.001 ) then
          write( p_strng( 35: 44 ), '(f10.5)' ) cps_piled
        else
          write( p_strng( 35: 44 ), '(1p,e10.3)' ) cps_piled
        end if
        p_strng( 45: 71 ) = ' (3.2 s frame time assumed)'
        write( p_strng( 47: 49 ), '(f3.1)' ) frame_def
        call PWRITE( p_strng )
        p_strng = '  (or a pile-up fraction of        %)'
        write( p_strng( 28: 34 ), '(f7.3)' ) pileup_frac * 100.0
        call PWRITE( p_strng )

        end


	subroutine DO_PILEUP
     &                 ( cps_given, frame_time, cps_piled, pileup_frac )

        implicit none

	real area_ratio, eef_exp
	parameter( area_ratio = 0.125 )
*       ratio of area of inner detect cell (1.0) to that of outer cells (8.0)
        parameter( eef_exp = 0.8860 )
*       Expected encircled energy fraction in inner detect cell

        real cps_given
*                          Naively calculated count rate (counts/sec)
        real frame_time
*                          Frame time in second
        real cps_piled
*                          Actual detected count rate with pile-up
        real pileup_frac
*                          Pileup fraction

* Approximate pileup test values (accurate to 10-20%):
*
*     input rate  ct_rate   detected rate    pileup fraction
*       0.03        0.10      0.095            0.038  (3.8%)
*       0.09        0.30      0.25             0.12   (12%)
*       0.24        0.80      0.45             0.28   (28%)
*       0.61        2.0       0.52             0.57   (57%)
*
*      Note:  model starts to break above 0.5 for an input rate.

        real ot_e_weight
*                          Weighting factor for outer cell
        real cpf_given
*                          Counts per Frame
        real cpf_exp_in, cpf_exp_ot
*             Expected count rates in inner detect cell and in outer cells
        real cpf_piled
*                        Total single photon detection rates
        real cpf_piled_in, cpf_piled_ot
*                        Single photon detection rates in inner & outer cells
        real zrate_pin, zrate_pot
*                        Zero photon detection rates in inner & outer cells
	real pfrac_in, pfrac_ot
*                        Pile-up fractions in inner & outer cells


* Compute weighting factor for encircled energy in outer detect cell
        ot_e_weight = ( 1.0 - eef_exp ) * area_ratio

* Desired quantity is counts/frame.  Assume all photons land
* in one pixel (this is the problem of pileup, after all).  Then
* the input count rate must be multiplied by the frame time, so
*
*	counts/frame = counts/sec * sec/frame
        cpf_given = cps_given * frame_time

* Calculate expected count rates in each detect cell.
        cpf_exp_in = cpf_given * eef_exp
        cpf_exp_ot = cpf_given * ot_e_weight

* Calculate the single photon detection rate in each detect cell.
        cpf_piled_in = cpf_exp_in * exp( -cpf_exp_in )
        cpf_piled_ot = cpf_exp_ot * exp( -cpf_exp_ot )

* Calculate the detected (piled) count rate.
        cpf_piled = cpf_piled_in + cpf_piled_ot / area_ratio
        cps_piled = cpf_piled / frame_time

* Calculate pileup fraction.

* Calculate the zero photon detection rate in each detect cell.
        zrate_pin = exp( -cpf_exp_in )
        zrate_pot = exp( -cpf_exp_ot )

* Calculate pileup fraction in each detect cell.  Weight by the 
* encircled energy fractions.

        pfrac_in = ( 1.0 - ( zrate_pin + cpf_piled_in ) )
     &                                        / ( 1.0 - zrate_pin )
        pfrac_in = pfrac_in * eef_exp
        pfrac_ot = ( 1.0 - ( zrate_pot + cpf_piled_ot ) )
     &                                        / ( 1.0 - zrate_pot )
        pfrac_ot = pfrac_ot * ( 1.0 - eef_exp )

        pileup_frac = pfrac_in + pfrac_ot

        end
