*+PMS_DOCRT

        subroutine PMS_DOCRT( in_rate, f_mssn, f_dtct, f_fltr, f_unit,
     &          f_lo, f_hi, t_mssn, t_dtct, t_fltr, t_unit, t_lo, t_hi )

        implicit none

        integer m_bin, m_band
        parameter( m_bin = 32768 )
        parameter( m_band = 10 )

        real in_rate
        integer f_mssn, f_dtct, f_fltr, f_unit
        real f_lo, f_hi
        integer t_mssn, t_dtct, t_fltr, t_unit
        real t_lo, t_hi

*       Description:
*         Integrates the model over the range, perhaps convoluted with
*         the instrument effective area curve, twice --- once for
*         an instrument (or none) for which the count rate (or flux)
*         is known, the second time for the target instrument.
*         This therefore does the count rate conversion in one step.
*
*       Arguments:
*         in_rate          (i) : known count rate/flux
*         f_mssn           (i) :   mission number for above (0 for flux)
*         f_dtct           (i) :   detector number for above
*         f_fltr           (i) :   filter number for above
*         f_unit           (i) :   unit code for above (if flux)
*         f_lo, f_hi       (i) :   integration range for above
*         t_mssn           (i) : * mission number for target (0 for flux)
*         t_dtct           (i) :   detector number for target
*         t_fltr           (i) :   filter number for target
*         t_unit           (i) :   unit code for target (if flux)
*         t_lo, t_hi       (i) :   integration range for target
*
*       Dependencies:
*         Other PIMMS routines
*
*       Origin:
*         Created for PIMMS by KM
*
*       Author:
*         Koji Mukai   (1992)  Test version for non-XPI PIMMS
*         Koji Mukai (1993 Mar) First official version
*         Koji Mukai (1995 Jul) Introduced restr in call to PMS_SPECL
*         Koji Mukai, 1997 Dec - renamed abs to column to avoid
*                                conflict with the Fortran intrinsic abs
*                              - Allowed energy range expressed in negative
*                                numbers, signifying that they were originally
*                                input as wavelength range in Angstroms
*         Koji Mukai, 2002 Nov - now accepts flux density for PIMMS v3.3
*         Koji Mukai, 2008 Oct - changed energy boundaries for results(1:n)
*-PMS_DOCRT

        integer sv_mssn( 2 ), sv_dtct( 2 ), sv_fltr( 2 )
        real e( m_bin, 2 ), a( m_bin, 0: m_band, 2 )
        real sv_lo( 0: m_band, 2 ), sv_hi( 0: m_band, 2 )
        integer n_bin( 2 ), n_band( 2 )
        save sv_mssn, sv_dtct, sv_fltr
        save e, a, sv_lo, sv_hi, n_bin, n_band

        real lo, hi
        logical restr
        integer go_bin, b
        real f_integ, t_integ, results( 0: m_band )
        character*80 xw_strng

        real PMS_INTEG

        if( f_mssn .eq. 65536 ) then
*         Input is model normalization
          f_integ = 1.0

        else
          if( f_mssn .gt. 0 .and. f_mssn .lt. 8192 ) then
*           Known rate is for an actual mission, not flux or flux density

            if( f_mssn .ne. sv_mssn( 1 ) .or. f_dtct .ne. sv_dtct( 1 )
     &                              .or. f_fltr .ne. sv_fltr( 1 ) ) then
*             Must read new effective area curve into buffer
              call PMS_RAREA( f_mssn, f_dtct, f_fltr, m_bin, m_band,
     &                      e( 1, 1 ), a( 1, 0, 1 ), n_bin( 1 ),
     &                      n_band( 1 ), sv_lo( 0, 1 ), sv_hi( 0, 1 ) )
              sv_mssn( 1 ) = f_mssn
              sv_dtct( 1 ) = f_dtct
              sv_fltr( 1 ) = f_fltr
            end if

            if( f_lo .eq. 0.0 .and. f_hi .eq. 0.0 ) then
*             User did not specify integration boundaries; use the limits
*             of the effective area curve
              lo = sv_lo( 0, 1 )
              hi = sv_hi( 0, 1 )
            else
*             Use the user-supplied value
*                Note negative values indicate that the user original
*                input the range in A; they have been converted to -keV
              lo = abs( f_lo )
              hi = abs( f_hi )
            end if
            go_bin = n_bin( 1 )

          else if( f_mssn .eq. 131072 ) then
*           Flux Density; there must be user-supplied boundaries (see above)
            lo = abs( f_lo )
            if( f_lo .lt. 0.0 ) then
*             Flux density per Angstrom
              hi = -1
            else
*             Flux density per keV
              hi = +1
            end if
            go_bin = 16

          else
*           Flux; there must be user-supplied boundaries (see above)
            lo = abs( f_lo )
            hi = abs( f_hi )
            go_bin = 16
          end if

          f_integ = PMS_INTEG( f_mssn, f_unit, lo, hi,
     &                                 e( 1, 1 ), a( 1, 0, 1 ), go_bin )
        end if

        call PMS_WTMDL( 2 )
        call PMS_WTINT
     &        ( f_mssn, f_dtct, f_fltr, f_unit, f_lo, f_hi, in_rate, 1 )
        if( f_mssn .gt. 0 .and. f_mssn .lt. 8132 ) then
          call PMS_SPEC2( f_mssn, f_dtct, f_fltr )
        end if

        if( f_integ .eq. 0.0 ) then
          call PWRITE(
     &    'PIMMS predicts 0.0 flux/count rate for original instrument' )
          call PWRITE( 'Conversion is therefore impossible' )
          return
        else
          write( xw_strng,
     &        '(''  (Internal model normalization = '',1pe10.3,'')'')' )
     &                                                 in_rate / f_integ
          call PWRITE( xw_strng( : 46 ) )
        end if

        if( t_mssn .eq. 65536 ) return
*       Model normalization has already been reported

        if( t_mssn .gt. 0 .and. t_mssn .lt. 8192 ) then
*         Target is an actual mission, not flux

          if( t_mssn .ne. sv_mssn( 2 ) .or. t_dtct .ne. sv_dtct( 2 )
     &                              .or. t_fltr .ne. sv_fltr( 2 ) ) then
*           Must read new effective area curve into buffer
            call PMS_RAREA( t_mssn, t_dtct, t_fltr, m_bin, m_band,
     &                      e( 1, 2 ), a( 1, 0, 2 ), n_bin( 2 ),
     &                      n_band( 2 ), sv_lo( 0, 2 ), sv_hi( 0, 2 ) )
            sv_mssn( 2 ) = t_mssn
            sv_dtct( 2 ) = t_dtct
            sv_fltr( 2 ) = t_fltr
          end if

          if( t_lo .eq. 0.0 .and. t_hi .eq. 0.0 ) then
*           User did not specify integration boundaries; use tho limits
*           of the effective area curve
            lo = sv_lo( 0, 2 )
            hi = sv_hi( 0, 2 )
            restr = .false.
          else
*           Use the user-supplied value --- see comments on f_lo, f_hi
            lo = abs( t_lo )
            hi = abs( t_hi )
            restr = .true.
          end if
          go_bin = n_bin( 2 )

        else if( t_mssn .eq. 131072 ) then
*         Flux Density; there must be user-supplied boundaries
          lo = abs( t_lo )
          if( t_lo .lt. 0.0 ) then
*           Flux density per Angstrom
            hi = -1
          else
*           Flux density per keV
            hi = +1
          end if
          go_bin = 16
          n_band( 2 ) = 0

        else
*         Flux; there must be user-supplied boundaries
          lo = abs( t_lo )
          hi = abs( t_hi )
          go_bin = 16
          n_band( 2 ) = 0
        end if

        t_integ = PMS_INTEG( t_mssn, t_unit, lo, hi,
     &                                e( 1, 2 ), a( 1, 0, 2 ), go_bin )
        results( 0 ) = t_integ * in_rate / f_integ
        call PMS_WTINT
     &       ( t_mssn, t_dtct, t_fltr, t_unit, t_lo, t_hi, results, 2 )
        do b = 1, n_band( 2 )
*          OLD: results( 1: n ) always in full band
*          NEW: results( 1: n ) can be in user-specified band
*          lo = sv_lo( b, 2 )
*          hi = sv_hi( b, 2 )
          t_integ = PMS_INTEG( t_mssn, t_unit, lo, hi,
     &                                e( 1, 2 ), a( 1, b, 2 ), go_bin )
          results( b ) = t_integ * in_rate / f_integ
        end do
        if( t_mssn .gt. 0 .and. t_mssn .lt. 8132 ) then
          call PMS_SPECL( t_mssn, t_dtct, t_fltr, restr,
     &                                            results, n_band( 2 ) )
        end if

        end

*+PMS_INTEG
        real function PMS_INTEG( mission, unit, lo, hi, e, a, n_bin )

        implicit none

        real kev2erg
        parameter( kev2erg = 1.60218e-09 )

        include 'pimms.inc'
        include 'pms_index.inc'

        integer mission, unit
        real lo, hi
        integer n_bin
        real e( n_bin ), a( n_bin )

*       Description:
*         Multi-component version, April 2000.  Parcels out the actual
*         job to different subroutines, such as SIMP2K.
*
*       Arguments:
*         mission      (i) : mission number (or 0 for flux)
*         unit         (i) : unit code (for flux)
*         lo, hi       (i) : integration boundaries
*         e, a         (i) : effective area curve (energy, area)
*         n_bin        (i) : number of points for above
*         <PMS_INTEG>  (r) : Result of integration
*
*       Dependencies:
*         SIMP2K, PMS_ADMDL (integration routines)
*
*       Origin:
*         Created for PIMMS by KM; based on old EXOSAT program at MSSL
*
*       Author:
*         Koji Mukai
*
*       Modified on 2001 September 28th to deal with oversensitive
*         warning from Sun's implimentation of IEEE arithmetic. (v3.2a)
*       Modified again on 2001 October 3rd (v3.2b) to prevent a crash
*         on alphas.
*       Modified one more time on October 15th (v3.2c) to make it
*         robust when two rows with identical energies are encountered
*         (as with old XTE HEXTE files, where it was a bug, but also
*         as in EXOSAT LE files, where instrumental edges are deliberately
*         treated this way).
*         
*-PMS_INTEG

        real temp, counts, o_e, o_spec, n_e, n_spec, average, small
        real nh_save( m_mdls ), nh_gsave, zlo, zhi
        real work, zpeak
        integer i, m, abs_unit

        real SPEC, ESPEC, SIMP2K, PMS_ADMDL, PMS_MDLAR, TRANMM
        external SPEC, ESPEC

        small = 1.0e-30
        zlo = lo * z1
        zhi = hi * z1
        if( mission .eq. 0 ) then
*         Flux
          if( unit .lt. 0 ) then
*           Unabsorbed flux is required
            do m = 1, n_comp
              nh_save( m ) = nh( m )
              nh( m ) = 0.0
            end do
            nh_gsave = nh_g
            nh_g = 0.0
            abs_unit = -unit
          else
            abs_unit = unit
          end if

*         Loop over model components
          temp = 0.0
          do m = 1, n_comp

            if( ( model_( m ) .ge. 1 .and. model_( m ) .le. 3 )
     &   .or. ( model_( m ) .eq. 5 .and. param( 2, m ) .gt. 0.0 ) ) then
*             Analytical (and smooth) model, use Simpson's rule integration
*             routine
              m_crrnt = m
              if( abs_unit .eq. 2 ) then
*               Unit is photons/cm/cm/s
                temp = temp + SIMP2K( SPEC, lo, hi, 1.0e-04 )
              else
*               Start off calculating in ergs/cm/cm/s
*                          (ESPEC returns kev/cm/cm/s)
                temp = temp
     &                  + SIMP2K( ESPEC, lo, hi, 1.0e-04 ) * kev2erg
              end if
              m_crrnt = 0

            else if( model_( m ) .eq. 5
     &                               .and. param( 2, m ) .eq. 0.0 ) then
*             "Gaussian" but a delta function in reality.
              if( param( 1, m ) .ge. zlo
     &                               .and. param( 1, m ) .le. zhi ) then
                if( abs_unit .eq. 2 ) then
*                   Unit is photons/cm/cm/s
                  temp = temp + norm( m ) / z1
     &                              * TRANMM( param( 1, m ), nh( m ) )
     &                              * TRANMM( param( 1, m ) / z1, nh_g )
                else
*                   Start off calculating in ergs/cm/cm/s
*                            (ESPEC returns kev/cm/cm/s)
                  temp = temp + norm( m ) / z1 * param( 1, m ) * kev2erg
     &                              * TRANMM( param( 1, m ), nh( m ) )
     &                              * TRANMM( param( 1, m ) / z1, nh_g )
                end if
              end if
            else
*             Model is a file (6-8 for plasma models, negative for others)
*             Probably not very smooth, unsuitable for Simpson's rule,
*             just add up trusting the energy steps in the file
*             ADMDL must know to handle component normalization and redshift.
              if( abs_unit .eq. 2 ) then
*               mode=1 of PMS_ADMDL is for photon/cm/cm/s
                temp = temp + PMS_ADMDL( lo, hi, 1, m )
              else
*               mode=2 is for keV/cm/cm/s
                temp = temp + PMS_ADMDL( lo, hi, 2, m ) * kev2erg
              end if
            end if
          end do
          PMS_INTEG = temp

          if( unit .lt. 0 ) then
            do m = 1, n_comp
              nh( m ) = nh_save( m )
            end do
            nh_g = nh_gsave
          end if

        else if( mission .eq. 131072 ) then
*         Flux Density
          if( unit .lt. 0 ) then
*           Unabsorbed flux is required
            do m = 1, n_comp
              nh_save( m ) = nh( m )
              nh( m ) = 0.0
            end do
            nh_gsave = nh_g
            nh_g = 0.0
            abs_unit = -unit
          else
            abs_unit = unit
          end if

          if( abs_unit .eq. 2 ) then
*           Unit right now is photon/cm/cm/s/keV
            PMS_INTEG = SPEC( lo )
          else
*           Unit right now is ergs/cm/cm/s/keV
            PMS_INTEG = ESPEC( lo ) * kev2erg / lo * lo
          end if

          if( unit .lt. 0 ) then
*           Undo the unabsorbed flux bit
            do m = 1, n_comp
              nh( m ) = nh_save( m )
            end do
            nh_g = nh_gsave
          end if

          if( hi .eq. -1 ) then
*           Convert to per angstrom
*                   (hi should be +1 if per keV)
            PMS_INTEG = PMS_INTEG * ang_kev / ( lo * lo ) 
          end if

        else
*         Integration of model convoluted with the instrument

*         Loop over model components
          temp = 0.0
          do m = 1, n_comp
            if( ( model_( m ) .ge. 1 .and. model_( m ) .le. 3 )
     &   .or. ( model_( m ) .eq. 5 .and. param( 2, m ) .gt. 0.0 ) ) then
*             Analytical model, smooth (including true Gaussian)
              m_crrnt = m
              counts = 0.0
              i = 1
              do while( e( i ) .le. lo .and. i .lt. n_bin )
                i = i + 1
              end do
              if( i .eq. n_bin )then
                goto 120
              else if( i .eq. 1 ) then
                if( e( 1 ) .ne. lo ) then
                  call PWRITE(
     &        'Warning:: integration from energies < calibration data' )
                end if
                i = i + 1
              end if
              o_e = max( e( i - 1 ), lo )
              o_spec = SPEC( o_e )
              if( o_spec .lt. small ) o_spec = 0.0
              do while( e( i ) .le. hi .and. i .le. n_bin )
                n_e = min( e( i ), hi )
                if( n_e .gt. o_e ) then
                  small = 1.0e-36 / ( n_e - o_e )
                  n_spec = SPEC( n_e )
                  if( n_spec .lt. small ) n_spec = 0.0
                  average =
     &                   ( o_spec * a( i - 1 ) + n_spec * a( i ) ) * 0.5
                  counts = counts + average * ( n_e - o_e )
                end if
                i = i + 1
                o_e = n_e
                o_spec = n_spec
              end do
              if( i .le. n_bin ) then
                n_e = hi
                n_spec = SPEC( n_e )
                average =
     &                   ( o_spec * a( i - 1 ) + n_spec * a( i ) ) * 0.5
                counts = counts + average * ( n_e - o_e )
              else
                if( e( n_bin ) .ne. hi ) then
                  call PWRITE(
     &          'Warning:: integration to energies > calibration data' )
                end if
              end if
              temp = temp + counts
 120          continue
              m_crrnt = 0

            else if( model_( m ) .eq. 5
     &                               .and. param( 2, m ) .eq. 0.0 ) then
*             It's actually a delta function
              if( param( 1, m ) .ge. zlo
     &                               .and. param( 1, m ) .le. zhi ) then
                zpeak = param( 1, m ) / z1
                i = 1
                do while( e( i + 1 ) .lt. zpeak )
                  i = i + 1
                end do
                work = ( zpeak - e( i ) ) / ( e( i + 1 ) / e( i ) )
                average = ( 1.0 - work ) * a( i ) + work * a( i + 1 )
                temp = temp + average * norm( m ) / z1
     &             * TRANMM( param( 1, m ), nh ) * TRANMM( zpeak, nh_g )
              end if

            else
*             Model came from a file, not smooth; assume sampling of the model
*             is more critical than that for the effective area curve.
              temp = temp + PMS_MDLAR( lo, hi, m, e, a, n_bin )
            end if
          end do
          PMS_INTEG = temp
        end if

        end



*+PMS_ADMDL
        real function PMS_ADMDL( lo, hi, mode, m )

        implicit none

        include 'pimms.inc'

        real lo, hi
        integer mode, m

*       Description:
*         Simple summation of file-based models between lo and hi
*
*       Arguments:
*         lo, hi      (i) : "Integration" boundary
*         mode        (i) : 1 for photon, 2 for energy (photon*keV)
*         m           (i) : Model component number
*
*       Dependencies:
*         None
*
*       Origin:
*         A quick cure by KM as Simpson's rule integration takes way too
*         long for spikey files like the Raymond-Smith.
*
*       Author:
*         Koji Mukai (1993 Mar 23) Original version
*-PMS_ADMDL
        real sum, factor, e_av, f_av
        integer fm, to, j
        real TRANMM

        PMS_ADMDL = 0.0
        j = beg_m( m )
        do while( e_in( j ) / z1 .lt. lo .and. j .lt. fin_m( m ) )
          j = j + 1
        end do
        if( j .eq. fin_m( m ) ) return
        if( j .eq. beg_m( m ) ) then
          call PWRITE(
     &             'Warning:: integration from energies < model table' )
          j = beg_m( m ) + 1
        end if
        fm = j
        do while( e_in( j ) / z1 .le. hi .and. j .lt. fin_m( m ) )
          j = j + 1
        end do
        if( j .eq. fin_m( m ) .and. e_in( j ) / z1 .lt. hi ) then
          call PWRITE(
     &               'Warning:: integration to energies > model table' )
        end if
        to = j

        sum = 0.0
        do j = fm, to
          e_av = ( e_in( j ) + e_in( j - 1 ) ) * 0.5
          f_av = norm( m ) * ( f_in( j ) + f_in( j - 1 ) ) * 0.5
     &                  * TRANMM( e_av, nh ) * TRANMM( e_av / z1, nh_g )
          factor = 1.0
          if( mode .eq. 2 ) then
            factor = factor * e_av / z1
          end if
          sum = sum + factor * f_av * ( e_in( j ) - e_in( j - 1 ) )
        end do
        
        PMS_ADMDL = sum / z1

        end


*+PMS_MDLAR
        real function PMS_MDLAR( lo, hi, m, e, a, n_bin )

        implicit none

        include 'pimms.inc'

        real lo, hi
        integer m
        integer n_bin
        real e( n_bin ), a( n_bin )

*       Description:
*         Counts from a file-based model as observed with a specified
*           instrument, between lo and hi.
*
*       Arguments:
*         lo, hi      (i) : "Integration" boundary
*         m           (i) : Model component number
*         e, a        (i) : The effective area curve
*         n_bin       (i) : The size of the array of the above
*
*       Dependencies:
*         None
*
*       Origin:
*         Numerical integration of model convolved with the instrument
*         response.
*
*       Author:
*         Koji Mukai
*-PMS_MDLAR

        real counts, average, o_e, n_e, ez, work
        integer i, j, n
        real TRANMM

        counts = 0.0
        i = 1
        do while( e( i ) .lt. lo .and. i .lt. n_bin )
          i = i + 1
        end do
        if( i .eq. n_bin ) then
          PMS_MDLAR = 0.0
          call PWRITE( 'ERROR:: Inconsistent Integration Boudaries' )
          return
        else if( i .eq. 1 ) then
          if( e( 1 ) .ne. lo ) then
            call PWRITE(
     &        'Warning:: integration from energies < calibration data' )
          end if
          i = i + 1
        end if
        o_e = max( e( i - 1 ), lo )
        j = beg_m( m )
        do while( e_in( j ) / z1 .lt. o_e .and. j .lt. fin_m( m ) )
          j = j + 1
        end do
        if( j .eq. fin_m( m ) ) then
          PMS_MDLAR = 0.0
          call PWRITE(
     &'ERROR:: non-overlapping energy ranges of model and calibration' )
          return
        end if
        ez = e_in( j ) / z1

        do while( e( i ) .le. hi .and. i .le. n_bin )
          n_e = min( e( i ), hi )
          average = 0.0
          n = 0
          do while( ez .lt. n_e .and. j .le. fin_m( m ) )
            ez = e_in( j ) / z1
            average = average + norm( m ) * f_in( j )
     &                    * TRANMM( e_in( j ), nh ) * TRANMM( ez, nh_g )
            n = n + 1
            j = j + 1
          end do
          if( n .eq. 0 .and. j .gt. fin_m( m ) ) then
            i = n_bin
          else if( n .eq. 0 ) then
            work = ( e_in( j ) + e_in( j - 1 ) ) * 0.5
            ez = work / z1
            average = norm( m ) * ( f_in( j ) + f_in( j - 1 ) ) * 0.5
     &                         * TRANMM( work, nh ) * TRANMM( ez, nh_g )
            counts = counts + average * ( n_e - o_e )
     &                                   * ( a( i - 1 ) + a( i ) ) * 0.5
          else
            average = average / real( n )
            counts = counts + average * ( n_e - o_e )
     &                                   * ( a( i - 1 ) + a( i ) ) * 0.5
          end if
          i = i + 1
          o_e = n_e
        end do

        if( i .le. fin_m( m ) ) then
          n_e = hi
          average = 0.0
          n = 0
          average = 0.0
          n = 0
          do while( ez .lt. n_e .and. j .le. fin_m( m ) )
            ez = e_in( j ) / z1
            average = average + norm( m ) * f_in( j )
     &                    * TRANMM( e_in( j ), nh ) * TRANMM( ez, nh_g )
            n = n + 1
            j = j + 1
          end do
          if( n .eq. 0 .and. j .gt. fin_m( m ) ) then
            i = n_bin
          else if( n .eq. 0 ) then
            work = ( e_in( j ) + e_in( j - 1 ) ) * 0.5
            ez = work / z1
            average = norm( m ) * ( f_in( j ) + f_in( j - 1 ) ) * 0.5
     &                         * TRANMM( work, nh ) * TRANMM( ez, nh_g )
            counts = counts + average * ( n_e - o_e )
     &                                   * ( a( i - 1 ) + a( i ) ) * 0.5
          else
            average = average / real( n )
            counts = counts + average * ( n_e - o_e )
     &                                   * ( a( i - 1 ) + a( i ) ) * 0.5
          end if
        else
          if( e( n_bin ) .ne. hi ) then
            call PWRITE(
     &          'Warning:: integration to energies > calibration data' )
          end if
        end if
        PMS_MDLAR = counts

        end
