*+PMS_WTSLC
        subroutine PMS_WTSLC( mission, detector, filter, unit, lo, hi )

        implicit none

        include 'pms_index.inc'

        integer mission, detector, filter, unit
        real lo, hi

*       Description:
*         Prints out information on the current selection on screen
*
*       Arguments:
*         mission      (i) : Mission of choice
*         detector     (i) : detector of choice
*         filter       (i) : filter of choice
*         unit         (i) : unit code for flux
*         lo, hi       (i) : energy range of interest
*
*       Dependencies:
*         PIMMS Index system, PWRITE
*
*       Origin:
*         Created by KM to help users
*
*       Author
*         Koji Mukai, 1993 March, first official version
*         Koji Mukai, 1997 December, Added Angstron support
*
*         Minor format changes in 2000 for Integral
*
*         Added flux density option for PIMMS v3.3, 2002 November
*
*         Added more digits for the energy/wavelength ranges, 2010 April
*
*-PMS_WTSLC

        character*80 xw_strng, string
        real lo_out, hi_out
        integer xw_len, st_len, u_len

        character*16 unit_st( 4 )
        integer unit_len( 4 )
        data unit_st / 'ergs/cm/cm/s', 'photons/cm/cm/s',
     &                           'mCrab', 'micro-Jansky' /
        data unit_len / 12, 15, 5, 12 /

        if( mission .eq. 0 ) then
*         FLUX
          if( unit .lt. 0 ) then
            if( lo .gt. 0.0 ) then
              write( xw_strng, 101 ) lo, hi, unit_st( -unit )
101           format( ' Unabsorbed flux (',
     &                               f11.3, '-', f11.3, ' keV) in ', a )
              xw_len = 51 + unit_len( -unit )
            else
              lo_out = -ang_kev / hi
              hi_out = -ang_kev / lo
              write( xw_strng, 102 ) lo_out, hi_out, unit_st( -unit )
102           format( ' Unabsorbed flux (',
     &                                 f11.5, '-', f11.5, ' A) in ', a )
              xw_len = 47 + unit_len( -unit )
            end if
          else
            if( lo .gt. 0.0 ) then
              write( xw_strng, 106 ) lo, hi, unit_st( unit )
106           format( ' Flux (', f11.3, '-', f11.3, ' keV) in ', a )
              xw_len = 40 + unit_len( unit )
            else
              lo_out = -ang_kev / hi
              hi_out = -ang_kev / lo
              write( xw_strng, 107 ) lo_out, hi_out, unit_st( unit )
107           format( ' Flux (', f11.5, '-', f11.5, ' A) in ', a )
              xw_len = 37 + unit_len( unit )
            end if
          end if

        else if( mission .eq. 131072 ) then
*         Flux Density
          if( unit .lt. 0 ) then
            u_len = unit_len( -unit )
            if( lo .gt. 0.0 ) then
              write( xw_strng, 111 ) lo, unit_st( -unit )( : u_len )
111           format( ' Unabsorbed flux density @',
     &                                     f11.3, 'keV in ', a, '/keV' )
              xw_len = 48 + u_len
            else
              lo_out = -ang_kev / lo
              write( xw_strng, 112 ) lo_out, unit_st( -unit )( : u_len )
112           format( ' Unabsorbed flux density @',
     &                                         f11.5, 'A in ', a, '/A' )
              xw_len = 44 + u_len
            end if
          else
	    u_len = unit_len( unit )
            if( lo .gt. 0.0 ) then
              write( xw_strng, 116 ) lo, unit_st( unit )( : u_len )
116           format( ' Flux density @', f11.3, 'keV in ', a, '/keV' )
              xw_len = 37 + u_len
            else
              lo_out = -ang_kev / lo
              write( xw_strng, 117 ) lo_out, unit_st( unit )( : u_len )
117           format( ' Flux Density @', f11.5, 'A in ', a, '/A' )
              xw_len = 33 + u_len
            end if
          end if

        else if( mission .eq. 65536 ) then
*         Model Normalization
          xw_strng = 'Model normalization'
        else
          call PDX_STRNG( mission, detector, filter, string, st_len )
          call UPCASE( string )
          if( st_len .lt. 0 ) then
            xw_strng
     &        = 'ERROR in PMS_WTSLC:: instrument name cannot be found'
            xw_len = 53
          else
            if( lo .eq. 0.0 .and. hi .eq. 0.0 ) then
              write( xw_strng, 121 ) string( : st_len )
121           format( ' Count rate in ', a )
              xw_len = 16 + st_len
            else if( lo .gt. 0.0 ) then
              write( xw_strng, 141 ) string( : st_len ), lo, hi
141           format( ' Count rate in ', a,
     &                             ' (', 0p, f9.3, '-', f9.3, ' keV)' )
              xw_len = 41 + st_len
            else
              lo_out = -ang_kev / hi
              hi_out = -ang_kev / lo
              write( xw_strng, 142 ) string( : st_len ), lo_out, hi_out
142           format( ' Count rate in ', a,
     &                             ' (', 0p, f9.5, '-', f9.5, ' A)' )
              xw_len = 39 + st_len
            end if
          end if
        end if

        call PWRITE( xw_strng( : xw_len ) )

        end


*+PMS_WTINT
        subroutine PMS_WTINT( mission, detector, filter, unit, lo, hi,
     &                                             number, case )

        implicit none

        include 'pms_index.inc'

        integer mission, detector, filter, unit, case
        real lo, hi, number

*       Write some info about the integration PIMMS just did to the screen

        character*128 xw_strng, string
        integer xw_len, st_len, u_len
        real lo_out, hi_out

        character*16 unit_st( 4 )
        integer unit_len( 4 )
        data unit_st / 'ergs/cm/cm/s', 'photons/cm/cm/s',
     &                           'mCrab', 'micro-Jansky' /
        data unit_len / 12, 15, 5, 12 /

        if( mission .eq. 0 ) then
*         FLUX
          if( unit .lt. 0 ) then
            if( case .eq. 1 ) then
              if( lo .gt. 0.0 ) then
                write( xw_strng, 101 ) lo, hi, number, unit_st( -unit )
101             format( '   and an unabsorbed flux (', f9.3, '-',
     &                             f9.3, 'keV) of ', 1p, e10.3, ' ', a )
                xw_len = 65 + unit_len( -unit )
              else
                lo_out = -ang_kev / hi
                hi_out = -ang_kev / lo
                write( xw_strng, 102 )
     &                          lo_out, hi_out, number, unit_st( -unit )
102             format( '   and an unabsorbed flux (', f11.5, '-',
     &                              f11.5, 'A) of ', 1p, e10.3, ' ', a )
                xw_len = 67 + unit_len( -unit )
              end if
            else
              if( lo .gt. 0.0 ) then
                write( xw_strng, 106 ) lo, hi, number, unit_st( -unit )
106             format( '* PIMMS predicts an unabsorbed flux (', f11.3,
     &                       '-', f11.3, 'keV) of ', 1p, e10.3, ' ', a )
                xw_len = 80 + unit_len( -unit )
              else
                lo_out = -ang_kev / hi
                hi_out = -ang_kev / lo
                write( xw_strng, 107 )
     &                          lo_out, hi_out, number, unit_st( -unit )
107             format( '* PIMMS predicts an unabsorbed flux (', f11.5,
     &                         '-', f11.5, 'A) of ', 1p, e10.3, ' ', a )
                xw_len = 78 + unit_len( -unit )
              end if
            end if
          else
            if( case .eq. 1 ) then
              if( lo .gt. 0.0 ) then
                write( xw_strng, 111 ) lo, hi, number, unit_st( unit )
111             format( '   and a flux (', f11.3, '-', f11.3,
     &                                   'keV) of ', 1p, e10.3, ' ', a )
                xw_len = 59 + unit_len( unit )
              else
                lo_out = -ang_kev / hi
                hi_out = -ang_kev / lo
                write( xw_strng, 112 )
     &                           lo_out, hi_out, number, unit_st( unit )
112             format( '   and a flux (', f11.5, '-', f11.5,
     &                                     'A) of ', 1p, e10.3, ' ', a )
                xw_len = 57 + unit_len( unit )
              end if
            else
              if( lo .gt. 0.0 ) then
                write( xw_strng, 116 ) lo, hi, number, unit_st( unit )
116             format( '* PIMMS predicts a flux (', f11.3,
     &                      '-', f11.3, 'keV) of ', 1p, e10.3, ' ', a )
                xw_len = 68 + unit_len( unit )
              else
                lo_out = -ang_kev / hi
                hi_out = -ang_kev / lo
                write( xw_strng, 117 )
     &                           lo_out, hi_out, number, unit_st( unit )
117             format( '* PIMMS predicts a flux (', f11.5,
     &                         '-', f11.5, 'A) of ', 1p, e10.3, ' ', a )
                xw_len = 66 + unit_len( unit )
              end if
            end if
          end if

        else if( mission .eq. 131072 ) then
*         Flux Density
          if( unit .lt. 0 ) then
            u_len = unit_len( -unit )
            if( case .eq. 1 ) then
              if( lo .gt. 0.0 ) then
                write( xw_strng, 201 ) lo, number,
     &                                       unit_st( -unit )( : u_len )
201             format( '   and an unabsorbed flux density @', f11.3,
     &                            'keV of ', 1p, e10.3, ' ', a, '/keV' )
                xw_len = 68 + u_len
              else
                lo_out = -ang_kev / lo
                write( xw_strng, 202 )
     &                       lo_out, number, unit_st( -unit )( : u_len )
202             format( '   and an unabsorbed flux density @', f9.5,
     &                                'A of ', 1p, e10.3, ' ', a, '/A' )
                xw_len = 62 + u_len
              end if
            else
              if( lo .gt. 0.0 ) then
                write( xw_strng, 206 ) lo, number,
     &                                       unit_st( -unit )( : u_len )
206             format( '* PIMMS predicts an unabsorbed flux density @',
     &                      f9.3, 'keV of ', 1p, e10.3, ' ', a, '/keV' )
                xw_len = 76 + u_len
              else
                lo_out = -ang_kev / lo
                write( xw_strng, 207 )
     &                       lo_out, number, unit_st( -unit )( : u_len )
207             format( '* PIMMS predicts an unabsorbed flux density @',
     &                       f9.5, 'A of ', 1p, e10.3, ' ', a, '/A' )
                xw_len = 72 + u_len
              end if
            end if
          else
            u_len = unit_len( unit )
            if( case .eq. 1 ) then
              if( lo .gt. 0.0 ) then
                write( xw_strng, 211 ) lo, number,
     &                                        unit_st( unit )( : u_len )
211             format( '   and a flux density @', f11.3,
     &                            'keV of ', 1p, e10.3, ' ', a, '/keV' )
                xw_len = 56 + u_len
              else
                lo_out = -ang_kev / lo
                write( xw_strng, 212 )
     &                        lo_out, number, unit_st( unit )( : u_len )
212             format( '   and a flux density @', f9.5,
     &                                'A of ', 1p, e10.3, ' ', a, '/A' )
                xw_len = 50 + u_len
              end if
            else
              if( lo .gt. 0.0 ) then
                write( xw_strng, 216 ) lo, number,
     &                                        unit_st( unit )( : u_len )
216             format( '* PIMMS predicts a flux density @', f9.3,
     &                            'keV of ', 1p, e10.3, ' ', a, '/keV' )
                xw_len = 64 + u_len
              else
                lo_out = -ang_kev / lo
                write( xw_strng, 217 )
     &                        lo_out, number, unit_st( unit )( : u_len )
217             format( '* PIMMS predicts a flux density @', f9.5,
     &                                'A of ', 1p, e10.3, ' ', a, '/A' )
                xw_len = 60 + u_len
              end if
            end if
          end if

        else if( mission .eq. 65536 ) then
*         Don't produce output for 'norm' instrument
          xw_strng = ' '
          xw_len = 0

        else
          call PDX_STRNG( mission, detector, filter, string, st_len )
          call UPCASE( string )
          if( st_len .lt. 0 ) then
            xw_strng
     &        = 'ERROR in PMS_WTINT:: instrument name cannot be found'
            xw_len = 53
          else
            if( lo .eq. 0.0 .and. hi .eq. 0.0 ) then
              if( case .eq. 1 ) then
                write( xw_strng, 121 ) number, string( : st_len )
121             format( '  and ', 1p, e10.3, ' cps in ', a )
                xw_len = 25 + st_len
              else
                write( xw_strng, 122 ) number, string( : st_len )
122             format( '* PIMMS predicts ', 1p, e10.3,
     &                                               ' cps with ', a )
                xw_len = 38 + st_len
              end if
            else if( lo .gt. 0.0 ) then
              if( case .eq. 1 ) then
                write( xw_strng, 141 )
     &                                number, string( : st_len ), lo, hi
141             format( '  and ', 1p, e10.3, ' cps in ', a,
     &                             ' (', 0p, f9.3, '-', f9.3, 'keV)' )
                xw_len = 50 + st_len
              else
                write( xw_strng, 142 )
     &                                number, string( : st_len ), lo, hi
142             format( '* PIMMS predicts ', 1p, e10.3, ' cps with ',
     &                            a, ' (', 0p, f9.3, '-', f9.3, 'keV)' )
                xw_len = 63 + st_len
              end if
            else
              lo_out = -ang_kev / hi
              hi_out = -ang_kev / lo
              if( case .eq. 1 ) then
                write( xw_strng, 161 )
     &                        number, string( : st_len ), lo_out, hi_out
161             format( '  and ', 1p, e10.3, ' cps in ', a,
     &                                 ' (', 0p, f9.5, '-', f9.5, 'A)' )
                xw_len = 48 + st_len
              else
                write( xw_strng, 162 )
     &                        number, string( : st_len ), lo_out, hi_out
162             format( '* PIMMS predicts ', 1p, e10.3, ' cps with ',
     &                              a, ' (', 0p, f9.5, '-', f9.5, 'A)' )
                xw_len = 61 + st_len
              end if
            end if
          end if
        end if

        if( xw_len .gt. 0 ) then
          call PWRITE( xw_strng( : xw_len ) )
        end if

        end


*+PMS_WTMDL
        subroutine PMS_WTMDL( case )

        implicit none

        include 'pimms.inc'

        integer case

*       Outputs to the screen a few simple information on the model
*          case=1, general info
*          case=2, with integration results
*       New version for the multi-component model
*       Added Galactic Nh output for z=0 but G_Nh>0, 2010 April

        character*80 xw_strng
        integer xw_len, len_mf, m

        integer LENTRIM

        do m = 1, n_comp
*         Loop over components

          if( model_( m ) .eq. 1 ) then
*           Blackbody
            if( m .eq. 1 ) then
              if( case .eq. 1 ) then
                write( xw_strng, 101 ) param( 1, m ), nh( m )
 101            format( '* Current model is BLACKBODY, kT= ',
     &                                  f8.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 64
              else
                write( xw_strng, 102 ) param( 1, m ), nh( m )
 102            format( '* For Blackbody model with kT=', f8.4,
     &                                        ' keV; NH = ', 1p, e10.3 )
                xw_len = 60
              end if
            else
              if( case .eq. 1 ) then
                write( xw_strng, 106 ) param( 1, m ), nh( m )
 106            format( '              plus BLACKBODY, kT= ',
     &                                  f8.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 64
              else
                write( xw_strng, 107 ) param( 1, m ), nh( m )
 107            format( '    + Blackbody model with kT=', f8.4,
     &                                        ' keV; NH = ', 1p, e10.3 )
                xw_len = 60
              end if
            end if

          else if( model_( m ) .eq. 2 ) then
*           Bremsstrahlung
            if( m .eq. 1 ) then
              if( case .eq. 1 ) then
                write( xw_strng, 111 ) param( 1, m ), nh( m )
 111            format( '* Current model is BREMSSTRAHLUNG, kT= ',
     &                                  f8.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 69
              else
                write( xw_strng, 112 ) param( 1, m ), nh( m )
 112            format( '* For thermal Bremsstrahlung model with kT=',
     &                                  f8.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 73
              end if
            else
              if( case .eq. 1 ) then
                write( xw_strng, 116 ) param( 1, m ), nh( m )
 116            format( '              plus BREMSSTRAHLUNG, kT= ',
     &                                  f8.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 69
              else
                write( xw_strng, 117 ) param( 1, m ), nh( m )
 117            format( '    + thermal Bremsstrahlung model with kT=',
     &                                  f8.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 73
              end if
            end if

          else if( model_( m ) .eq. 3 ) then
*           Power Law Family
            if( param( 2, m ) .eq. 0.0 ) then
*             Simple Power Law
              if( m .eq. 1 ) then
                if( case .eq. 1 ) then
                  write( xw_strng, 121 ) -param( 1, m ), nh( m )
 121              format( '* Current model is POWER LAW, Photon ',
     &                          'Index = ', f7.4, '; NH = ', 1p, e10.3 )
                  xw_len = 79
                else
                  write( xw_strng, 122 ) -param( 1, m ), nh( m )
 122              format( '* For power law model with photon index =',
     &                                      f7.4, '; NH = ', 1p, e10.3 )
                  xw_len = 66
                end if
              else
                if( case .eq. 1 ) then
                  write( xw_strng, 126 ) -param( 1, m ), nh( m )
 126              format( '              plus POWER LAW, Photon ',
     &                          'Index = ', f7.4, '; NH = ', 1p, e10.3 )
                  xw_len = 79
                else
                  write( xw_strng, 127 ) -param( 1, m ), nh( m )
 127              format( '    + power law model with photon index =',
     &                                      f7.4, '; NH = ', 1p, e10.3 )
                  xw_len = 66
                end if
              end if

            else if( param( 3, m ) .eq. 0.0 ) then
*             Cut-off Power Law
              if( m .eq. 1 ) then
                if( case .eq. 1 ) then
                  write( xw_strng, 131 )
     &                            -param( 1, m ), param( 2, m ), nh( m )
 131              format( '* Current model is CUTOFF PL, Index = ',
     &               f5.2, ', Ecut ', f7.2, ' keV; NH = ', 1p, e10.3 )
                  xw_len = 79
                else
                  write( xw_strng, 132 )
     &                            -param( 1, m ), param( 2, m ), nh( m )
 132              format( '* For cutoff pl model with index =', f5.2,
     &                     ', Ecut ', f7.2, ' keV; NH = ', 1p, e10.3 )
                  xw_len = 75
                end if
              else
                if( case .eq. 1 ) then
                  write( xw_strng, 136 )
     &                            -param( 1, m ), param( 2, m ), nh( m )
 136              format( '              plus CUTOFF PL, Index = ',
     &               f5.2, ', Ecut ', f7.2, ' keV; NH = ', 1p, e10.3 )
                  xw_len = 79
                else
                  write( xw_strng, 137 )
     &                            -param( 1, m ), param( 2, m ), nh( m )
 137              format( '    + cutoff pl model with index =', f5.2,
     &                     ', Ecut ', f7.2, ' keV; NH = ', 1p, e10.3 )
                  xw_len = 75
                end if
              end if

            else
*             Power Law with High-Energy (exponential) cut-off
              if( m .eq. 1 ) then
                if( case .eq. 1 ) then
                  call PWRITE( '* Current model is power-law with ' //
     &                                           'high-energy cut-off' )
                  write( xw_strng, 141 ) -param( 1, m ), param( 2, m ),
     &                                            param( 3, m ), nh( m )
 141              format( '  Index = ', f5.2, ', Ecut ', f7.2,
     &       ' keV, E(e-folding) ', f7.2, ' keV; NH = ', 1p, e10.3 )
                  xw_len = 76
                else
                  call PWRITE( '* For power-law model with ' //
     &                                      'high-energy cut-off with' )
                  write( xw_strng, 141 ) -param( 1, m ), param( 2, m ),
     &                                            param( 3, m ), nh( m )
                  xw_len = 78
                end if
              else
                if( case .eq. 1 ) then
                  call PWRITE( '              plus power-law with ' //
     &                                           'high-energy cut-off' )
                  write( xw_strng, 141 ) -param( 1, m ), param( 2, m ),
     &                                            param( 3, m ), nh( m )
                  xw_len = 78
                else
                  call PWRITE( '    + power-law model with ' //
     &                                      'high-energy cut-off with' )
                  write( xw_strng, 141 ) -param( 1, m ), param( 2, m ),
     &                                            param( 3, m ), nh( m )
                  xw_len = 78
                end if
              end if
            end if

          else if( model_( m ) .eq. 4 ) then
*           Old Raymond-Smith model --- should not be here
            if( m .eq. 1 ) then
              if( case .eq. 1 ) then
                write( xw_strng, 151 )
     &                             param( 1, m ), param( 2, m ), nh( m )
 151            format( '* Current model is RAYMOND SMITH, ',
     &        'kT=', f7.4, ' keV (logT=', f5.2, '); NH = ', 1p, e10.3 )
                xw_len = 79
              else
                write( xw_strng, 152 )
     &                             param( 1, m ), param( 2, m ), nh( m )
 152            format( '* For Raymond Smith model with kT=', f7.4,
     &                      ' keV (logT=', f5.2, '); NH = ', 1p, e10.3 )
                xw_len = 76
              end if
            else
              if( case .eq. 1 ) then
                write( xw_strng, 156 )
     &                             param( 1, m ), param( 2, m ), nh( m )
 156            format( '              plus RAYMOND SMITH, ',
     &        'kT= ', f8.4, ' keV (logT=', f3.1, '); NH = ', 1p, e10.3 )
                xw_len = 79
              else
                write( xw_strng, 157 )
     &                             param( 1, m ), param( 2, m ), nh( m )
 157            format( '    + Raymond Smith model with kT=', f7.4,
     &                      ' keV (logT=', f5.2, '); NH = ', 1p, e10.3 )
                xw_len = 75
              end if
            end if

          else if( model_( m ) .eq. 5 ) then
*           Gaussian
            if( m .eq. 1 ) then
              if( case .eq. 1 ) then
                write( xw_strng, 161 )
     &                             param( 1, m ), param( 2, m ), nh( m )
 161            format( '* Current model is GAUSSIAN, E=',
     &            f8.4, ' keV; sigma=', f7.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 79
              else
                write( xw_strng, 162 )
     &                             param( 1, m ), param( 2, m ), nh( m )
 162            format( '* For Gaussian model with E=', f8.4,
     &                  ' keV; sigma=', f7.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 76
              end if
            else
              if( case .eq. 1 ) then
                write( xw_strng, 166 )
     &                             param( 1, m ), param( 2, m ), nh( m )
 166            format( '              plus GAUSSIAN, E=',
     &            f8.4, ' keV; sigma=', f7.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 79
              else
                write( xw_strng, 167 )
     &                             param( 1, m ), param( 2, m ), nh( m )
 167            format( '    + Gaussian model with E=', f8.4,
     &                  ' keV; sigma=', f7.4, ' keV; NH = ', 1p, e10.3 )
                xw_len = 76
              end if
            end if

          else if( model_( m ) .eq. 6 ) then
*           New Raymond-Smith model
            if( m .eq. 1 ) then
              if( case .eq. 1 ) then
                write( xw_strng, 171 )
 171            format(
     &                '* Current model is PLASMA (Raymond-Smith) with' )
                xw_len = 46
              else
                write( xw_strng, 172 )
 172            format( '* For PLASMA (Raymond-Smith) model with' )
                xw_len = 39
              end if
            else
              if( case .eq. 1 ) then
                write( xw_strng, 173 )
 173            format(
     &                '              plus PLASMA (Raymond-Smith) with' )
                xw_len = 46
              else
                write( xw_strng, 174 )
 174            format( '    + PLASMA (Raymond-Smith) model with' )
                xw_len = 39
              end if
            end if
            call PWRITE( xw_strng )
            write( xw_strng, 179 )
     &             param( 1, m ), param( 2, m ), param( 3, m ), nh( m )
 179        format( '                        kT=', f7.4, 'keV (logT=',
     &                   f5.2, '), Abund=', f3.1, '; NH = ', 1p, e10.3 )
            xw_len = 78

          else if( model_( m ) .eq. 7 ) then
*           New APEC model
            if( m .eq. 1 ) then
              if( case .eq. 1 ) then
                write( xw_strng, 181 )
 181            format( '* Current model is PLASMA (APEC) with' )
                xw_len = 37
              else
                write( xw_strng, 182 )
 182            format( '* For PLASMA (APEC) model with' )
                xw_len = 30
              end if
            else
              if( case .eq. 1 ) then
                write( xw_strng, 183 )
 183            format( '              plus PLASMA (APEC) with' )
                xw_len = 37
              else
                write( xw_strng, 184 )
 184            format( '    + PLASMA (APEC) model with' )
                xw_len = 30
              end if
            end if
            call PWRITE( xw_strng )
            write( xw_strng, 179 )
     &             param( 1, m ), param( 2, m ), param( 3, m ), nh( m )
            xw_len = 78

          else if( model_( m ) .eq. 8 ) then
*           New Mekal model
            if( m .eq. 1 ) then
              if( case .eq. 1 ) then
                write( xw_strng, 191 )
 191            format( '* Current model is PLASMA (mekal) with' )
                xw_len = 38
              else
                write( xw_strng, 192 )
 192            format( '* For PLASMA (mekal) model with' )
                xw_len = 31
              end if
            else
              if( case .eq. 1 ) then
                write( xw_strng, 193 )
 193            format( '              plus PLASMA (mekal) with' )
                xw_len = 38
              else
                write( xw_strng, 194 )
 194            format( '    + PLASMA (mekal) model with' )
                xw_len = 31
              end if
            end if
            call PWRITE( xw_strng )
            write( xw_strng, 179 )
     &             param( 1, m ), param( 2, m ), param( 3, m ), nh( m )
            xw_len = 78

          else if( model_( m ) .eq. -1 ) then
*           Model read in from a file
            len_mf = LENTRIM( mf_name( m ) )
            if( nh( m ) .eq. 0.0 ) then
*             ...with no additional Nh
              if( m .eq. 1 ) then
                if( case .eq. 1 ) then
                  write( xw_strng, 201 ) mf_name( m )( : len_mf )
 201              format( '* Model from ', a )
                  xw_len = len_mf + 14
                else
                  write( xw_strng, 202 ) mf_name( m )( : len_mf )
 202              format( '* For model ', a )
                  xw_len = len_mf + 13
                end if
              else
                if( case .eq. 1 ) then
                  write( xw_strng, 206 ) mf_name( m )( : len_mf )
 206              format( '        plus ', a )
                  xw_len = len_mf + 14
                else
                  write( xw_strng, 207 ) mf_name( m )( : len_mf )
 207              format( '          + ', a )
                  xw_len = len_mf + 13
                end if
              end if
            else
*             There is additional Nh
              if( m .eq. 1 ) then
                if( case .eq. 1 ) then
                  write( xw_strng, 211 )
     &                                 mf_name( m )( : len_mf ), nh( m )
 211              format( '* Model from ', a, '; NH = ', 1p, e10.3 )
                  xw_len = len_mf + 31
                else
                  write( xw_strng, 212 )
     &                                 mf_name( m )( : len_mf ), nh( m )
 212              format( '* For model ', a, '; NH = ', 1p, e10.3 )
                  xw_len = len_mf + 30
                end if
              else
                if( case .eq. 1 ) then
                  write( xw_strng, 216 )
     &                                 mf_name( m )( : len_mf ), nh( m )
 216              format( '        plus ', a, '; NH = ', 1p, e10.3 )
                  xw_len = len_mf + 31
                else
                  write( xw_strng, 217 )
     &                                 mf_name( m )( : len_mf ), nh( m )
 217              format( '          + ', a, '; NH = ', 1p, e10.3 )
                  xw_len = len_mf + 30
                end if
              end if
            end if
          end if

          call PWRITE( xw_strng( : xw_len ) )

          if( m .ge. 2 ) then
            if( model_( m ) .eq. 5 ) then
              write( xw_strng, 301 ) m_rat( m )
 301          format( '                          (Eq.W=', f8.4, ' eV)' )
              xw_len = 44
            else
              write( xw_strng, 302 ) m_rat( m ), m_e( m )
 302          format( '             (', f8.4, ' times component 1 at ',
     &                                                  f10.4, ' keV)' )
              xw_len = 59
            end if
            call PWRITE( xw_strng( : xw_len ) )
          end if
        end do

        if( z1 .ne. 1.0 ) then
          if( nh_g .eq. 0.0 ) then
            write( xw_strng, 401 ) z1 - 1.0
 401	    format( ' ...redshifted with z=', f8.4 )
            xw_len = 30
          else
            write( xw_strng, 402 ) z1 - 1.0, nh_g
 402	    format( ' ...redshifted with z=', f8.4,
     &               ' and a Galactic Nh=', 1p, e10.3 )
            xw_len = 59
          end if
          call PWRITE( xw_strng( : xw_len ) )
        else if( nh_g .gt. 0.0 ) then
          write( xw_strng, 403 ) nh_g
 403      format( ' and a Galactic Nh=', 1p, e10.3 )
          xw_len = 29
          call PWRITE( xw_strng( : xw_len ) )
        end if

        end
