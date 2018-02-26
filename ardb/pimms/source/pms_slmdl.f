        subroutine PMS_SLMDL( flag )

        implicit none

        integer m_rs, m_abund
        real t_kev
        parameter( m_rs = 128, m_abund = 16 )
        parameter( t_kev = 1.16048e+07 )

        include 'pimms.inc'

*       Description:
*         Reads in model definition from user via ARK-FCI
*
*       Arguments:
*
*       Dependencies:
*         ARK-FCI, LENTRIM
*
*       Origin:
*         Created by KM
*
*       Author:
*         Koji Mukai, 1993 March, first official XPI version
*         Koji Mukai, 1993 March, first official ARK-FCI version
*         Koji Mukai, 1993 May, modified cope with increased number
*                     of Raymond-Smith models and 2nd parameter for
*                     power law; also some bug fixes etc.
*         Koji Mukai, 1997 December.  Allows positive power law index
*                     Power law now defined as E**( -param( 1 ) )
*         Koji Mukai, 2000 April.  Multi-component version
*                     Corrected a power law bug with v3.1a release.
*         Koji Mukai, 2002 January.  Corrected a bug in initialization
*                     of 'highest', which was preventing re-use of
*                     same memory areas for file-based models.
*         Koji Mukai, 2004 August.  Corrected a bug so that a model
*                     file used as 2nd-nth component can have its own Nh
*         Koji Mukai, 2006 August.  Protecting against core dumps
*                     on Linux (etc.) compilers due to internal read
*         Koji Mukai, 2010 January.  A major update to turn off old
*                     "raymondsmith" model and allow the new "plasma"
*                     model, which could be APEC, mekal, or R-S.
*-PMS_SLMDL

        integer flag

        character*128 one_line
        character*80 xw_strng
        character*60 model_st, unit_st, last_st
        integer len_model, i, j, mf_unit, unit_code
        integer num_c, num_f, num_i, m, nc_tgt, mm
        integer np_num, np_beg, next, extra
        integer lflag
        integer n_abu, n_logt, k, l, l_proot
        real logt( m_rs ), kt( m_rs ), av_t, val_1, val_m
        real abun( m_abund )
        character*8 p_root
        character*3 logt_str( m_rs )
        character*2 ab_str( m_abund )
        character*32 file_st
        logical warn_rs, grid_warn

        integer LENTRIM
        real SPEC

        character*32 bb_st, tb_st, pl_st, rs_st, ga_st

        include 'sitespec.inc'

        data bb_st, tb_st / 'blackbody', 'bremsstrahlung' /
        data pl_st, rs_st / 'powerlaw', 'plasma' /
        data ga_st / 'gaussian' /

        
        call INQ_PARAM( num_c, num_f, num_i )
        if( num_c .ge. 1 ) then
          call GET_PARAM_C( 1, model_st )
          call LOCASE( model_st )
          len_model = LENTRIM( model_st )
        end if
        if( model_st( : len_model ) .eq. '?'
     &                   .or. model_st( : len_model ) .eq. 'help' ) then
*         No model names were given
          call PWRITE( 'PIMMS knows about the following' )
          call PWRITE( '     Blackbody (BB)' )
          call PWRITE( '     Powerlaw (PL)' )
          call PWRITE( '     (Thermal) Bremsstrahlung (TB)' )
          call PWRITE( '     Plasma' )
          call PWRITE(
     &     '  For the above, enter name, temperature/slope and Nh' )
          call PWRITE( '                 e.g., model bb 2.4 2e21' )
          call PWRITE( 'also Gaussian (GA) model is available' )
          call ARKOPN( mf_unit, mdir_name, 'model.idx', 'idx',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
          if( flag .lt. 0 ) then
            call PWRITE( 'SEVERE ERROR:: failed to open model.idx' )
            stop
          end if
          call PWRITE(
     &    'PIMMS also knows about following precalculated models:' )
          do while( .true. )
            read( mf_unit, '(a)', end = 100 ) one_line
            call PWRITE( '  ' // one_line( : LENTRIM( one_line ) ) )
          end do
100       continue
          close( mf_unit )
          call PWRITE(
     &                '  For these, type MODEL <modelname> [<Nh>]' )
          call PWRITE(
     &            'PIMMS can also use any Ascii files containing:' )
          call PWRITE( '       energy (kev), flux (photons/cm/cm/s) '
     &                                          // 'pairs as model' )
          call PWRITE( '  For this, type <filename> [<Nh>]' )
          call PWRITE( ' ' )
          call PWRITE( 'For a full usage of the model command, '
     &                     // '(including multi-component models),' )
          call PWRITE( 'please consult the User''s Guide' )
          flag = -1
          return
        end if

        if( num_c .gt. 1 ) then
          call GET_PARAM_C( num_c, last_st )
          call LOCASE( last_st )
          if( last_st .eq. 'z' ) then
*         Read redshift and (optional) Galactic NH
            call PARAM_SUB( num_c, np_num, np_beg )
            if( np_num .eq. 0 ) then
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter redshift >', one_line )
                call RD_REAL( one_line, z1, lflag )
              end do
              z1 = 1.0 + z1
              nh_g = 0.0
            else
              call GT2_PARAM_N( np_beg, z1 )
              z1 = 1.0 + z1
              if( np_num .ge. 2 ) then
                if( np_num .ge. 3 ) then
                  call PWRITE(
     &             'WARNING:: Extra numerical parameter(s) - ignoring' )
                end if
                call GT2_PARAM_N( np_beg + 1, nh_g )
              else
                nh_g = 0.0
              end if
            end if
            nc_tgt = num_c - 1
          else
            nc_tgt = num_c
            z1 = 1.0
            nh_g = 0.0
          end if
        else
          nc_tgt = 1
          z1 = 1.0
          nh_g = 0.0
        end if

*       nc_tgt is the number of character strings to be parsed
*       (after taking care of the last z, if present, which signifies
*       the user specified the redshift)

*       m is the model component index, mm is the index for the
*       next character string parameter to be parsed as the model name
*       PARAM_SUB is used to figure out that the character string mm
*       is followed by np_num numerical parameters starting with index
*       np_beg

        m = 1
        mm = 1
        extra = 0
        warn_rs = .true.
        do while( mm .le. nc_tgt .and. m .le. m_mdls )
          if( mm .gt. 1 ) then
            call GET_PARAM_C( mm, model_st )
            call LOCASE( model_st )
            len_model = LENTRIM( model_st )
          end if
          call PARAM_SUB( mm, np_num, np_beg )
          if( np_num .lt. 0 ) then
*           FCI Error
            call PWRITE( 'ERROR parsing parameters:: FCI is confused' )
            flag = np_num
            return
          end if

          if( ( len_model .eq. 2 .and. model_st( : 2 ) .eq. 'bb' )
     &   .or. ( len_model .eq. 2 .and. model_st( : 2 ) .eq. 'bl' )
     &     .or. model_st( : len_model ) .eq. bb_st( : len_model ) ) then 
*           BLACKBODY model
            if( np_num .eq. 0 ) then
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter blackbody temperature (keV) >',
     &                                                        one_line )
                call RD_REAL( one_line, param( 1, m ), lflag )
              end do
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter Nh >', one_line )
                call RD_REAL( one_line, nh( m ), lflag )
              end do
            else if( np_num .eq. 1 ) then
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter Nh >', one_line )
                call RD_REAL( one_line, nh( m ), lflag )
              end do
            else
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              call GT2_PARAM_N( np_beg + 1, nh( m ) )
            end if
            model_( m ) = 1
            param( 2, m ) = 0.0
            param( 3, m ) = 0.0
            next = 2

          else if( ( len_model .eq. 2 .and. model_st( : 2 ) .eq. 'tb' )
     &        .or. ( len_model .eq. 2 .and. model_st( : 2 ) .eq. 'br' )
     &     .or. model_st( : len_model ) .eq. tb_st( : len_model ) ) then 
*           (thermal) BREMSSTRAHLUNG model
            if( np_num .eq. 0 ) then
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR(
     &            'Enter bremsstrahlung temperature (keV) >', one_line )
                call RD_REAL( one_line, param( 1, m ), lflag )
              end do
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter Nh >', one_line )
                call RD_REAL( one_line, nh( m ), lflag )
              end do
            else if( np_num .eq. 1 ) then
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter Nh >', one_line )
                call RD_REAL( one_line, nh( m ), lflag )
              end do
            else
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              call GT2_PARAM_N( np_beg + 1, nh( m ) )
            end if
            model_( m ) = 2
            param( 2, m ) = 0.0
            param( 3, m ) = 0.0
            next = 2

          else if( ( len_model .eq. 2 .and. model_st( : 2 ) .eq. 'pl' )
     &     .or. model_st( : len_model ) .eq. pl_st( : len_model ) ) then
*             POWERLAW model
            if( mm .eq. 1 )then
              next = 0
            else
              next = 2
            end if
            if( np_num .eq. 0 ) then
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter photon index >', one_line )
                call RD_REAL( one_line, param( 1, m ), lflag )
              end do
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter Nh >', one_line )
                call RD_REAL( one_line, nh( m ), lflag )
              end do
              param( 2, m ) = 0.0
              param( 3, m ) = 0.0
              next = 2
            else if( np_num .eq. 1 ) then
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter Nh >', one_line )
                call RD_REAL( one_line, nh( m ), lflag )
              end do
              param( 2, m ) = 0.0
              param( 3, m ) = 0.0
              next = 2
            else if( np_num - next .le. 2 ) then
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              call GT2_PARAM_N( np_beg + 1, nh( m ) )
              param( 2, m ) = 0.0
              param( 3, m ) = 0.0
              next = 2
            else if( np_num - next .eq. 3 ) then
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              call GT2_PARAM_N( np_beg + 1, param( 2, m ) )
              call GT2_PARAM_N( np_beg + 2, nh( m ) )
              param( 3, m ) = 0.0
              next = 3
            else
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              call GT2_PARAM_N( np_beg + 1, param( 2, m ) )
              call GT2_PARAM_N( np_beg + 2, param( 3, m ) )
              call GT2_PARAM_N( np_beg + 3, nh( m ) )
              next = 4
            end if
            param( 1, m ) = -param( 1, m )
            model_( m ) = 3

          else if( ( len_model .eq. 2 .and. model_st( : 2 ) .eq. 'ga' )
     &     .or. model_st( : len_model ) .eq. ga_st( : len_model ) ) then
*           GAUSSIAN model
            next = 2
            if( m .gt. 1 ) then
              nh( m ) = nh( 1 )
            else
              nh( m ) = 0.0
            end if
            if( np_num .eq. 0 ) then
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter line center energy >', one_line )
                call RD_REAL( one_line, param( 1, m ), lflag )
              end do
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter physical width (keV) >',
     &                                                       one_line )
                call RD_REAL( one_line, param( 2, m ), lflag )
              end do
            else if( np_num .eq. 1 ) then
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter physical width (keV) >',
     &                                                       one_line )
                call RD_REAL( one_line, param( 2, m ), lflag )
              end do
            else if( np_num .eq. 2 ) then
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              call GT2_PARAM_N( np_beg + 1, param( 2, m ) )
            else
              call GT2_PARAM_N( np_beg, param( 1, m ) )
              call GT2_PARAM_N( np_beg + 1, param( 2, m ) )
              if( ( m .eq. 1 .and. np_num .eq. 3 )
     &                      .or. ( m .gt. 1 .and. np_num .eq. 4 ) ) then
                call GT2_PARAM_N( np_beg + 2, nh( m ) )
                next = 3
              end if
            end if
            model_( m ) = 5
            param( 3, m ) = 0.0

*          else if( ( len_model .eq. 2 .and. model_st( : 2 ) .eq. 'rs' )
*     &     .or. model_st( : len_model ) .eq. rs_st( : len_model ) ) then
**           RAYMONDSMITH model
*            OLD CODE DELETED

          else if(
     &          model_st( : len_model ) .eq. rs_st( : len_model ) ) then
*           PLASMA model
            if( p_model .eq. 0 ) then
              call ARKOPN( mf_unit, mdir_name, 'rs.idx', 'idx',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
            else if( p_model .eq. 1 ) then
              call ARKOPN( mf_unit, mdir_name, 'apec.idx', 'idx',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
            else if( p_model .eq. 2 ) then
              call ARKOPN( mf_unit, mdir_name, 'mekal.idx', 'idx',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
            else
              call PWRITE( 'ERROR:: Plasma code unrecognized' )
              flag = -101
              return
            end if
            if( flag .lt. 0 ) then
              call PWRITE( 'ERROR: Failed to opne plasma index file' )
              flag = -102
              return
            end if
            read( mf_unit, * ) p_root
            l_proot = LENTRIM( p_root )
            read( mf_unit, * ) n_abu
            do k = 1, n_abu
              read( mf_unit, * ) abun( k ), ab_str( k )
            end do
            read( mf_unit, * ) n_logt
            do l = 1, n_logt
              read( mf_unit, * ) logt( l ), logt_str( l )
              kt( l ) = 10.0 ** logt( l ) / t_kev
            end do
            close( mf_unit )

            if( mm .lt. nc_tgt ) then
              call GET_PARAM_C( mm + 1, unit_st )
              call LOCASE( unit_st )
              if( unit_st .eq. 'logt' ) then
                unit_code = 2
*               There was a unit string - temperature should have preceded it
                if( np_num .eq. 1 ) then
                  call GT2_PARAM_N( np_beg, param( 1, m ) )
                else if( np_num .eq. 0 ) then
                  one_line = ' '
                  lflag = 0
                  do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                    call GET_PSTR(
     &                   'Enter log10(plasma temperature) >', one_line )
                    call RD_REAL( one_line, param( 1, m ), lflag )
                  end do
                else
                  call PWRITE( 'ERROR:: incorrect plasma model syntax' )
                  flag = -np_num
                  return
                end if
                mm = mm + 1
                call PARAM_SUB( mm, np_num, np_beg )
                if( np_num .lt. 0 ) then
*                 FCI Error
                  call PWRITE(
     &                    'ERROR parsing parameters:: FCI is confused' )
                  flag = np_num
                  return
                end if
                np_num = np_num + 1
                np_beg = np_beg - 1
                extra = extra + 1
              else if( unit_st .eq. 'kev' ) then
                unit_code = 1
*               There was a unit string - temperature should have preceded it
                if( np_num .eq. 1 ) then
                  call GT2_PARAM_N( np_beg, param( 1, m ) )
                else if( np_num .eq. 0 ) then
                  one_line = ' '
                  lflag = 0
                  do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                    call GET_PSTR(
     &                   'Enter plasma temperature (keV) >', one_line )
                    call RD_REAL( one_line, param( 1, m ), lflag )
                  end do
                else
                  call PWRITE( 'ERROR:: incorrect plasma model syntax' )
                  flag = -np_num
                  return
                end if
                mm = mm + 1
                call PARAM_SUB( mm, np_num, np_beg )
                if( np_num .lt. 0 ) then
*                 FCI Error
                  call PWRITE(
     &                    'ERROR parsing parameters:: FCI is confused' )
                  flag = np_num
                  return
                end if
                np_num = np_num + 1
                np_beg = np_beg - 1
                extra = extra + 1
              else
*               Forget I even checked.  Assume next string is a model name
                unit_code = 1
                if( np_num .ge. 1 ) then
                  call GT2_PARAM_N( np_beg, param( 1, m ) )
                else
                  one_line = ' '
                  lflag = 0
                  do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                    call GET_PSTR(
     &                   'Enter plasma temperature (keV) >', one_line )
                    call RD_REAL( one_line, param( 1, m ), lflag )
                  end do
                end if
              end if
            else
              unit_code = 1
              if( np_num .ge. 1 ) then
                call GT2_PARAM_N( np_beg, param( 1, m ) )
              else
                one_line = ' '
                lflag = 0
                do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                  call GET_PSTR(
     &                   'Enter plasma temperature (keV) >', one_line )
                  call RD_REAL( one_line, param( 1, m ), lflag )
                end do
              end if
            end if

            if( np_num .le. 1 ) then
              one_line = ' '
              lflag = 0
              do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                call GET_PSTR( 'Enter abundance (solar) >', one_line )
                call RD_REAL( one_line, param( 3, m ), lflag )
              end do
            else
              call GT2_PARAM_N( np_beg + 1, param( 3, m ) )
            end if

            grid_warn = .false.
            if( unit_code .eq. 1 ) then
              do l = 1, n_logt
                if( param( 1, m ) .eq. kt( l ) ) goto 530
              end do
            else if( unit_code .eq. 2 ) then
              do l = 1, n_logt
                if( param( 1, m ) .eq. logt( l ) ) goto 530
              end do
            else
              stop 'Confusion in unit_code'
            end if
            grid_warn = .true.

            if( unit_code .eq. 1 ) then
              l = 1
              do while( l .lt. n_logt )
                av_t = ( kt( l ) + kt( l + 1 ) ) * 0.5
                if( param( 1, m ) .lt. av_t ) goto 530
                l = l + 1
              end do
            else
              l = 1
              do while( l .lt. n_logt )
                av_t = ( logt( l ) + logt( l + 1 ) ) * 0.5
                if( param( 1, m ) .lt. av_t ) goto 530
                l = l + 1
              end do
            end if
530         continue
*           The cases where temperature exactly matched rejoin here

            do k = 1, n_abu
                if( param( 3, m ) .eq. abun( k ) ) goto 540
            end do
            grid_warn = .true.
            k = 1
            do while( k .lt. n_abu )
              av_t = ( abun( k ) + abun( k + 1 ) ) * 0.5
              if( param( 3, m ) .lt. av_t ) goto 540
              k = k + 1
            end do
540         continue

            if( grid_warn ) then
              if( warn_rs ) then
                if( p_model .eq. 0 ) then
                  write( xw_strng, 551 ) n_logt, n_abu
551               format( 'NOTE: This version of PIMMS has a grid of ',
     &                    i3, 'x', i2, ' grid of Raymond-Smith models' )
                  call PWRITE( xw_strng( : 77 ) )
                else if( p_model .eq. 1 ) then
                  write( xw_strng, 552 ) n_logt, n_abu
552               format( 'NOTE: This version of PIMMS has a grid of ',
     &                    i3, 'x', i2, ' grid of APEC models' )
                  call PWRITE( xw_strng( : 68 ) )
                else if( p_model .eq. 2 ) then
                  write( xw_strng, 553 ) n_logt, n_abu
553               format( 'NOTE: This version of PIMMS has a grid of ',
     &                    i3, 'x', i2, ' grid of mekal models' )
                  call PWRITE( xw_strng( : 69 ) )
                end if
                write( xw_strng, 554 ) kt( 1 ), logt( 1 ),
     &                                 kt( n_logt ), logt( n_logt )
554             format( '      from kT=', f6.3, ' keV (logT=', f5.2,
     &                  ') to kT=', f6.3, ' keV (logT=', f5.2, ')' )
                call PWRITE( xw_strng( : 67 ) )
                write( xw_strng, 555 ) abun( 1 ), abun ( n_abu )
555             format( '      and abundances from ', f4.2,
     &                                                ' to ', f4.2 )
                call PWRITE( xw_strng( : 38 )  )
                warn_rs = .false.
              end if
              write( xw_strng, 561 ) kt( l ), logt( l )
561           format( '         Selected temperature is ', f6.3,
     &                                ' keV (log T is ', f5.2, ')' )
              call PWRITE( xw_strng( : 60 ) )
              write( xw_strng, 562 ) abun( k )
562           format( '         and selected abundance is ', f3.1 )
              call PWRITE( xw_strng( : 38 ) )
            end if

            file_st = p_root( :l_proot ) // ab_str( k )
     &                                           // '_' // logt_str( l )
            call LOCASE( file_st )
            call ARKOPN( mf_unit, mdir_name, file_st, 'mdl',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
            if( flag .eq. 0 ) then
*             Successfully opend the Raymond-Smith files
              if( m .eq. 1 ) highest = 0
              j_old( m ) = highest
              j = highest + 1
              beg_m( m ) = j
              do while( .true. )
                read( mf_unit, *, end = 590 ) e_in( j ), f_in( j )
                j = j + 1
                if( j .gt. m_pix ) then
                  call PWRITE( 'ERROR: Model file buffer full!!!' )
                  call PWRITE( '       (Increase m_pix in pimms.inc)' )
                  flag = -50
                  return
                end if
              end do
590           continue
              close( mf_unit )
              fin_m( m ) = j - 1
              highest = fin_m( m )
              e_old( m ) = 0.0
              outob( m ) = .false.
              if( np_num .le. 2 ) then
                one_line = ' '
                lflag = 0
                do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                  call GET_PSTR( 'Enter Nh >', one_line )
                  call RD_REAL( one_line, nh( m ), lflag )
                end do
              else
                call GT2_PARAM_N( np_beg + 2, nh( m ) )
              end if
              next = 3
              model_( m ) = 6 + p_model
              param( 1, m ) = kt( l )
              param( 2, m ) = logt( l )
              param( 3, m ) = abun( k )
            else
              call PWRITE( 'ERROR:: Failed to open the plasma file' )
              flag = -88
              return
            end if

          else
*           Was not a known model name --- try for a file
            call ARKOPN( mf_unit, ' ', model_st, 'mdl',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
            if( flag .ne. 0 ) then
              call ARKOPN( mf_unit, mdir_name, model_st, 'mdl',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
            end if
              
            if( flag .eq. 0 ) then
*             Was a file
              if( m .eq. 1 ) highest = 0
              j_old( m ) = highest
              j = highest + 1
              beg_m( m ) = j
              do while( .true. )
                read( mf_unit, *, end = 610 ) e_in( j ), f_in( j )
                j = j + 1
                if( j .gt. m_pix ) then
                  call PWRITE( 'ERROR: Model file buffer full!!!' )
                  call PWRITE( '       (Increase m_pix in pimms.inc)' )
                  flag = -50
                  return
                end if
              end do
610           continue
              close( mf_unit )
              fin_m( m ) = j - 1
              highest = fin_m( m )
              e_old( m ) = 0.0
              outob( m ) = .false.
              model_( m ) = -1
              mf_name( m ) = model_st
              if( m .eq. 1 ) then
                if( np_num .eq. 1 ) then
                  call GT2_PARAM_N( np_beg, nh( m ) )
                  next = 1
                else
                  nh( m ) = 0.0
                  next = 0
                end if
              else
                if( np_num .eq. 3 ) then
                  call GT2_PARAM_N( np_beg, nh( m ) )
                  next = 1
                else
                  nh( m ) = 0.0
                  next = 0
                end if
              end if
              param( 1, m ) = 0.0
              param( 2, m ) = 0.0
              param( 3, m ) = 0.0
            else
              call PWRITE( 'ERROR: Unknown model' )
              return
            end if
          end if

          if( nh( m ) .lt. 30.0 .and. nh( m ) .ne. 0.0 ) then
            call PWRITE( 'Value of Nh is rather small ---' )
            call PWRITE( 'Assuming that log10(Nh) was given ' )
            nh( m ) = 10.0 ** nh( m )
          end if

          if( m .gt. 1 ) then
*           Find normalization factors.
            if( np_num .le. next ) then
              if( model_( m ) .eq. 5 ) then
                one_line = ' '
                lflag = 0
                do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                  call GET_PSTR(
     &   'Enter equivalent width (eV) relative to comp. 1 >', one_line )
                  call RD_REAL( one_line, m_rat( m ), lflag )
                end do
                m_e( m ) = param( 1, m )
              else
                one_line = ' '
                lflag = 0
                do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                  call GET_PSTR(
     &                    'Enter flux relative to comp. 1 >', one_line )
                  call RD_REAL( one_line, m_rat( m ), lflag )
                end do
                one_line = ' '
                lflag = 0
                do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                  call GET_PSTR( '...at what energy ?', one_line )
                  call RD_REAL( one_line, m_e( m ), lflag )
                end do
              end if
            else if( np_num .eq. next + 1 ) then
              call GT2_PARAM_N( np_beg + next, m_rat( m ) )
              if( model_( m ) .eq. 5 ) then
                m_e( m ) = param( 1, m )
              else
                one_line = ' '
                lflag = 0
                do while( one_line .eq. ' ' .or. lflag .lt. 0 )
                  call GET_PSTR( '...at what energy ?', one_line )
                  call RD_REAL( one_line, m_e( m ), lflag )
                end do
              end if
            else
              call GT2_PARAM_N( np_beg + next, m_rat( m ) )
              if( model_( m ) .eq. 5 ) then
                call PWRITE(
     &              'Warning:: Ignoring extra numerical parameter(s)' )
                m_e( m ) = param( 1, m )
              else
                if( np_num .gt. next + 2 ) then
                  call PWRITE(
     &              'Warning:: Ignoring extra numerical parameter(s)' )
                end if
                call GT2_PARAM_N( np_beg + next + 1, m_e( m ) )
              end if
            end if
            n_comp = m
            norm( m ) = 1.0
            m_crrnt = 1
            val_1 = SPEC( m_e( m ) )
            m_crrnt = m
            val_m = SPEC( m_e( m ) )
            m_crrnt = 0
            if( model_( m ) .eq. 5 ) then
              if( val_1 .eq. 0.0 ) then
                n_comp = m - 1
                call PWRITE(
     &           'ERROR:: Cannot normalize components at this energy' )
                flag = -123
                return
              end if
              norm( m ) = val_1 * m_rat( m ) * 0.001
            else
              if( val_1 .eq. 0.0 .or. val_m .eq. 0.0 ) then
                n_comp = m - 1
                call PWRITE(
     &           'ERROR:: Cannot normalize components at this energy' )
                flag = -123
                return
              end if
              norm( m ) = ( val_1 * m_rat( m ) ) / val_m
            end if
          else
            if( np_num .gt. next ) then
              call PWRITE(
     &              'Warning:: Ignoring extra numerical parameter(s)' )
            end if
            norm( m ) = 1.0
          end if

          m = m + 1
          mm = mm + 1
        end do

        n_comp = nc_tgt - extra

        end



*+SEL_PLASM
        subroutine SEL_PLASM( flag )

        implicit none

        integer m_rs, m_abund
        real t_kev
        parameter( m_rs = 128, m_abund = 16 )
        parameter( t_kev = 1.16048e+07 )

        include 'pimms.inc'

*       Description:
*         Reads in the choice of plasma model from user via ARK-FCI
*
*       Arguments:
*         flag (output): status
*
*       Dependencies:
*         ARK-FCI, LENTRIM, LOCASE
*
*       Origin:
*         Created by KM
*
*       Author:
*         Koji Mukai, 2010 January
*         Minor re-ordering to conform to the strict standard, 2010 August
*-PMS_SLMDL

        integer flag

        integer l_im, j, k, l, m, save_fin, mf_unit
        logical changed, warn_rs, grid_warn
        integer n_abu, n_logt, l_proot
        real logt( m_rs ), kt( m_rs ), av_t
        real abun( m_abund )
        character*8 p_root
        character*3 logt_str( m_rs )
        character*2 ab_str( m_abund )
        character*32 file_st
        character*16 in_model
        character*80 xw_strng

        integer LENTRIM

        character*16 m_names( 3 )
        include 'sitespec.inc'
        data m_names / 'raymond-smith', 'apec', 'mekal' /

        call PARAM_C( 1, in_model, 'Enter plasma model name >' )
        call LOCASE( in_model )
        l_im = LENTRIM( in_model )
        flag = 0
        if( in_model .eq. 'rs'
     &          .or. in_model( :l_im ) .eq. m_names( 1 )( :l_im ) ) then
*         Raymond-Smith selected
          if( p_model .eq. 0 ) flag = 1
          p_model = 0
        else if( in_model( :l_im ) .eq. m_names( 2 )( :l_im ) ) then
*         APEC selected
          if( p_model .eq. 1 ) flag = 1
          p_model = 1
        else if( in_model( :l_im ) .eq. m_names( 3 )( :l_im ) ) then
*         mekal selected
          if( p_model .eq. 2 ) flag = 1
          p_model = 2
        else
*         Invalid plasma code name
          flag = -99
        end if
*       Flag=1 (non-error) means this plasma command was redundant ---
*              the user tried to re-choose the current plasma model
        changed = .false.
        if( flag .eq. 0 ) then
*         Check if the current model has a plasma component
          warn_rs = .true.
          do m = 1, n_comp
            if( model_( m ) .ge. 6 .and. model_( m ) .le. 8 ) then
*              YES
              changed = .true.
              if( p_model .eq. 0 ) then
                call ARKOPN( mf_unit, mdir_name, 'rs.idx', 'idx',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
              else if( p_model .eq. 1 ) then
                call ARKOPN( mf_unit, mdir_name, 'apec.idx', 'idx',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
              else if( p_model .eq. 2 ) then
                call ARKOPN( mf_unit, mdir_name, 'mekal.idx', 'idx',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
              end if
              if( flag .lt. 0 ) then
                call PWRITE( 'ERROR: Failed to opne plasma index file' )
                flag = -102
                return
              end if
              read( mf_unit, * ) p_root
              l_proot = LENTRIM( p_root )
              read( mf_unit, * ) n_abu
              do k = 1, n_abu
                read( mf_unit, * ) abun( k ), ab_str( k )
              end do
              read( mf_unit, * ) n_logt
              do l = 1, n_logt
                read( mf_unit, * ) logt( l ), logt_str( l )
                kt( l ) = 10.0 ** logt( l ) / t_kev
              end do
              close( mf_unit )

              grid_warn = .false.
              do l = 1, n_logt
                if( param( 2, m ) .eq. logt( l ) ) goto 530
              end do
              grid_warn = .true.

              l = 1
              do while( l .lt. n_logt )
                av_t = ( logt( l ) + logt( l + 1 ) ) * 0.5
                if( param( 2, m ) .lt. av_t ) goto 530
                l = l + 1
              end do
530           continue
*             The cases where temperature exactly matched rejoin here

              do k = 1, n_abu
                if( param( 3, m ) .eq. abun( k ) ) goto 540
              end do
              grid_warn = .true.
              k = 1
              do while( k .lt. n_abu )
                av_t = ( abun( k ) + abun( k + 1 ) ) * 0.5
                if( param( 3, m ) .lt. av_t ) goto 540
                k = k + 1
              end do
540           continue

              if( grid_warn ) then
                if( warn_rs ) then
                  if( p_model .eq. 0 ) then
                    write( xw_strng, 551 ) n_logt, n_abu
551                 format( 'NOTE: This version of PIMMS has a grid ',
     &             'of ', i3, 'x', i2, ' grid of Raymond-Smith models' )
                    call PWRITE( xw_strng( : 77 ) )
                  else if( p_model .eq. 1 ) then
                    write( xw_strng, 552 ) n_logt, n_abu
552                 format( 'NOTE: This version of PIMMS has a grid ',
     &                     'of ', i3, 'x', i2, ' grid of APEC models' )
                    call PWRITE( xw_strng( : 68 ) )
                  else if( p_model .eq. 2 ) then
                    write( xw_strng, 553 ) n_logt, n_abu
553                 format( 'NOTE: This version of PIMMS has a grid ',
     &                    'of ', i3, 'x', i2, ' grid of mekal models' )
                    call PWRITE( xw_strng( : 69 ) )
                  end if
                  write( xw_strng, 554 ) kt( 1 ), logt( 1 ),
     &                                 kt( n_logt ), logt( n_logt )
554               format( '      from kT=', f6.3, ' keV (logT=', f5.2,
     &                  ') to kT=', f6.3, ' keV (logT=', f5.2, ')' )
                  call PWRITE( xw_strng( : 67 ) )
                  write( xw_strng, 555 ) abun( 1 ), abun ( n_abu )
555               format( '      and abundances from ', f4.2,
     &                                                ' to ', f4.2 )
                  call PWRITE( xw_strng( : 38 )  )
                  warn_rs = .false.
                end if
                write( xw_strng, 561 ) kt( l ), logt( l )
561             format( '         Selected temperature is ', f6.3,
     &                                ' keV (log T is ', f5.2, ')' )
                call PWRITE( xw_strng( : 60 ) )
                write( xw_strng, 562 ) abun( k )
562             format( '         and selected abundance is ', f3.1 )
                call PWRITE( xw_strng( : 38 ) )
              end if

              file_st = p_root( :l_proot ) // ab_str( k )
     &                                           // '_' // logt_str( l )
              call LOCASE( file_st )
              call ARKOPN( mf_unit, mdir_name, file_st, 'mdl',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
              if( flag .eq. 0 ) then
*               Successfully opend the Plasma file
                if( m .eq. 1 ) highest = 0
                j_old( m ) = highest
                j = highest + 1
                save_fin = fin_m( m )
                beg_m( m ) = j
                do while( .true. )
                  read( mf_unit, *, end = 590 ) e_in( j ), f_in( j )
                  j = j + 1
                  if( j .gt. m_pix ) then
                    call PWRITE( 'ERROR: Model file buffer full!!!' )
                    call PWRITE( '      (Increase m_pix in pimms.inc)' )
                    flag = -50
                    return
                  end if
                end do
590             continue
                close( mf_unit )
                fin_m( m ) = j - 1
                if( save_fin .ne. fin_m( m ) ) then
*                 Replacement plsama model has a different # of points
                  call PWRITE(
     &       'ERROR: Mismatch in number of points among plasma models' )
                  flag = -60
                  return
                end if
                highest = fin_m( m )
                e_old( m ) = 0.0
                outob( m ) = .false.
                model_( m ) = 6 + p_model
              else
                call PWRITE( 'ERROR:: Failed to open the plasma file' )
                flag = -88
                return
              end if
            end if
          end do
          if( changed ) then
            call PWRITE( 'Current model has been changed' )
            call PMS_WTMDL( 1 )
          else
            if( p_model .eq. 0 ) then
              call PWRITE( 'Plasma code changed to Raymond-Smith' )
            else if( p_model .eq. 1 ) then
              call PWRITE( 'Plasma code changed to APEC' )
            else if( p_model .eq. 2 ) then
              call PWRITE( 'Plasma code changed to mekal' )
            end if
            call PWRITE(
     &        'This will take effect next time a plasma model is used' )
          end if
        end if

        end

