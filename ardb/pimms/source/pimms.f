        program PIMMS

*
*       Portable, Interactive, Multi-Mission Simulator
*

        implicit none

        include 'pimms.inc'

	integer command_no, status

	integer f_mssn, f_dtct, f_unit, f_fltr
	integer t_mssn, t_dtct, t_unit, t_fltr
	integer n_mssn, n_dtct, n_unit, n_fltr
        integer num_c, num_f, num_i
        integer b_unit, flag
        real f_lo, f_hi, t_lo, t_hi, n_lo, n_hi, input
        character*128 help_str, lf_name

        character*16, commands( 12 ), prompt

        integer LENTRIM

        include 'sitespec.inc'

        data commands / 'PLASMA', 'MODEL', 'FROM', 'INSTRUMENT', 'TO',
     &                  'GO', 'SHOW', 'DIRECTORY', 'LOG', 'OUTPUT',
     &                  'QUIT', 'EXIT' /
        data prompt / 'PIMMS >' /

        call INTCOM( commands, 12, prompt )

        call ARKOPN( b_unit, ddir_name, 'pms_banner.txt', 'txt',
     &               'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &               1, flag )
        if( flag .ne. 0 )  then
          call PWRITE(
     &             'FATAL ERROR:: failed to open PIMMS banner file' )
          call PWRITE( 'PIMMS was looking for ''pms_banner.txt'' in' )
          call PWRITE( ddir_name( : LENTRIM( ddir_name ) ) )
          stop
        end if

        do while( .true. )
          read( b_unit, '(a)', end=90 ) help_str
          call PWRITE( help_str( : LENTRIM( help_str ) ) )
        end do
 90     continue
        close( b_unit )

*       Initialize the PIMMS Index structure
        call PMS_INTLZ( )

*     HARDWIRED section begins
*       initial model is bremsstrahlung, 10 keV, nh=1E21
*       also, set p_model to 1, which means the plasma model is APEC
        n_comp = 1
        m_crrnt = 0
        p_model = 1
        norm( 1 ) = 1.0
        model_( 1 ) = 2
        param( 1, 1 ) = 10.0
        nh( 1 ) = 1.0e+21
        z1 = 1.0
        nh_g = 0.0

*       initially from 2-10 keV flux in ergs
        f_mssn = 0
        f_unit = 1
        f_lo = 2.0
        f_hi = 10.0
        t_lo = 0.0
        t_hi = 0.0
*     HARDWIRED section ends

*       initially to whatever mission is first on the pms_mssn.lst file
        call PDX_FIRST( t_mssn, t_dtct, t_fltr )

        call FORCE_COMM( 'SHOW' )
        
*       BEGIN infinite loop
        do while( .true. )

*         Get command from FCI
          call GETCOM( command_no )

          status = 0
********************************************************************************
          if( command_no .eq. 11 .or. command_no .eq. 12
     &                .or. command_no .eq. 13 ) then
*                        QUIT command
            goto 1000
********************************************************************************
          else if( command_no .eq. 1 ) then
*                        PLASMA command
            call SEL_PLASM( status )
            if( status .lt. 0 ) then
              call CLRCOM( )
              call PWRITE( 'ERROR:: Failed to choose the plasma model' )
            end if
********************************************************************************
          else if( command_no .eq. 2 ) then
*                        MODEL command
            call PMS_SLMDL( status )
            if( status .lt. 0 ) then
              call CLRCOM( )
              call PWRITE( 'ERROR selecting spectral model' )
            end if
********************************************************************************
          else if( command_no .eq. 3 ) then
c                        ! begin FROM command
            call PMS_SELCT
     &         ( 1, f_mssn, f_dtct, f_fltr, f_unit, f_lo, f_hi, status )
            if( status .lt. 0 ) then
              call CLRCOM( )
              call PWRITE( 'ERROR selecting output mission' )
            end if
********************************************************************************
          else if( command_no .eq. 4 .or. command_no .eq. 5 ) then
*                        INSTRUMENT command
            call PMS_SELCT
     &         ( 1, t_mssn, t_dtct, t_fltr, t_unit, t_lo, t_hi, status )
            if( status .lt. 0 ) then
              call CLRCOM( )
              call PWRITE( 'ERROR selecting output mission' )
            end if
********************************************************************************
          else if( command_no .eq. 6 ) then
*                        GO command
            call PARAM_N( 1, input, 1.0, 'Enter input rate >' )
            call INQ_PARAM( num_c, num_f, num_i )
            if( num_c .ge. 1 ) then
              call PMS_SELCT
     &         ( 2, n_mssn, n_dtct, n_fltr, n_unit, n_lo, n_hi, status )
              if( status .eq. 0 ) then
                call PMS_DOCRT( input, n_mssn, n_dtct, n_fltr, n_unit,
     &          n_lo, n_hi, t_mssn, t_dtct, t_fltr, t_unit, t_lo, t_hi )
              end if
            else
*               Use default mission for converting FROM
              call PMS_DOCRT( input, f_mssn, f_dtct, f_fltr, f_unit,
     &          f_lo, f_hi, t_mssn, t_dtct, t_fltr, t_unit, t_lo, t_hi )
            end if
********************************************************************************
          else if( command_no .eq. 7 ) then
*                        SHOW command
            call PMS_WTMDL( 1 )
            call PWRITE(
     &                   '   <--- Use ''MODEL'' command to change' )
            call PWRITE(
     &  '        and ''PLASMA'' command to switch among APEC/mekal/RS' )
            call PWRITE(
     &                       '* By default, input rate is taken to be' )
            call PMS_WTSLC( f_mssn, f_dtct, f_fltr, f_unit, f_lo, f_hi )
            call PWRITE(
     &            '   <--- Use ''FROM'' command to change the default' )
            call PWRITE( '* Simulation product will be' )
            call PMS_WTSLC( t_mssn, t_dtct, t_fltr, t_unit, t_lo, t_hi )
            call PWRITE( '   <--- Use ''INSTRUMENT'' command '
     &                            // 'to switch to another instrument' )
********************************************************************************
          else if( command_no .eq. 8 ) then
*                        DIRECTORY command
            call PMS_AFDIR( )
********************************************************************************
          else if( command_no .eq. 9 ) then
*                        LOG command
            call PARAM_C( 1, lf_name,
     &                           'Enter log file name (or ''close'') ' )
            call LOCASE( lf_name )
            if( lf_name .eq. 'close' ) then
              call PWCLOS( )
            else
              call PWOPEN( lf_name )
            end if
********************************************************************************
          else if( command_no .eq. 10 ) then
*                        OUTPUT command
            call PMS_OUTPT( status )
            if( status .lt. 0 ) then
              call CLRCOM( )
              call PWRITE( 'ERROR:: Failed to write the current model' )
            end if
********************************************************************************
          else if( command_no .eq. 14 ) then
*                        HELP command
            call INQ_PARAM( num_c, num_f, num_i )
            if( num_c .eq. 1 ) then
              call GET_PARAM_C( 1, help_str )
            else
              help_str = ' '
            end if
            call ARK_HELP( ddir_name, 'pimms.ahl', help_str )
********************************************************************************
          else
            call PWRITE( '* PIMMS commands are: ' )
            call PWRITE( '  PLASMA, MODEL, FROM, INSTRUMENT, SHOW, '
     &                      // 'DIRECTORY, LOG, OUTPUT, HELP and QUIT' )
          end if
c                        ! end of command name if block

        end do
c                        ! end of infinite loop

1000    continue

        end
