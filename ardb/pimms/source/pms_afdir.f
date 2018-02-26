*+PMS_AFDIR
        subroutine PMS_AFDIR( )

        implicit none

*       Description:
*         Shows the directory of missions
*
*       Arguments:
*         None
*
*       Dependencies:
*         PIMMS Index structure, ARK/FCI
*
*       Origin:
*         Created by KM
*
*       Author:
*         Koji Mukai, 1993 March, original version
*-PMS_AFDIR

        integer num_c, num_f, num_i
        character*128 mssn_st, dtct_st
        integer mssn_in, dtct_in

        integer PDX_GTMIN, PDX_GTDIN

        call INQ_PARAM( num_c, num_f, num_i )
        if( num_c .eq. 0 ) then
*         DIRECTORY called without arguments, do a full listing
          call PDX_SHLST( 0, 0 )
        else
*         At least mission name is given
          call GET_PARAM_C( 1, mssn_st )
          call LOCASE( mssn_st )
          mssn_in = PDX_GTMIN( mssn_st )
          if( mssn_in .le. 0 ) then
*           Some mistake in mission name, do a full listing
            call PWRITE(
     &           'Mission name unknown, giving full directory listing' )
            call PDX_SHLST( 0, 0 )
          else
*           Mission name understood
            if( num_c .eq. 2 ) then
*             Detector name also given
              call GET_PARAM_C( 2, dtct_st )
              call LOCASE( dtct_st )
              dtct_in = PDX_GTDIN( mssn_in, dtct_st )
              if( dtct_in .le. 0 ) then
*               Some mistake in detector name, do a full listing for mission
                call PWRITE( 'Detector name unknown' )
                call PDX_SHLST( mssn_in, 0 )
              else
*               Mission and detector names understood
                call PDX_SHLST( mssn_in, dtct_in )
              end if
            else
*             Only mission name given, list all detector/filter combinations
              call PDX_SHLST( mssn_in, 0 )
            end if
          end if
        end if

        end
