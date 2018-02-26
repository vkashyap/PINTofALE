*+PDX_SETDT
        subroutine PDX_SETDT( mission, detector, entry )

        implicit none

        include 'pms_index.inc'

        integer mission, detector
        integer entry

*       Description:
*         Creates index entry for mission, detector
*
*       Arguments:
*         mission      (i) : Mission number
*         detector     (i) : Detector number
*         entry        (i) : Index number for the above combination
*
*       Dependencies:
*         None
*
*       Origin:
*         Cooked up by KM to handle the PIMMS-specific database.
*
*       Author:
*         Koji Mukai, 1993 Feb 26, Original version
*
*         Koji Mukai, 2010 Jan 15, bug fix for mission with no choice
*                                  of detectors or filters
*-PDX_SETDT

        integer offset

        offset = ind_mssn( mission )
        ind_dtct( offset + detector ) = entry

        end



*+PDX_GETDT
        integer function PDX_GETDT( mission, detector )

        implicit none

        include 'pms_index.inc'

        integer mission, detector

*       Description:
*         Retrieve index entry for mission, detector
*
*       Arguments:
*         mission      (i) : Mission number
*         detector     (i) : Detector number
*         <PDX_GETDT>  (r) : Index number for the above combination
*
*       Dependencies:
*         None
*
*       Origin:
*         Cooked up by KM to handle the PIMMS-specific database.
*
*       Author:
*         Koji Mukai, 1993 Feb 26, Original version
*-PDX_GETDT

        integer offset, n_dtct

        if( mission .ge. 1 .and. mission .le. n_mssn ) then
          offset = ind_mssn( mission )
          n_dtct = st_array( offset )

          if( n_dtct .gt. 0 .and.
     &        ( detector .ge. 1 .and. detector .le. n_dtct ) ) then
            PDX_GETDT = ind_dtct( offset + detector )
          else if( n_dtct .lt. 0 .and. detector .eq. 0 ) then
            PDX_GETDT = 0
          else
            PDX_GETDT = -2
          end if
        else
          PDX_GETDT = -1
        end if

        end



*+PDX_SETFL
        subroutine PDX_SETFL( mission, detector, filter, entry )

        implicit none

        include 'pms_index.inc'

        integer mission, detector, filter
        integer entry

*       Description:
*         Creates index entry for mission, detector, filter
*
*       Arguments:
*         mission      (i) : Mission number
*         detector     (i) : Detector number
*         filter       (i) : Filter number
*         entry        (i) : Index number for the above combination
*
*       Dependencies:
*         None
*
*       Origin:
*         Cooked up by KM to handle the PIMMS-specific database.
*
*       Author:
*         Koji Mukai, 1993 Feb 26, Original version
*-PDX_SETFL

        integer offset
        integer PDX_GETDT

        offset = PDX_GETDT( mission, detector )
        ind_fltr( offset + filter ) = entry

        end



*+PDX_GETFL
        integer function PDX_GETFL( mission, detector, filter )

        implicit none

        include 'pms_index.inc'

        integer mission, detector, filter

*       Description:
*         Retrieve index entry for mission, detector, filter
*
*       Arguments:
*         mission      (i) : Mission number
*         detector     (i) : Detector number
*         filter       (i) : Filter number
*         <PDX_GETFL>  (r) : Index number for the above combination
*
*       Dependencies:
*         None
*
*       Origin:
*         Cooked up by KM to handle the PIMMS-specific database.
*
*       Author:
*         Koji Mukai, 1993 Feb 26, Original version
*-PDX_GETFL

        integer offset, n_fltr
        integer PDX_GETDT

        offset = PDX_GETDT( mission, detector )
        if( offset .lt. 0 ) then
          PDX_GETFL = offset
        else if( offset .eq. 0 ) then
          if( filter .eq. 0 ) then
            PDX_GETFL = 0
          else
            PDX_GETFL = -3
          end if
        else
          n_fltr = st_array( offset )
          if( n_fltr .le. 0 ) then
            if( filter .eq. 0 ) then
              PDX_GETFL = 0
            else
              PDX_GETFL = -3
            end if
          else if( filter .ge. 1 .and. filter .le. n_fltr ) then
            PDX_GETFL = ind_fltr( offset + filter )
          else
            PDX_GETFL = -3
          end if
        end if

        end



*+PDX_CKALL
        subroutine PDX_CKALL( mission, detector, filter, status )

        implicit none

        integer mission, detector, filter
        integer status

*       Description
*         Checks for effective area file (.area) etc. and returns
*         the encoded info.
*
*       Arguments
*         mission      (i) : Input mission number
*         detector     (i) : Input detector number
*         filter       (i) : Input filter number
*         status       (o) : Return status
*
*       Dependencies:
*         PDX_MKNAM and other PIMMS index routines
*
*       Origin:
*         Created by KM to combat boredom
*
*       Author:
*         Koji Mukai, 1993 Mar 08, the original version
*-PDX_CKALL

        character*256 temp_name
        integer temp_unit
        integer flag

        include 'sitespec.inc'

        status = 0
        call PDX_MKNAM( mission, detector, filter, 'AREA',
     &                                            temp_name, flag )
        if( flag .ge. 0 ) then
          call ARKOPN( temp_unit, ddir_name, temp_name, 'area',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
        end if
        if( flag .eq. 0 ) then
          status = status + 1
          close( temp_unit )
        end if

        call PDX_MKNAM( mission, detector, filter, 'SPECIAL',
     &                                            temp_name, flag )
        if( flag .ge. 0 ) then
          call ARKOPN( temp_unit, ddir_name, temp_name, 'SPECIAL',
     &                 'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                 1, flag )
        end if
        if( flag .eq. 0 ) then
          status = status + 100
          close( temp_unit )
        end if

        end



*+PDX_CKMOD
        logical function PDX_CKMOD( record )

        implicit none

        integer record

*       Description:
*         Sanity check
*
*       Arguments:
*         record   (i) : Index record
*
*       Dependencies:
*         None
*
*       Origin:
*         Created by KM as part of the PIMMS Index routines
*
*       Author:
*         Koji Mukai (1993 Mar 22) Original Version
*         KM, 1994 Dec, small PIMMS version
*-PDX_CKMOD

        if( mod( abs( record ), 10 ) .eq. 1 ) then
          PDX_CKMOD = .true.
        else
          PDX_CKMOD = .false.
        end if

        end

*+PDX_MKNAM
        subroutine PDX_MKNAM( mission, detector, filter, ext_name,
     &                                                f_name, flag )

        implicit none

        include 'pms_index.inc'

        integer mission, detector, filter
        character*( * ) ext_name
        character*( * ) f_name
        integer flag

*       Description:
*         Concatnates disk and directory names, mission, detector and
*         filter numbers and returns a XANLIB-compatible complete
*         file name specification
*
*       Arguments:
*         mission      (i) : Mission number, as stored in index
*         detector     (i) : Detector number, as stored in index
*         filter       (i) : Filter number, as stored in index
*         ext_name     (i) : Extension name.
*         f_name       (o) : Completed file name
*         flag         (o) : Negative if error occurred
*
*       Dependencies:
*         None?
*
*       Origin:
*         Created by KM for managing files
*
*       Author:
*         Koji Mukai, 1993 Mar 01, Original Version
*-PDX_MKNAM

        integer len_out, len_seg, last, index
        integer LENTRIM, PDX_GETDT, PDX_GETFL

        len_out = len( f_name )
        f_name = ' '
        last = 0

        if( mission .le. 0 .or. mission .gt. n_mssn ) then
          flag = -11
          goto 900
        end if
        index = ind_mssn( mission )
        len_seg = ln_array( index ) + 1
        if( last + len_seg .gt. len_out ) then
          flag = -1
          goto 900
        end if
        f_name( last + 1: last + len_seg )
     &                = nm_array( index )( : len_seg - 1 ) // '_'
        last = last + len_seg

        index = PDX_GETDT( mission, detector )
        if( index .lt. 0 ) then
          flag = -10 + index
          goto 900
        else if( index .eq. 0 ) then
          f_name( last + 1: last + 1 ) = '_'
          last = last + 1
        else if( index .gt. 0 ) then
          len_seg = ln_array( index ) + 1
          if( last + len_seg .gt. len_out ) then
            flag = -1
            goto 900
          end if
          f_name( last + 1: last + len_seg )
     &                = nm_array( index )( : len_seg - 1 ) // '_'
          last = last + len_seg

          index = PDX_GETFL( mission, detector, filter )
          if( index .lt. 0 ) then
            flag = -10 + index
            goto 900
          else if( index .gt. 0 ) then
            len_seg = ln_array( index )
            if( last + len_seg .gt. len_out ) then
              flag = -1
              goto 900
            end if
            f_name( last + 1: last + len_seg )
     &                = nm_array( index )( : len_seg )
            last = last + len_seg
          end if
        end if

        len_seg = LENTRIM( ext_name ) + 1
        if( last + len_seg .gt. len_out ) then
          flag = -1
          goto 900
        end if
        f_name( last + 1: last + len_seg )
     &                               = '.' // ext_name( : len_seg )
        call LOCASE( f_name( : last + len_seg ) )

        flag = 0
        return

900     continue

        end


*+PDX_STRNG
        subroutine PDX_STRNG( mission, detector, filter, string, s_len )

        implicit none

        include 'pms_index.inc'

        integer mission, detector, filter
        character*( * ) string
        integer s_len

*       Description:
*         Concatnates mission, detector and filter names.
*
*       Arguments:
*         mission      (i) : Mission number, as stored in index
*         detector     (i) : Detector number, as stored in index
*         filter       (i) : Filter number, as stored in index
*         string       (o) : Output of the form 'EXOSAT LE ALP'
*         s_len        (o) : Effective length of the string,
*                            or negative if error occured.
*
*       Dependencies:
*         None?
*
*       Origin:
*         Created by KM for managing files
*
*       Author:
*         Koji Mukai, 1993 Mar 01, Original Version
*-PDX_STRNG

        integer len_out, len_seg, last, index
        integer PDX_GETDT, PDX_GETFL

        len_out = len( string )
        string = ' '
        last = 0

        if( mission .le. 0 .or. mission .gt. n_mssn ) then
          s_len = -11
          goto 900
        end if
        index = ind_mssn( mission )
        len_seg = ln_array( index ) + 1
        if( last + len_seg .gt. len_out ) then
          s_len = -1
          goto 900
        end if
        string( last + 1: last + len_seg )
     &                = nm_array( index )( : len_seg - 1 ) // ' '
        last = last + len_seg

        index = PDX_GETDT( mission, detector )
        if( index .lt. 0 ) then
          s_len = -10 + index
          goto 900
        else if( index .gt. 0 ) then
          len_seg = ln_array( index ) + 1
          if( last + len_seg .gt. len_out ) then
            s_len = -1
            goto 900
          end if
          string( last + 1: last + len_seg )
     &                = nm_array( index )( : len_seg - 1 ) // ' '
          last = last + len_seg

          index = PDX_GETFL( mission, detector, filter )
          if( index .lt. 0 ) then
            s_len = -10 + index
            goto 900
          else if( index .gt. 0 ) then
            len_seg = ln_array( index )
            if( last + len_seg .gt. len_out ) then
              s_len = -1
              goto 900
            end if
            string( last + 1: last + len_seg )
     &                = nm_array( index )( : len_seg )
            last = last + len_seg
          end if
        end if

        s_len = last
        return

900     continue

        end


*+PDX_SHLST
        subroutine PDX_SHLST( in_mission, in_detector )

        implicit none

        include 'pms_index.inc'

        integer in_mission, in_detector

*       Description:
*         Writes out a list of missions, detector, filters and what can
*         be done.
*
*       Arguments:
*         in_mission   (i) : input mission number
*         in_detector  (i) : input detector number
*
*       Dependencies:
*         Other PIMMS Index routines, and LENTRIM, AH_MORE
*
*       Origin:
*         Created by KM as part of PIMMS Index routines
*
*       Author:
*         Koji Mukai (1993 Mar 19) Original version
*-PDX_SHLST

        character*( max_len ) the_name
        character*80 xw_strng, dummy
        integer name_len, status, n_dtct, n_fltr
        integer mission, detector, filter, j, k, l, n_char, l_page

        integer PDX_GETDT, PDX_GETFL, LENTRIM

        l_page = 0

        if( in_mission .eq. 0 ) then
          do j = 1, n_mssn
            mission = ind_mssn( j )
            the_name = nm_array( mission )
            call UPCASE( the_name )
            name_len = ln_array( mission )
            status = st_array( mission )
            if( status .lt. 0 ) then
*             Dead end --- just list status here
              call PWRITE( 'Mission ' // the_name( : name_len ) )
              call AH_MORE( l_page, dummy )
              call PDX_WTSTS( -status, 1, l_page )
            else
              write( xw_strng,
     &               '(''Mission '',a,'' with'',i2,''  detectors'')' )
     &                                    the_name( : name_len ), status
              n_char = LENTRIM( xw_strng )
              call PWRITE( xw_strng( : n_char ) )
              call AH_MORE( l_page, dummy )
              n_dtct = status
              do k = 1, n_dtct
                detector = PDX_GETDT( j, k )
                the_name = nm_array( detector )
                call UPCASE( the_name )
                name_len = ln_array( detector )
                status = st_array( detector )
                if( status .lt. 0 ) then
*                 Dead end --- just list status here
                  call PWRITE( '  Detector '
     &                                // the_name( : name_len ) )
                  call AH_MORE( l_page, dummy )
                  call PDX_WTSTS( -status, 2, l_page )
                else
                  write( xw_strng,
     &         '(''  Detector '',a,'' with '',i2,'' filters'')' )
     &                                    the_name( : name_len ), status
                  n_char = LENTRIM( xw_strng )
                  call PWRITE( xw_strng( : n_char ) )
                  call AH_MORE( l_page, dummy )
                  n_fltr = status
                  do l = 1, n_fltr
                    filter = PDX_GETFL( j, k, l )
                    the_name = nm_array( filter )
                    call UPCASE( the_name )
                    name_len = ln_array( filter )
                    call PWRITE( '    Filter '
     &                                // the_name( : name_len ) )
                    call AH_MORE( l_page, dummy )
                    status = st_array( filter )
                    call PDX_WTSTS( status, 3, l_page )
                  end do
                end if
              end do
            end if
          end do
        else if( in_mission .lt. 0 .or. in_mission .gt. n_mssn ) then
          call PWRITE( 'ERROR in PDX_SHLST:: PIMMS is confused' )
        else
          mission = ind_mssn( in_mission )
          the_name = nm_array( mission )
          call UPCASE( the_name )
          name_len = ln_array( mission )
          call PWRITE( 'Mission ' // the_name( : name_len ) )
          call AH_MORE( l_page, dummy )
          status = st_array( mission )
          if( in_detector .eq. 0 ) then
            if( status .lt. 0 ) then
*             Dead end --- just list status here
              call PDX_WTSTS( -status, 1, l_page )
            else
              write( xw_strng,
     &               '(''       has '', I2, ''  detectors.'')' ) status
              call PWRITE( xw_strng( : 26 ) )
              call AH_MORE( l_page, dummy )
              n_dtct = status
              do k = 1, n_dtct
                detector = PDX_GETDT( in_mission, k )
                the_name = nm_array( detector )
                call UPCASE( the_name )
                name_len = ln_array( detector )
                call PWRITE( '  Detector '
     &                                // the_name( : name_len ) )
                call AH_MORE( l_page, dummy )
                status = st_array( detector )
                if( status .lt. 0 ) then
*                 Dead end --- just list status here
                  call PDX_WTSTS( -status, 2, l_page )
                else
                  write( xw_strng,
     &               '(''         has '', I2, ''  filters.'')' ) status
                  call PWRITE( xw_strng( : 26 ) )
                  call AH_MORE( l_page, dummy )
                  n_fltr = status
                  do l = 1, n_fltr
                    filter = PDX_GETFL( in_mission, k, l )
                    the_name = nm_array( filter )
                    call UPCASE( the_name )
                    name_len = ln_array( filter )
                    call PWRITE( '    Filter '
     &                                // the_name( : name_len ) )
                    call AH_MORE( l_page, dummy )
                    status = st_array( filter )
                    call PDX_WTSTS( status, 3, l_page )
                  end do
                end if
              end do
            end if
          else if( in_detector .lt. 0
     &                            .or. in_detector .gt. status ) then
            call PWRITE( 'ERROR in PDX_SHLST:: PIMMS is confused' )
          else
            detector = PDX_GETDT( in_mission, in_detector )
            the_name = nm_array( detector )
            call UPCASE( the_name )
            name_len = ln_array( detector )
            call PWRITE( '  Detector '
     &                                // the_name( : name_len ) )
            call AH_MORE( l_page, dummy )
            status = st_array( detector )
            if( status .lt. 0 ) then
*             Dead end --- just list status here
              call PDX_WTSTS( -status, 2, l_page )
            else
              write( xw_strng,
     &               '(''         has '', I2, ''  filters.'')' ) status
              call PWRITE( xw_strng( : 26 ) )
              call AH_MORE( l_page, dummy )
              n_fltr = status
              do l = 1, n_fltr
                filter = PDX_GETFL( in_mission, in_detector, l )
                the_name = nm_array( filter )
                call UPCASE( the_name )
                name_len = ln_array( filter )
                call PWRITE( '    Filter '
     &                                // the_name( : name_len ) )
                call AH_MORE( l_page, dummy )
                status = st_array( filter )
                call PDX_WTSTS( status, 3, l_page )
              end do
            end if
          end if
        end if

        end

*+PDX_WTSTS
        subroutine PDX_WTSTS( status, space, l_page )

        implicit none

        integer status, space, l_page

*       Description:
*         Given the status code, PWRITEs a suitable comment
*
*       Arguments:
*         status       (i) : mission/detector/filter status
*         space        (i) : precedes the message proper with varying
*                               number of spaces
*         l_page      (i/o): Necessary for AH_MORE
*
*       Dependencies:
*         PWRITE, AH_MORE
*
*       Origin:
*         Created as a part of the PIMMS Index structure
*
*       Author:
*         Koji Mukai (1993 Mar 19) original version
*-PDX_WTSTS

        character*80 dummy

        character space_1*4, space_2*6, space_3*8
        data space_1, space_2, space_3 / ' ', ' ', ' ' /

*        if( mod( status, 10 ) .eq. 1 ) then
*          if( space .eq. 1 ) then
*            call PWRITE( space_1 //
*     &                   'Count rate information available' )
*          else if( space .eq. 2 ) then
*            call PWRITE( space_2 //
*     &                   'Count rate information available' )
*          else if( space .eq. 3 ) then
*            call PWRITE( space_3 //
*     &                   'Count rate information available' )
*          end if
*        end if

        if( mod( status / 10, 10 ) .eq. 1 ) then
          if( space .eq. 1 ) then
            call PWRITE( space_1 //
     &                   'Images can be simulated' )
            call AH_MORE( l_page, dummy )
          else if( space .eq. 2 ) then
            call PWRITE( space_2 //
     &                   'Images can be simulated' )
            call AH_MORE( l_page, dummy )
          else if( space .eq. 3 ) then
            call PWRITE( space_3 //
     &                   'Images can be simulated' )
            call AH_MORE( l_page, dummy )
          end if
        end if

        if( mod( status / 100, 10 ) .eq. 1 ) then
          if( space .eq. 1 ) then
            call PWRITE( space_1 //
     &              'Instrument specific information is available' )
            call AH_MORE( l_page, dummy )
          else if( space .eq. 2 ) then
            call PWRITE( space_2 //
     &              'Instrument specific information is available' )
            call AH_MORE( l_page, dummy )
          else if( space .eq. 3 ) then
            call PWRITE( space_3 //
     &              'Instrument specific information is available' )
            call AH_MORE( l_page, dummy )
          end if
        end if
        
        end

*+PDX_GTMIN
        integer function PDX_GTMIN( ms_name )

        implicit none

        include 'pms_index.inc'

        character*( * ) ms_name

*       Description:
*         Given the mission name, returns the mission number
*
*       Arguments:
*         ms_name       (i) : Mission name string
*         <PDX_GTMIN>   (r) : 0 if given name string does not match any
*                               mission names in index,
*                             negative if given name string is ambiguous
*                               (matches -PDX_GTMIN possibilities)
*                             positive if matches one (mission number is
*                               returned in this case)
*
*       Dependencies:
*         LENTRIM
*
*       Origin:
*         Created as a part of the PIMMS Index structure
*
*       Author:
*         Koji Mukai (1993 Mar 19) original version
*-PDX_GTMIN

        integer msa_len, ms_len, j, pointer, exact, matches
        integer match_index( max_mssn )
        integer LENTRIM

        ms_len = LENTRIM( ms_name )
*       Assume it has been converted to lower case alerady, now we've
*       found the actual length of the character.

        matches = 0
        exact = 0
        do j = 1, n_mssn
          pointer = ind_mssn( j )
          if( ms_name( : ms_len )
     &                       .eq. nm_array( pointer )( : ms_len ) ) then
            matches = matches + 1
            match_index( matches ) = j
            msa_len = LENTRIM( nm_array( pointer ) )
            if( ms_len .eq. msa_len ) exact = matches
          end if
        end do

        if( matches .eq. 0 ) then
          PDX_GTMIN = 0
        else if( matches .eq. 1 ) then
          PDX_GTMIN = match_index( matches )
        else
          if( exact .eq. 0 ) then
            PDX_GTMIN = -matches
          else
            PDX_GTMIN = match_index( exact )
          end if
        end if

        end


*+PDX_GTDIN
        integer function PDX_GTDIN( mission, dt_name )

        implicit none

        include 'pms_index.inc'

        integer mission
        character*( * ) dt_name

*       Description:
*         Given the mission number and detector name, returns the detector
*         number.
*
*       Arguments:
*         mission       (i) : Mission number
*         dt_name       (i) : Detector name string
*         <PDX_GTDIN>   (r) : 0 if given name string does not match any
*                               detector names in index,
*                             -1 if error in mission number/no detectors
*                             otherwise negative if ambiguous
*                               (matches -PDX_GTDIN possibilities)
*                             positive if matches one (detector number is
*                               returned in this case)
*
*       Dependencies:
*         LENTRIM and PDX_GETDT
*
*       Origin:
*         Created as a part of the PIMMS Index structure
*
*       Author:
*         Koji Mukai (1993 Mar 20) original version
*-PDX_GTDIN

        integer dta_len, dt_len, j, pointer, exact, matches, n_dtct
        integer match_index( max_mssn )
        integer LENTRIM, PDX_GETDT

        if( mission .le. 0 .or. mission .gt. n_mssn ) then
          PDX_GTDIN = -1
        else
          pointer = ind_mssn( mission )
          n_dtct = st_array( pointer )
          if( n_dtct .lt. 0 ) then
            PDX_GTDIN = -1
          else

            dt_len = LENTRIM( dt_name )
*           Assume it has been converted to lower case alerady, now we've
*           found the actual length of the character.

            matches = 0
            exact = 0
            do j = 1, n_dtct
              pointer = PDX_GETDT( mission, j )
              if( dt_name( : dt_len )
     &                       .eq. nm_array( pointer )( : dt_len ) ) then
                matches = matches + 1
                match_index( matches ) = j
                dta_len = LENTRIM( nm_array( pointer ) )
                if( dta_len .eq. dt_len ) exact = matches
              end if
            end do

            if( matches .eq. 0 ) then
              PDX_GTDIN = 0
            else if( matches .eq. 1 ) then
              PDX_GTDIN = match_index( matches )
            else
              if( exact .eq. 0 ) then
                PDX_GTDIN = -matches
              else
                PDX_GTDIN = match_index( exact )
              end if
            end if

          end if
        end if

        end


*+PDX_GTFIN
        integer function PDX_GTFIN( mission, detector, fl_name )

        implicit none

        include 'pms_index.inc'

        integer mission, detector
        character*( * ) fl_name

*       Description:
*         Given the mission and detector numbers and filter name, returns
*         the filter number.
*
*       Arguments:
*         mission       (i) : Mission number
*         detector      (i) : Detector number
*         fl_name       (i) : Filter name string
*         <PDX_GTFIN>   (r) : 0 if given name string does not match any
*                               filter names in index,
*                             -1 if error in mission or detector number
*                               or no detectors or filters
*                             otherwise negative if ambiguous
*                               (matches -PDX_GTFIN possibilities)
*                             positive if matches one (filter number is
*                               returned in this case)
*
*       Dependencies:
*         LENTRIM, PDX_GETDT and PDX_GETFL
*
*       Origin:
*         Created as a part of the PIMMS Index structure
*
*       Author:
*         Koji Mukai (1993 Mar 20) original version
*-PDX_GTFIN

        integer fla_len, fl_len, j, pointer
        integer exact, matches, n_dtct, n_fltr
        integer match_index( max_mssn )
        integer LENTRIM, PDX_GETDT, PDX_GETFL

        if( mission .le. 0 .or. mission .gt. n_mssn ) then
          PDX_GTFIN = -1
        else
          pointer = ind_mssn( mission )
          n_dtct = st_array( pointer )
          if( n_dtct .lt. 0 .or.
     &               detector .le. 0 .or. detector .gt. n_dtct ) then
            PDX_GTFIN = -1
          else
            pointer = PDX_GETDT( mission, detector )
            n_fltr = st_array( pointer )
            if( n_fltr .lt. 0 ) then
              PDX_GTFIN = -1
            else

              fl_len = LENTRIM( fl_name )
*             Assume it has been converted to lower case alerady, now we've
*             found the actual length of the character.

              matches = 0
              exact = 0
              do j = 1, n_fltr
                pointer = PDX_GETFL( mission, detector, j )
                if( fl_name( : fl_len )
     &                       .eq. nm_array( pointer )( : fl_len ) ) then
                  matches = matches + 1
                  match_index( matches ) = j
                  fla_len = LENTRIM( nm_array( pointer ) )
                  if( fl_len .eq. fla_len ) exact = matches
                end if
              end do

              if( matches .eq. 0 ) then
                PDX_GTFIN = 0
              else if( matches .eq. 1 ) then
                PDX_GTFIN = match_index( matches )
              else
                if( exact .eq. 0 ) then
                  PDX_GTFIN = -matches
                else
                  PDX_GTFIN = match_index( exact )
                end if
              end if

            end if
          end if
        end if

        end



*+PDX_FIRST
        subroutine PDX_FIRST( mission, detector, filter )

        implicit none

        include 'pms_index.inc'

        integer mission, detector, filter

*       Description:
*         Retrieve index entry for the first mission, detector, filter
*         combination in the pms_mssn.lst file
*
*       Arguments:
*         mission      (o) : Mission number
*         detector     (o) : Detector number
*         filter       (o) : Filter number
*
*       Dependencies:
*         PDX_GETDT
*
*       Origin:
*         Cooked up by KM to handle the PIMMS-specific database.
*
*       Author:
*         Koji Mukai, 1995 Aug, Original version
*-PDX_GETDT

        integer offset, n_dtct, n_fltr
        integer PDX_GETDT

*       Assume that there is at least one mission
        mission = 1
        offset = ind_mssn( mission )
        n_dtct = st_array( offset )

        if( n_dtct .ge. 1 ) then
*         There is at least 1 detector, choose it
          detector = 1
          offset = PDX_GETDT( mission, detector )
          n_fltr = st_array( offset )
          if( n_fltr .ge. 1 ) then
*           Looks like there are filter for this mission+detector
            filter = 1
          else
*           Looks like there is no filter for this mission+detector
            filter = 0
          end if
        else
*         This mission doesn't have detector entries
          detector = 0
          filter = 0
        end if

        end
