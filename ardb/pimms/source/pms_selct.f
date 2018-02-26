*+PMS_SELCT

        subroutine PMS_SELCT
     &            ( nn, mission, detector, filter, unit, lo, hi, flag )

        implicit none

        include 'pms_index.inc'

        integer n_unit
        parameter( n_unit = 2 )

        integer nn
        integer mission, detector, filter
        integer unit
        real lo, hi
        integer flag

*       Description:
*         Looks up the ARK-FCI parameters and decides what mission etc.
*         combinations have been selected by the user
*
*       Arguments:
*         nn            (i) : Next numerical argument
*         mission       (o) : Number for chosen mission
*         detector      (o) : Number for chosen detector
*         filter        (o) : Number for chosen filter
*         unit          (o) : Code for unit, if flux is chosen
*         lo, hi        (o) : Range of flux
*         flag          (o) : Output flag
*
*       Dependencies:
*         XPI routines, PIMMS Index routines, LENTRIM
*
*       Origin:
*         Created by KM to automatically select mission using the
*         PIMMS Index structure.
*
*       Author:
*         Koji Mukai (1993 Mar 22) Original XPI version
*         Koji Mukai, 1993 March, ARK-FCI version
*         Koji Mukai, 1997 December, fudged for flux range in Angstroms
*         Koji Mukai, 2002 November, added flux density option
*-PMS_SELCT

        integer mssn_in, dtct_in, fltr_in, unit_in
        integer len_mssn, len_unit, len_spare
        integer status, pointer, record
        integer matches, matched, k
        integer next, num_c, num_f, num_i
        real lo_in, hi_in, temp
        character*( max_len ) unit_st, rang_st, spar_st
        character*( max_len ) mssn_st, dtct_st, fltr_st

        integer LENTRIM, PDX_GETDT, PDX_GETFL
        integer PDX_GTMIN, PDX_GTDIN, PDX_GTFIN
        logical PDX_CKMOD

        character*( max_len ) unit_name( n_unit ), norm_str
        character*( max_len ) flux_str, dens_str, unab_str, angs_str
        data unit_name / 'ergs', 'photons' /
        data norm_str, flux_str, dens_str, unab_str, angs_str
     &/ 'normalization', 'flux', 'density', 'unabsorbed', 'angstroms'  /

*       Initialize tentative outputs
        mssn_in = 0
        dtct_in = 0
        fltr_in = 0
        unit_in = 0
        lo_in = 0.0
        hi_in = 0.0

        call INQ_PARAM( num_c, num_f, num_i )
        next = 1

*       Obtain first parameter, should be a mission name
        if( num_c .lt. next ) then
*         No mission name given
          flag = -1

        else
*         Success --- now check what mission number, if any, this corresponds to
          call GET_PARAM_C( next, mssn_st )
          next = next + 1
          call LOCASE( mssn_st )
          mssn_in = PDX_GTMIN( mssn_st )

          if( mssn_in .lt. 0 ) then
*           Given string matches more than one misison names
            flag = -2

          else if( mssn_in .eq. 0 ) then
*           No match --- now check if user said "flux"
            len_mssn = LENTRIM( mssn_st )

*           Begin FLUX special
            if( mssn_st( : len_mssn ) .eq. flux_str( : len_mssn ) ) then
              mssn_in = 0
*             Okay, flux chosen as mission; mode sanity check

              call PARAM_C( next, unit_st,
     &                              ' Enter units (ergs or photons) >' )
              next = next + 1

*             Okay, string given, see if it matches the list of units.
              call LOCASE( unit_st )
              len_unit = LENTRIM( unit_st )
              matches = 0
              do k = 1, n_unit
                if( unit_st( : len_unit )
     &                    .eq. unit_name( k )( : len_unit ) ) then
                  matches = matches + 1
                  matched = k
                end if
              end do

              if( matches .eq. 0 ) then
*               Doesn't look like a unit name to me
                flag = -12

              else if( matches .ge. 2 ) then
*               Unit name ambiguous (shouldn't be here)
                flag = -11

              else
                unit_in = matched
*               Okay, unit correctly read --- now find energy range

                call PARAM_C( next, rang_st,
     &                 ' Enter integration range in keV (e.g. 2-10) >' )
                next = next + 1
*               Now decode range string
                call DCD_RRANG( rang_st, lo_in, hi_in, status )

                if( status .lt. 0 ) then
*                 Decode failed
                  flag = -102

                else
*                 Decode successfull --- now check for 4th parameter

                  if( num_c .lt. next ) then
*                   No 4th parameter given, okay to proceed
                    flag = 0

                  else
*                   4th/5th parameters given, check "unabsorbed" or "Angstrom"
                    do while( num_c .ge. next )
                      call GET_PARAM_C( next, spar_st )
                      call LOCASE( spar_st )
                      len_spare = LENTRIM( spar_st )
                      if( spar_st( : len_spare )
     &                              .eq. unab_str( : len_spare ) ) then
*                       Yes, user did specify unabsorbed flux
*                       Signal this with a negative code for unit
                        unit_in = -unit_in
                        flag = 0

                      else if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                       Yes, user did specify range in Angstroms
*                       Signal this with negative range numbers
                        temp = ang_kev / lo_in
                        lo_in = -ang_kev / hi_in
                        hi_in = -temp
                        flag = 0

                      else
*                       No, the string was not "unabsrobed"
                        flag = -103
                        next = num_c

*                     end of if block, did user say "unabsorbed"?
                      end if
                      next = next + 1
                    end do

*                 end of if block, was there a 4th parameter?
                  end if

*               end of if block, was range string successfully decoded?
                end if

*             end of if block, was a valid unit given?
              end if

*           End of FLUX special

*           Begin FLUX DENSITY special
            else if( mssn_st( : len_mssn )
     &                                .eq. dens_str( : len_mssn ) ) then
              mssn_in = 131072
*             Okay, flux density chosen as mission; mode sanity check

              call PARAM_C( next, unit_st,
     &                              ' Enter units (ergs or photons) >' )
              next = next + 1

*             Okay, string given, see if it matches the list of units.
              call LOCASE( unit_st )
              len_unit = LENTRIM( unit_st )
              matches = 0
              do k = 1, n_unit
                if( unit_st( : len_unit )
     &                    .eq. unit_name( k )( : len_unit ) ) then
                  matches = matches + 1
                  matched = k
                end if
              end do

              if( matches .eq. 0 ) then
*               Doesn't look like a unit name to me
                flag = -12

              else if( matches .ge. 2 ) then
*               Unit name ambiguous (i.e., 'm' )
                flag = -11

              else
                unit_in = matched
*               Okay, unit correctly read --- now find the energy

                call PARAM_N( nn, lo_in, 1.0, 'Enter energy in keV >' )
                hi_in = -999.9

*               Now check for 4th parameter

                if( num_c .lt. next ) then
*                 No 4th parameter given, okay to proceed
                  flag = 0

                else
*                 4th/5th parameters given, check "unabsorbed" or "Angstrom"
                  do while( num_c .ge. next )
                    call GET_PARAM_C( next, spar_st )
                    call LOCASE( spar_st )
                    len_spare = LENTRIM( spar_st )
                    if( spar_st( : len_spare )
     &                              .eq. unab_str( : len_spare ) ) then
*                     Yes, user did specify unabsorbed flux density
*                     Signal this with a negative code for unit
                      unit_in = -unit_in
                      flag = 0

                    else if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                     Yes, user did specify energy in Angstroms
*                     Signal this with negative range numbers
                      lo_in = -ang_kev / lo_in
                      flag = 0

                    else
*                     No, the string was not "unabsrobed"
                      flag = -103
                      next = num_c

*                   end of if block, did user say "unabsorbed"?
                    end if
                    next = next + 1
                  end do

*               end of if block, was there a 4th parameter?
                end if

*             end of if block, was a valid unit given?
              end if

*           End of FLUX DENSITY special

            else if( mssn_st( : len_mssn )
     &                                .eq. norm_str( : len_mssn ) ) then
*             User said "normalization"
              mssn_in = 65536

            else
*             Mission name did not match any in index, and was not "flux"
*             or "normalization"
              flag = -3

*           end of if block, does the mission name actually say "flux"?
            end if

          else
*           Mission name did match one in index --- now what's in index?
            pointer = ind_mssn( mssn_in )
            record = st_array( pointer )

            if( record .lt. 0 ) then
*             This mission has no subdivisions into detectors or filters !
*             Do a sanity check on modes now.

              if( PDX_CKMOD( record ) ) then
*               Mode and Index record are compatible --- check for range

                if( num_c .lt. next ) then
*                 Energy range not given --- this is fine
                  flag = 0

                else
*                 Energy range was given --- decode it
                  call GET_PARAM_C( next, rang_st )
                  next = next + 1
                  call DCD_RRANG( rang_st, lo_in, hi_in, status )

                  if( status .ne. 0 ) then
*                   Range decode failed
                    flag = -102

                  else
*                   Range decoded successfully --- everything's A-OK
                    if( num_c .lt. next ) then
*                     No extra parameter given, okay to proceed
                      flag = 0
                    else
*                     extra parameter given, check if "Angstrom"
                      call GET_PARAM_C( next, spar_st )
                      call LOCASE( spar_st )
                      len_spare = LENTRIM( spar_st )
                      if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                       Yes, user did specify range in Angstroms
*                       Signal this with negative range numbers
                        temp = ang_kev / lo_in
                        lo_in = -ang_kev / hi_in
                        hi_in = -temp
                        flag = 0

                      else
*                       No, the string was not "Angstroms"
                        flag = -103

*                     end of if block, did user say "Angstroms"?
                      end if

*                   end of if block, was there an extra parameter?
                    end if

*                 end of if block, was range successfully decoded?
                  end if

*               end of if block, was a range string given?
                end if

              else
*               mode and Index record are incompatible
                flag = -22

*             end of if block, are mode and index record compatible?
              end if

            else if( record .eq. 1 ) then
*             This mission has 1 detector --- check 2nd parameter

              if( num_c .lt. next ) then
*               No second parameter given, which may be okay
                dtct_in = 1
                pointer = PDX_GETDT( mssn_in, dtct_in )
                record = st_array( pointer )

                if( record .gt. 0 ) then
*                 Ah, but there are more than 1 filters
                  flag = -7

                else
*                 Do a sanity check on modes

                  if( PDX_CKMOD( record ) ) then
*                   Mode and Index record are compatible; everything A-OK
                    flag = 0

                  else
*                   Mode and Index record are incompatible.
                    flag = -22

*                 end of if block, is mode compatible with index record?
                  end if

*               end of if block, are there multiple filters
                end if

              else
*               There was a second string --- was it the detector name?
                call GET_PARAM_C( next, dtct_st )
                next = next + 1
                call LOCASE( dtct_st )
                pointer = PDX_GETDT( mssn_in, 1 )
                record = st_array( pointer )
                dtct_in = PDX_GTDIN( mssn_in, dtct_st )

                if( dtct_in .lt. 0 ) then
*                 Error/confusion
                  flag = -104

                else if( dtct_in .eq. 0 ) then
*                 Was not the detector name; should the user given the
*                 detector or filter name?

                  if( record .gt. 1 ) then
*                   Yes, he should have because there are more than 1 filters
*                   ...but let's allow a shorthand of not giving the detector
*                   name and going straight to the filter name!
                    dtct_in = 1
                    fltr_st = dtct_st
                    fltr_in = PDX_GTFIN( mssn_in, dtct_in, fltr_st )

                    if( fltr_in .lt. -1 ) then
*                     Ambiguous filter name
                      flag = -8

                    else if( fltr_in .eq. -1 ) then
*                     Index seems to be confused
                      flag = -104

                    else if( fltr_in .eq. 0 ) then
*                     No match for filter name either
                      flag = -7

                    else
*                     Okay, this does match a filter name; Sanity check
                      pointer = PDX_GETFL( mssn_in, dtct_in, fltr_in )
                      record = st_array( pointer )

                      if( PDX_CKMOD( record ) ) then
*                       Mode and Index record are compatible

                        if( num_c .lt. next ) then
*                         No range string given, which is fine.
                          flag = 0

                        else
*                         There was a range string; now try to read it
                          call GET_PARAM_C( next, rang_st )
                          next = next + 1
                          call DCD_RRANG
     &                                ( rang_st, lo_in, hi_in, status )

                          if( status .lt. 0 ) then
*                           No, could not decode the range string
                            flag = -102

                          else
*                           Yes, it was a valid range string
                            if( num_c .lt. next ) then
*                             No extra parameter given, okay to proceed
                              flag = 0
                            else
*                             extra parameter given, check if "Angstrom"
                              call GET_PARAM_C( next, spar_st )
                              call LOCASE( spar_st )
                              len_spare = LENTRIM( spar_st )
                              if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                               Yes, user did specify range in Angstroms
*                               Signal this with negative range numbers
                                temp = ang_kev / lo_in
                                lo_in = -ang_kev / hi_in
                                hi_in = -temp
                                flag = 0

                              else
*                               No, the string was not "Angstroms"
                                flag = -103

*                             end of if block, did user say "Angstroms"?
                              end if

*                           end of if block, was there an extra parameter?
                            end if

*                         end of if block, was it a valid range string
                          end if

*                       end of if block, did the user give a range string?
                        end if

                      else
*                       No, the mode and Index record are incompatible
                        flag = -22

*                     end of if block, was mode and Index record compatible?
                      end if

*                   end of if block, does it match a filter name
                    end if

                  else if( record .eq. 1 ) then
*                   There is a single filter for this single detector
                    dtct_in = 1
                    fltr_st = dtct_st
                    fltr_in = PDX_GTFIN( mssn_in, dtct_in, fltr_st )

                    if( fltr_in .lt. -1 ) then
*                     Ambiguous filter name
                      flag = -8

                    else if( fltr_in .eq. -1 ) then
*                     Index seems to be confused
                      flag = -104

                    else if( fltr_in .eq. 0 ) then
*                     No, it's not filter name either, but then the user need
*                     not input this.  Check if it's the range string instead
                      fltr_in = 1
                      rang_st = dtct_st
                      call DCD_RRANG( rang_st, lo_in, hi_in, status )

                      if( status .lt. 0 ) then
*                       It's not a valid range string
                        flag = -102

                      else
*                       Yes, it is --- but wait, have to do the sanity check
                        pointer = PDX_GETFL( mssn_in, dtct_in, fltr_in )
                        record = st_array( pointer )

                        if( PDX_CKMOD( record ) ) then
*                         Yes, mode and Index record are compatible
                          if( num_c .lt. next ) then
*                           No extra parameter given, okay to proceed
                            flag = 0
                          else
*                           extra parameter given, check if "Angstrom"
                            call GET_PARAM_C( next, spar_st )
                            call LOCASE( spar_st )
                            len_spare = LENTRIM( spar_st )
                            if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                             Yes, user did specify range in Angstroms
*                             Signal this with negative range numbers
                              temp = ang_kev / lo_in
                              lo_in = -ang_kev / hi_in
                              hi_in = -temp
                              flag = 0

                            else
*                             No, the string was not "Angstroms"
                              flag = -103

*                           end of if block, did user say "Angstroms"?
                            end if

*                         end of if block, was there an extra parameter?
                          end if

                        else
*                         No, mode and Index record are incompatible
                          flag = -22

*                       end of if block, was mode and Index record compatible?
                        end if

*                     end of if block, was it a valid range string?
                      end if

                    else
*                     Yes, it was the filter name --- funny but valid I suppose
*                     Now let's do the sanity check
                      pointer = PDX_GETFL( mssn_in, dtct_in, fltr_in )
                      record = st_array( pointer )

                      if( PDX_CKMOD( record ) ) then
*                       Yes, mode and Index record are compatible

                        if( num_c .lt. next ) then
*                         No range string give --- everything is okay
                          flag = 0

                        else
*                         There was a range string --- now try to decode it
                          call GET_PARAM_C( next, rang_st )
                          next = next + 1
                          call DCD_RRANG
     &                                 ( rang_st, lo_in, hi_in, status )

                          if( status .lt. 0 ) then
*                           Failed to decode
                            flag = -102

                          else
*                           Succeeded to decode
                            if( num_c .lt. next ) then
*                             No extra parameter given, okay to proceed
                              flag = 0
                            else
*                             extra parameter given, check if "Angstrom"
                              call GET_PARAM_C( next, spar_st )
                              call LOCASE( spar_st )
                              len_spare = LENTRIM( spar_st )
                              if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                               Yes, user did specify range in Angstroms
*                               Signal this with negative range numbers
                                temp = ang_kev / lo_in
                                lo_in = -ang_kev / hi_in
                                hi_in = -temp
                                flag = 0

                              else
*                               No, the string was not "Angstroms"
                                flag = -103

*                             end of if block, did user say "Angstroms"?
                              end if

*                           end of if block, was there an extra parameter?
                            end if

*                         end of if block, was the range string valid?
                          end if

*                       end of if block, was there a range string?
                        end if

                      else
*                       No, mode and Index record are incompatible
                        flag = -22

*                     end of if block, mode and Index record compatible?
                      end if

*                   end of if block, was it the filter name?
                    end if

                  else
*                   No need whatsoever for detector or filter names;
*                   interpret the second string as a range.
                    rang_st = dtct_st
                    call DCD_RRANG( rang_st, lo_in, hi_in, status )

                    if( status .ne. 0 ) then
*                     It wasn't a range, either --- output range decode error
                      flag = -102

                    else
*                     Yes, it was a range --- so it's detector 1 of 1 okay
                      dtct_in = 1
*                     Do a sanity check on modes

                      if( PDX_CKMOD( record ) ) then
*                       Mode and Index record are compatible; everything A-OK
                        if( num_c .lt. next ) then
*                         No extra parameter given, okay to proceed
                          flag = 0
                        else
*                         extra parameter given, check if "Angstrom"
                          call GET_PARAM_C( next, spar_st )
                          call LOCASE( spar_st )
                          len_spare = LENTRIM( spar_st )
                          if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                           Yes, user did specify range in Angstroms
*                           Signal this with negative range numbers
                            temp = ang_kev / lo_in
                            lo_in = -ang_kev / hi_in
                            hi_in = -temp
                            flag = 0

                          else
*                           No, the string was not "Angstroms"
                            flag = -103

*                         end of if block, did user say "Angstroms"?
                          end if

*                       end of if block, was there an extra parameter?
                        end if

                      else
*                       Mode and Index record are incompatible.
                        flag = -22

*                     end of if block, is mode compatible with index record?
                      end if

*                   Is it a valid range string?
                    end if

*                 How many filters does this single detector have?
                  end if

                else
*                 Was the detector name; how many filters does it have?
                  pointer = PDX_GETDT( mssn_in, dtct_in )
                  record = st_array( pointer )

                  if( record .gt. 0 ) then
*                   Multiple filters for this mission/detector

                    if( num_c .lt. next ) then
*                     No filter name given
                      flag = -7

                    else
*                     Filter name was given; check validity
                      call GET_PARAM_C( next, fltr_st )
                      next = next + 1
                      call LOCASE( fltr_st )
                      fltr_in = PDX_GTFIN( mssn_in, dtct_in, fltr_st )

                      if( fltr_in .lt. -1 ) then
*                       Ambiguous
                        flag = -8

                      else if( fltr_in .eq. -1 ) then
*                       PIMMS Index is confused
                        flag = -104

                      else if( fltr_in .eq. 0 ) then
*                       Filter name not valid
                        flag = -9

                      else
*                       A correct filter name was given; Sanity check
                        pointer
     &                   = PDX_GETFL( mssn_in, dtct_in, fltr_in )
                        record = st_array( pointer )

                        if( PDX_CKMOD( record ) ) then
*                         Mode and Index record are compatible, range?

                          if( num_c .lt. next ) then
*                           No range was given --- everything is A-OK.
                            flag = 0

                          else
*                           Range was given; try to decode it
                            call GET_PARAM_C( next, rang_st )
                            next = next + 1
                            call DCD_RRANG
     &                                 ( rang_st, lo_in, hi_in, status )

                            if( status .lt. 0 ) then
*                             Range decode failed
                              flag = -102

                            else
*                             Range successfully decoded
                              if( num_c .lt. next ) then
*                               No extra parameter given, okay to proceed
                                flag = 0
                              else
*                               extra parameter given, check if "Angstrom"
                                call GET_PARAM_C( next, spar_st )
                                call LOCASE( spar_st )
                                len_spare = LENTRIM( spar_st )
                                if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                                 Yes, user did specify range in Angstroms
*                                 Signal this with negative range numbers
                                  temp = ang_kev / lo_in
                                  lo_in = -ang_kev / hi_in
                                  hi_in = -temp
                                  flag = 0

                                else
*                                 No, the string was not "Angstroms"
                                  flag = -103

*                               end of if block, did user say "Angstroms"?
                                end if

*                             end of if block, was there an extra parameter?
                              end if

*                           end of if block, was the range successfully decoded?
                            end if

*                         end of if block, was the range string given?
                          end if

                        else
*                         Mode and Index record are incompatible
                          flag = -22

*                       end of if block, mode and Index record compatible?
                        end if

*                     end of if block, was the filter name valid?
                      end if

*                   end of if block, was the filter name given?
                    end if

                  else
*                   No filters; sanity check first

                    if( PDX_CKMOD( record ) ) then
*                     Mode and Index record are compatible; was there range?

                      if( num_c .lt. next ) then
*                       No range string given, which is fine
                        flag = 0

                      else
*                       Range string was found; is it decodable?
                        call GET_PARAM_C( next, rang_st )
                        next = next + 1
                        call DCD_RRANG( rang_st, lo_in, hi_in, status )

                        if( status .lt. 0 ) then
*                         DCD_RRANG failed
                          flag = -102

                        else
*                         Successful decode - everything is fine
                          if( num_c .lt. next ) then
*                           No extra parameter given, okay to proceed
                            flag = 0
                          else
*                           extra parameter given, check if "Angstrom"
                            call GET_PARAM_C( next, spar_st )
                            call LOCASE( spar_st )
                            len_spare = LENTRIM( spar_st )
                            if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                             Yes, user did specify range in Angstroms
*                             Signal this with negative range numbers
                              temp = ang_kev / lo_in
                              lo_in = -ang_kev / hi_in
                              hi_in = -temp
                              flag = 0

                            else
*                             No, the string was not "Angstroms"
                              flag = -103

*                           end of if block, did user say "Angstroms"?
                            end if

*                         end of if block, was there an extra parameter?
                          end if

*                       end of if block, was range string valid?
                        end if

*                     end of if block, was range string given?
                      end if

                    else
*                     Mode and Index record are incompatible
                      flag = -22

*                   end of if block, mode and Index record compatible?
                    end if

*                 end of if block, any filters?
                  end if

*               end of if block, was the 2nd parameter the detector name?
                end if

*             end of if block, was there a 2nd parameter?
              end if

            else
*             This mission has a healthy number of detectors

              if( num_c .lt. next ) then
*               No detector name given
                flag = -4

              else
*               Yes, there is a detector string; is it valid?
                call GET_PARAM_C( next, dtct_st )
                next = next + 1
                call LOCASE( dtct_st )
                dtct_in = PDX_GTDIN( mssn_in, dtct_st )

                if( dtct_in .lt. -1 ) then
*                 Ambiguous name
                  flag = -5

                else if( dtct_in .eq. -1 ) then
*                 Index structure is confused
                  flag = -104

                else if( dtct_in .eq. 0 ) then
*                 No match for the detector name given
                  flag = -6

                else
*                 Valid detector name given; how many filters are there?
                  pointer = PDX_GETDT( mssn_in, dtct_in )
                  record = st_array( pointer )

                  if( record .lt. 0 ) then
*                   There are no filters; sanity check first

                    if( PDX_CKMOD( record ) ) then
*                     Mode and Index record are compatible; check for range

                      if( num_c .lt. next ) then
*                       There was no range string, which is fine
                        flag = 0

                      else
*                       Range string was found; now try to decode it
                        call GET_PARAM_C( next, rang_st )
                        next = next + 1
                        call DCD_RRANG( rang_st, lo_in, hi_in, status )

                        if( status .lt. 0 ) then
*                         Range string could not be decoded
                          flag = -102

                        else
*                         Range string was decoded --- everything's A-OK
                          if( num_c .lt. next ) then
*                           No extra parameter given, okay to proceed
                            flag = 0
                          else
*                           extra parameter given, check if "Angstrom"
                            call GET_PARAM_C( next, spar_st )
                            call LOCASE( spar_st )
                            len_spare = LENTRIM( spar_st )
                            if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                             Yes, user did specify range in Angstroms
*                             Signal this with negative range numbers
                              temp = ang_kev / lo_in
                              lo_in = -ang_kev / hi_in
                              hi_in = -temp
                              flag = 0

                            else
*                                 No, the string was not "Angstroms"
                              flag = -103

*                           end of if block, did user say "Angstroms"?
                            end if

*                         end of if block, was there an extra parameter?
                          end if

*                       end of if block, was range string decoded?
                        end if

*                     end of if block, was there a range string?
                      end if

                    else
*                     Mode and Index record are incompatible
                      flag = -22

*                   end of if block, mode and Index record consistent?
                    end if

                  else if( record .eq. 1 ) then
*                   There is 1 filter (a strange state of affiars, but
*                   there you go).  Check for compatibility first.
                    fltr_in = 1
                    pointer = PDX_GETFL( mssn_in, dtct_in, fltr_in )
                    record = st_array( pointer )

                    if( PDX_CKMOD( record ) ) then
*                     Mode and Index record are compatible; check for range

                      if( num_c .lt. next ) then
*                       There was no filter/range string, which is fine
                        flag = 0

                      else
*                       filter/range string was found; try to match it
                        call GET_PARAM_C( next, fltr_st )
                        next = next + 1
                        call LOCASE( fltr_st )
                        fltr_in = PDX_GTFIN( mssn_in, dtct_in, fltr_st )

                        if( fltr_in .le. -1 ) then
*                         PIMMS Index must be confused
                          flag = -103

                        else if( fltr_in .eq. 0 ) then
*                         That was not a filter name; reset fltr_in and
*                         try treating it as the range string
                          fltr_in = 1
                          rang_st = fltr_st
                          call DCD_RRANG
     &                                 ( rang_st, lo_in, hi_in, status )

                          if( status .lt. 0 ) then
*                           Range string could not be decoded
                            flag = -102

                          else
*                           Range string was decoded --- everything's A-OK
                            if( num_c .lt. next ) then
*                             No extra parameter given, okay to proceed
                              flag = 0
                            else
*                             extra parameter given, check if "Angstrom"
                              call GET_PARAM_C( next, spar_st )
                              call LOCASE( spar_st )
                              len_spare = LENTRIM( spar_st )
                              if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                               Yes, user did specify range in Angstroms
*                               Signal this with negative range numbers
                                temp = ang_kev / lo_in
                                lo_in = -ang_kev / hi_in
                                hi_in = -temp
                                flag = 0

                              else
*                               No, the string was not "Angstroms"
                                flag = -103

*                             end of if block, did user say "Angstroms"?
                              end if

*                           end of if block, was there an extra parameter?
                            end if

*                         end of if block, was range string decoded?
                          end if

                        else
*                         User did give the correct filter name; now try range

                          if( num_c .lt. next ) then
*                           No range string given, which is fine
                            flag = 0

                          else
*                           There was a range string; now decode it
                            call GET_PARAM_C( next, rang_st )
                            next = next + 1
                            call DCD_RRANG
     &                                 ( rang_st, lo_in, hi_in, status )

                            if( status .lt. 0 ) then
*                             Range string could not be decoded
                              flag = -102

                            else
*                             Range string was decoded --- everything's A-OK
                              if( num_c .lt. next ) then
*                               No extra parameter given, okay to proceed
                                flag = 0
                              else
*                               extra parameter given, check if "Angstrom"
                                call GET_PARAM_C( next, spar_st )
                                call LOCASE( spar_st )
                                len_spare = LENTRIM( spar_st )
                                if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                                 Yes, user did specify range in Angstroms
*                                 Signal this with negative range numbers
                                  temp = ang_kev / lo_in
                                  lo_in = -ang_kev / hi_in
                                  hi_in = -temp
                                  flag = 0

                                else
*                                 No, the string was not "Angstroms"
                                  flag = -103

*                               end of if block, did user say "Angstroms"?
                                end if

*                             end of if block, was there an extra parameter?
                              end if

*                           end of if block, was range string decoded?
                            end if

*                         end of if block, was there a range string?
                          end if

*                       end of if block, is it the valid filter name?
                        end if

*                     end of if block, was the 3rd string given?
                      end if

                    else
*                     Mode and Index record are incompatible
                      flag = -22

*                   end of if block, mode and Index record consistent?
                    end if

                  else
*                   More than 1 filters, filter name is now compulsory

                    if( num_c .lt. next ) then
*                     But filter name was not given
                      flag = -7

                    else
*                     Okay, filter string was found; try to match it
                      call GET_PARAM_C( next, fltr_st )
                      next = next + 1
                      call LOCASE( fltr_st )
                      fltr_in = PDX_GTFIN( mssn_in, dtct_in, fltr_st )

                      if( fltr_in .lt. -1 ) then
*                       Filter name was ambiguous
                        flag = -8

                      else if( fltr_in .eq. -1 ) then
*                       Oh dear, PIMMS Index is confused
                        flag = -104

                      else if( fltr_in .eq. 0 ) then
*                       String does not match filter names
                        flag = -9

                      else
*                       Okay, found the filter name; now sanity check
                        pointer = PDX_GETFL( mssn_in, dtct_in, fltr_in )
                        record = st_array( pointer )

                        if( PDX_CKMOD( record ) ) then
*                         Yes, mode and Index record are consistent; range

                          if( num_c .lt. next ) then
*                           No range string, which is just fine
                            flag = 0

                          else
*                           Range string was found; try decoding
                            call GET_PARAM_C( next, rang_st )
                            next = next + 1
                            call DCD_RRANG
     &                                 ( rang_st, lo_in, hi_in, status )

                            if( status .lt. 0 ) then
*                             No, could not decode it
                              flag = -102

                            else
*                             Yes, range was found, everything is A-OK
                              if( num_c .lt. next ) then
*                               No extra parameter given, okay to proceed
                                flag = 0
                              else
*                               extra parameter given, check if "Angstrom"
                                call GET_PARAM_C( next, spar_st )
                                call LOCASE( spar_st )
                                len_spare = LENTRIM( spar_st )
                                if( spar_st( : len_spare )
     &                              .eq. angs_str( : len_spare ) ) then
*                                 Yes, user did specify range in Angstroms
*                                 Signal this with negative range numbers
                                  temp = ang_kev / lo_in
                                  lo_in = -ang_kev / hi_in
                                  hi_in = -temp
                                  flag = 0

                                else
*                                 No, the string was not "Angstroms"
                                  flag = -103

*                               end of if block, did user say "Angstroms"?
                                end if

*                             end of if block, was there an extra parameter?
                              end if

*                           end of if block, was it a valid range string?
                            end if

*                         end of if block, was there a range string?
                          end if

                        else
*                         No, mode and Index record are inconsistent
                          flag = -22

*                       end of if block, mode and Index record consistent?
                        end if

*                     end of if block, does the string match a filter name?
                      end if

*                   end of if block, was a filter name string given?
                    end if

*                 end of if block, how many filters does this detector have?
                  end if

*               end of if block, was the detector name valid?
                end if

*             end of if block, was detector name given?
              end if

*           end of if block, how many detectors does this mission have?
            end if

*         end of if block, was the first parameter a valid mission name?
          end if

*       end of if block, was a parameter given?
        end if

*       Now output error messages or assign value to the arguments

        if( flag .eq. -1 ) then
*         No mission name given
          call PWRITE( 'No mission name (or ''flux'') given' )
          call PWRITE(
     &          'Type ''directory'' for a list of valid mission names' )
        else if( flag .eq. -2 ) then
*         Given string matches more than one misison names
          call PWRITE( 'Ambiguousd mission name' )
          call PWRITE(
     &          'Type ''directory'' for a list of valid mission names' )
        else if( flag .eq. -3 ) then
*         Given string does not match any misison names
          call PWRITE( 'Unknown mission name' )
          call PWRITE(
     &          'Type ''directory'' for a list of valid mission names' )
        else if( flag .eq. -4 ) then
*         No detector name given, where one is necessary
          call PWRITE( 'No detector name given' )
          call PWRITE(
     &          'Type ''directory'' for a list of valid mission names' )
        else if( flag .eq. -5 ) then
*         Given string matches more than one detector names
          call PWRITE( 'Ambiguousd detector name' )
          call PWRITE(
     &          'Type ''directory'' for a list of valid mission names' )
        else if( flag .eq. -6 ) then
*         Given string does not match any detector names
          call PWRITE( 'Unknown detector name' )
          call PWRITE(
     &          'Type ''directory'' for a list of valid mission names' )
        else if( flag .eq. -7 ) then
*         No filter name given, where one is necessary
          call PWRITE( 'No filter name given' )
          call PWRITE(
     &          'Type ''directory'' for a list of valid mission names' )
        else if( flag .eq. -8 ) then
*         Given string matches more than one filter names
          call PWRITE( 'Ambiguousd filter name' )
          call PWRITE(
     &          'Type ''directory'' for a list of valid mission names' )
        else if( flag .eq. -9 ) then
*         Given string does not match any filter names
          call PWRITE( 'Unknown filter name' )
          call PWRITE(
     &          'Type ''directory'' for a list of valid mission names' )
        else if( flag .eq. -11 ) then
*         Unit name ambiguous (i.e., 'm' )
          call PWRITE( 'Ambiguous unit of flux; use ergs or photons' )
        else if( flag .eq. -12 ) then
*         Doesn't look like a unit name to me
          call PWRITE( 'Unknown unit of flux; use ergs or photons' )

        else if( flag .eq. -21 ) then
*         Image simulation with flux as output !!
          call PWRITE(
     &            'Image simulation needs a mission, not ''flux''' )
        else if( flag .eq. -22 ) then
*         Simulation type and index mismatch
          call PWRITE( 'PIMMS does not know how to calculate '
     &                   // 'this data product for this instrument' )

        else if( flag .eq. -101 ) then
*         XPI reports an error trying to read parameters
          call PWRITE(
     &       'ERROR in PMS_SELECT:: parameter read error from XPI' )
        else if( flag .eq. -102 ) then
*         DCD_RRANG reports error
          call PWRITE( 'ERROR in PMS_SELCT:: garbled energy range' )
        else if( flag .eq. -103 ) then
*         No, the string was not "unabsrobed"
          call PWRITE( 'ERROR in PMS_SELCT:: unexpected extra '
     &                                     // 'parameter is given' )
        else if( flag .eq. -104 ) then
*         Error/confusion in PIMMS index
          call PWRITE(
     &     'SEVERE ERROR in PMS_SELCT:: Mission Index is confused' )

        else if( flag .eq. 0 ) then
*         Was success --- note that on any error conditions, previous
*         values of mission etc. are kept.

          mission = mssn_in
          detector = dtct_in
          filter = fltr_in
          unit = unit_in
          lo = lo_in
          hi = hi_in

        end if

        end
