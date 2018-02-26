*+PMS_SPECL
        subroutine PMS_SPECL
     &              ( mission, detector, filter, restr, results, n_res )

        implicit none

        include 'pms_index.inc'

        integer mission, detector, filter
        logical restr
        integer n_res
        real results( 0: n_res )

*       Description:
*         Given the predicted count rate, this drives the very instrument
*         specific subroutines that displays limits etc.  The routine
*         down in the tree is the one that's hard-wired, this one is NOT.
*
*       Arguments:
*         mission   (i) : mission number
*         detector  (i) : detector number
*         filter    (i) : filter number
*         restr     (i) : whether the output energy range was specified
*         results   (i) : predicted count rate(s)
*         n_es      (i) : number of extra results passed onto it
*
*       Dependencies:
*         PIMMS Index routines
*
*       Origin:
*         Created by KM to minimize the hard-wiring within PIMMS
*
*       Author:
*         Koji Mukai 1993 Mar, original version
*              Always adding new missions.
*         Latest mode: 2001 July, updating XMM for AO-2
*-PMS_SPECL

        integer pointer, record, status, flag, temp_unit, spec_num
        character*256 temp_name

        integer PDX_GETDT, PDX_GETFL

        include 'sitespec.inc'

        pointer = ind_mssn( mission )
        record = st_array( pointer )
        if( record .gt. 0 ) then
          pointer = PDX_GETDT( mission, detector )
          record = st_array( pointer )
          if( record .gt. 0 ) then
            pointer = PDX_GETFL( mission, detector, filter )
            record = st_array( pointer )
          end if
        end if
        record = abs( record )
        if( mod( record / 100, 10 ) .eq. 1 ) then
          status = 0
          call PDX_MKNAM( mission, detector, filter, 'SPECIAL',
     &                                            temp_name, flag )
          if( flag .ge. 0 ) then
            call ARKOPN( temp_unit, ddir_name, temp_name, 'SPECIAL',
     &                   'OLD', 'READONLY', 'FORMATTED', 'SEQUENTIAL',
     &                   1, flag )
            read( temp_unit, * ) spec_num
            call PMS_HDWRD( spec_num, restr, results, n_res, temp_unit )
            close( temp_unit )
          end if
        end if

        end

*+PMS_HDWRD
        subroutine PMS_HDWRD( special, restr, results, n_res, s_lun )

        implicit none

        integer special
        logical restr
        integer n_res
        real results( 0: n_res )
        integer s_lun

*       This is the interface between general PIMMS routines
*       and instrument specific routines, via the hardwired code
*
*       Modified in 2005 Nov, to include update Suzaku support
*       Modified in 2006 Apr, to updated Swift/BAT support
*       Modified in 2006 May to amend Swift/XRT support
*       Modified in 2006 Oct to update Suzaku support
*       Modified in 2007 Sep now that EPIC pn rate is PATTERN=0 only
*       Modified in 2008 Sep adding warning about MOS1 timing mode data
*       Modified in 2008 Oct with Suzaku/HXD systematic detection limits
*       Modified in 2008 Oct to clarify (somewhat) "further info" message
*       Modified in 2010 Jan to correct the "ge" vs "eq" bug
*       Modified in 2010 Oct to delete "HXD nominal" comments from Suzaku
*       Modified in 2011 Feb to include Integral JEM-X, which now has
*                   operational units.
*       Last modified in 2001 Apr to include NuSTAR and ASTRO-H SXS

        character*78 pw_strng

*       Current "SPECIAL"s
*         101 & 102:	ASCA SIS & GIS
*	  201 - 215:	XTE (PCA, old HEXTE, ASM; XTE with LLDs)
*         304:          HEAO-1 A4
*         401 - 404:    SAX (LECS, MECS, HPGSPC and PDS)
*         501 - 531:    CHANDRA (ACIS pileup)
*         601, 612-614:	XMM (MOS and PN pile-up)
*         701-722:      ASTRO-E/Astro-E2/Suzaku
*         801 & 821:    Integral
*         901-913:      Swift
*         1001:         NuSTAR
*         1101-1104:    ASTRO-H SXS
*         1131:         ASTRO-H SGD

        if( special .eq. 101 ) then
*         This hard-wired code means we're dealing with ASCA SIS
          call PWRITE( '    [Count rate per single SIS, not per ' //
     &                            'pair, over an entire chip assuming' )
          call PWRITE( '    a point source at the 1-CCD mode ' //
     &                         'position; within the maximum circular' )
          call PWRITE( '    extraction region that fits on the ' //
     &                                'default chip, usable rates are' )
          pw_Strng = '              (3.2 arcmin radius, SIS-0) and' //
     &                                 '           (2.5 arcmin, SIS-1)]'
          if( results( 1 ) .gt. 1.0
     &                           .and. results( 1 ) .le. 500000.0 ) then
            write( pw_strng( 5: 13 ), '(f9.1)' ) results( 1 )
          else
            write( pw_strng( 5: 13 ), '(1p,e9.2)' ) results( 1 )
          end if
          if( results( 2 ) .gt. 1.0
     &                           .and. results( 2 ) .le. 500000.0 ) then
            write( pw_strng( 46: 54 ), '(f9.1)' ) results( 2 )
          else
            write( pw_strng( 46: 54 ), '(1p,e9.2)' ) results( 2 )
          end if
          call PWRITE( pw_strng )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call SIS_LIMIT( results( 0 ) )
          end if
        else if( special .eq. 102 ) then
*         This hard-wired code means we're dealing with ASCA GIS
          call PWRITE( '    [Count rate per single GIS, not per ' //
     &                           'pair, but over the entire detector.' )
          call PWRITE( '    For point source at the 1-CCD mode ' //
     &                            'position and a 24 pixel (6 arcmin)' )
          pw_Strng = '    extraction radius, usable rates are ' //
     &                          '          (GIS2) and           (GIS3)]'
          if( results( 1 ) .gt. 1.0
     &                           .and. results( 1 ) .le. 500000.0 ) then
            write( pw_strng( 41: 49 ), '(f9.1)' ) results( 1 )
          else
            write( pw_strng( 41: 49 ), '(1p,e9.2)' ) results( 1 )
          end if
          if( results( 2 ) .gt. 1.0
     &                           .and. results( 2 ) .le. 500000.0 ) then
            write( pw_strng( 62: 70 ), '(f9.1)' ) results( 2 )
          else
            write( pw_strng( 62: 70 ), '(1p,e9.2)' ) results( 2 )
          end if
          call PWRITE( pw_strng )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call GIS_LIMIT( results( 0 ) )
          end if

        else if( special .eq. 201 ) then
*         This hard-wired code means we're dealing with XTE PCA
          call PWRITE( '   (Count rate is per PCU)' )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call PCA_LIMIT( results, n_res )
          end if
        else if( special .eq. 202
     &           .or. ( special .ge. 210 .and. special .le. 215 ) ) then
*         This hard-wired code means we're dealing with XTE HEXTE
          if( restr ) then
            call PWRITE( '   (Source-only count rate in 1 cluster)' )
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call HXT_LIMIT( results, n_res, special )
          end if
        else if( special .eq. 203 ) then
*         This hard-wired code means we're dealing with XTE ASM
          call PWRITE( '   (Source-only count rate for a single ' //
     &                                 'SSC, with the source on-axis)' )



        else if( special .eq. 304 ) then
*         This hard-wired code means we're dealing with HEAO-1 A4
          call PWRITE( '   (Count rate in Levin et al (1984, ApJS ' //
     &                                          '54, 581) band A+B+C)' )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call HA4_LEVIN( results, n_res )
          end if


        else if ( special .eq. 401 ) then
*         This hard-wired code means we're dealing with SAX LECS
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call SAX_LIMIT( special, results, n_res )
          end if
        else if ( special .eq. 402 ) then
*         This hard-wired code means we're dealing with SAX MECS
          call PWRITE( '   (Source count rate for 2 MECS)' )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call SAX_LIMIT( special, results( 0 ), n_res )
          end if
        else if ( special .eq. 403 ) then
*         This hard-wired code means we're dealing with SAX PDS
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call SAX_LIMIT( special, results, n_res )
          end if
        else if ( special .eq. 404 ) then
*         This hard-wired code means we're dealing with SAX HPGSPC
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call SAX_LIMIT( special, results, n_res )
          end if

        else if( special .ge. 501 .and. special .le. 502 ) then
***        else if( special .ge. 501 .and. special .le. 503 ) then
*         Hard-wired code 501 means we're dealing with Chandra ACIS-I
*         Hard-wired code 502 means we're dealing with Chandra ACIS-S-BI
***       Hard-wired code 503 used to mean we're dealing with Chandra ACIS-S FI
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call PWRITE( '% Pileup estimate for ACIS:' )
            call ACIS_CHOICE( results( 0 ) )
          end if
        else if( special .eq. 507 .or. special .eq. 511 .or.
     &           special .eq. 527 .or. special .eq. 531 ) then
*         Hard-wired code 507: we're dealing with Chandra HETG-ACIS-S ORDER0
*         Hard-wired code 511: we're dealing with Chandra LETG-ACIS-S ORDER0
*         Hard-wired code 527: we're dealing with Chandra HETG-ACIS-I ORDER0
*         Hard-wired code 531: we're dealing with Chandra LETG-ACIS-I ORDER0
*         Just run PILEUP for the default frame time
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call PWRITE( '% Pileup estimate for ACIS:' )
            call ACIS_NOCHOICE( results( 0 ) )
          end if
        else if( special .eq. 504 .or. special .eq. 505 .or.
     &           special .eq. 506 .or. special .eq. 510 .or.
     &           special .eq. 520 ) then
*         Hard-wired code 504: we're dealing with Chandra HETG-ACIS-S HEG1
*         Hard-wired code 505: we're dealing with Chandra HETG-ACIS-S HEG1MEG1
*         Hard-wired code 506: we're dealing with Chandra HETG-ACIS-S MEG1
*         Hard-wired code 510: we're dealing with Chandra LETG-ACIS-S LETG1
*         Hard-wired code 520: we're dealing with Chandra LETG-HRC-S LETG1
          call PWRITE(
     &        '  (=Sum of source count rates in +1st and -1st orders)' )

        else if( special .eq. 601 ) then
*         Hard-wired code 601 means we're dealing with XMM EPIC MOS
*         Pile-up calculations are assumed to be independent of the
*           choice of filters (Open, Med, Thin or Thick)
          call PWRITE( '  for on-axis observation, PATTERN=0-12, ' //
     &                                   'before dead time correction' )
          call PWRITE(
     &         '  and assuming a 5 arcmin extraction radius ' //
     &                                    'enclosing ~100% of the PSF' )
          call PWRITE(
     &         '  (The count rate within a more reasonable point ' //
     &                                      'source extraction region' )
          call PWRITE( '   of ~15 arcsec radius would be roughly ' //
     &                                           '~70% of this value)' )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call EPIC_MOS( results, n_res )
            call PWRITE( 'WARNING: MOS1 Timing mode data now ' //
     &                   'suffers from an out-of-scale column' )
            call PWRITE( 'which results in 20-30% loss of ' //
     &                   'counts for sources at the nominal position.' )
          end if
c        else if( special .eq. 609 ) then
c*         Hard-wired code 609 means we're dealing with XMM EPIC MOS (OPEN)
c          if( restr ) then
c            call PWRITE( '% No instrument specific information ' //
c     &                   'available when energy range is specified.' )
c          else
c            call EPIC_MOS_OPEN( results )
c          end if
c        else if( special .eq. 611 ) then
c*         Hard-wired code 611 means we're dealing with XMM EPIC-PN OPEN
c          if( restr ) then
c            call PWRITE( '% No instrument specific information ' //
c     &                   'available when energy range is specified.' )
c          else
c            call EPIC_PN_OPEN( results, special, restr )
c          end if
        else if( special .ge. 612 .and. special .le. 614 ) then
*         Hard-wired codes 612-614 mean we're dealing with XMM EPIC-PN
*         612: thin, 613: medium, 614: thick
          call PWRITE( '  for on-axis observation, PATTERN=0-4, ' //
     &                                   'before dead time correction' )
          call PWRITE(
     &         '  and assuming a 5 arcmin extraction radius ' //
     &                                    'enclosing ~100% of the PSF' )
          call PWRITE(
     &         '  (The count rate within a more reasonable point ' //
     &                                      'source extraction region' )
          call PWRITE( '   of ~15 arcsec radius would be roughly ' //
     &                                           '~70% of this value)' )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call PWRITE(
     &                  'WARNING: The pile-up estimate is approximate' )
            call EPIC_PN( results, n_res, special )
          end if

        else if( special .eq. 701 ) then
*         Hard-wired code 701 means we're dealing with Suzaku XRS (no filter)
          call PWRITE( '% This estimate is appropriate ' //
     &                 'for a point source centerd on the array' )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call PWRITE( '% Event grade split estimates are:' )
            call XRS_PILEUP( results( 0 ), special )
          end if
        else if( special .eq. 702 ) then
*         Hard-wired code 702 means we're dealing with Suzaku XRS with BE
          call PWRITE( '% This estimate is appropriate ' //
     &                 'for a point source centerd on the array' )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call PWRITE( '% Event grade split estimates are:' )
            call XRS_PILEUP( results( 0 ), special )
          end if

        else if( special .eq. 711 ) then
*         Hard-wired code 711 means we're dealing with Suzaku XIS-FI
          call PWRITE( '  (count rate per FI XIS)' )
          call PWRITE( '  (for a point source at the "XIS nominal"' //
     &                        ' pointing position)' )
          call PWRITE( '  (correction for the 5.6% dead area' //
     &                       ' created by charge injection included)' )

        else if( special .eq. 712 ) then
*         Hard-wired code 712 means we're dealing with Suzaku XIS-BI
          call PWRITE( '  (for a point source at the "XIS nominal"' //
     &                        ' pointing position)' )
          call PWRITE( '  (correction for the 7.4% dead area' //
     &                       ' created by charge injection included)' )

        else if( special .eq. 721 ) then
*         Hard-wired code 721 means we're dealing with Suzaku HXD PIN
          call PWRITE( '  (for a point source at the "XIS nominal"' //
     &                        ' pointing position)' )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call PIN_LIMIT( results, n_res )
          end if
        else if( special .eq. 722 ) then
*         Hard-wired code 722 means we're dealing with Suzaku HXD GSO
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call GSO_LIMIT( results, n_res )
          end if

c
        else if( special .eq. 801 ) then
*         Hard-wired code 801 is reseved for Integral ISGRI
*                 (used to have 802 and 811, but these have been discontinued)
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call INT_LIMIT( special, results, n_res )
          end if
        else if( special .eq. 821 ) then
*         Hard-wired code 821 is for Integral JEM-X
          call PWRITE( '  (Note that 2 units of JEM-X have been' //
     &                        ' operational since October 2010.)' )
          call PWRITE( '  (The abovecount rate is per unit' //
     &                       ' and the two units are very similar.)' )

        else if( special .eq. 901 ) then
*         Hard-wired code 901 means we are dealing with Swift BAT
          call PWRITE( '  (Channel 3-78 count rate per ' //
     &                                   'fully illuminated detector)' )
          if( restr ) then
            call PWRITE( '% Further instrument-specific information ' //
     &                      'can be obtained by re-running PIMMS' )
            call PWRITE( '  without specifying a limited ' //
     &                                    '(non-default) energy range' )
          else
            call BAT_LIMIT( results, n_res )
          end if

        else if( special .eq. 911 ) then
*         Hard-wired code 911: Swift XRT pc (photon counting mode)
          call PWRITE( '  (Grade 0-12 on-axis count ' //
     &                      'rate for an infinite extraction region)'  )
          pw_strng = '   or            cps in Grade 0'
          write( pw_strng( 7: 16 ), '(1p,e10.3)' ) results( 1 )
          call PWRITE( pw_strng )
        else if ( special .eq. 912 ) then
*         Hard-wired code 911: Swift XRT wt (window timing mode)
          call PWRITE( '  (Greade 0-2 on-axis count ' //
     &                      'rate for an infinite extraction region)'  )
          pw_strng = '   or            cps in Grade 0'
          write( pw_strng( 7: 16 ), '(1p,e10.3)' ) results( 1 )
          call PWRITE( pw_strng )
        else if( special .eq. 913 ) then
*         Hard-wired code 911: Swift XRT pd (photodiode mode)
          call PWRITE( '  (Grade 0-5 on-axis count ' //
     &                      'rate for an infinite extraction region)'  )
          pw_strng = '   or            cps in Grade 0-2'
          write( pw_strng( 7: 16 ), '(1p,e10.3)' ) results( 1 )
          call PWRITE( pw_strng )
          pw_strng = '   or            cps in Grade 0'
          write( pw_strng( 7: 16 ), '(1p,e10.3)' ) results( 2 )
          call PWRITE( pw_strng )

        else if( special .eq. 1001 ) then
*          Hard-wired code 1001: NuSTAR
          call PWRITE( '  (For both modules for a 50% PSF ' //
     &                  'extraction without accounting for deadtime)'  )
          if( restr ) then
            call PWRITE( '% PIMMS provides NuSTAR background and ' //
     &                      'deadtime information if run without' )
            call PWRITE( '  specifying a limited output energy range' )
          else
            call NUSTAR_LIMIT( results, n_res )
          end if

        else if( special .ge. 1101 .and. special .le. 1104 ) then
*         Hard-wired code for ASTRO-H SXS (1101: open, 1102: Be25 filter;
*                                          1103: Be50 filter; 1104: OBF)
          call PWRITE( '  (For a point source ' //
     &                                'excatly centered on the array)' )
          call SXS_GRADES( results( 0 ) )

        else if( special .eq. 1131 ) then
*         Hard-wired call for ASTRO-H SGD
          if( restr ) then
            call PWRITE( '% PIMMS provides ASTRO-H SGD background ' //
     &                        'and exposure time for source detection' )
            call PWRITE(
     &     '  if run without specifying a limited output energy range' )
          else
            call SGD_LIMIT( results, n_res )
          end if

        end if

        end
