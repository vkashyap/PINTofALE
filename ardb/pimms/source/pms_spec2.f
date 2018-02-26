*+PMS_SPEC2
        subroutine PMS_SPEC2( mission, detector, filter )

        implicit none

        include 'pms_index.inc'

        integer mission, detector, filter

*       Description:
*         Gives the opportunity to print out mission specific information
*         when they are used as inputs.
*
*       Arguments:
*         mission   (i) : mission number
*         detector  (i) : detector number
*         filter    (i) : filter number
*
*       Dependencies:
*         PIMMS Index routines
*
*       Origin:
*         Special warning when XTE PCA is used as input.
*         Since expanded for many other missions.
*
*       Author:
*         Koji Mukai 1999 July, original version
*-PMS_SPEC2

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
            call PMS_HDWD2( spec_num, temp_unit )
            close( temp_unit )
          end if
        end if

        end

*+PMS_HDWD2
        subroutine PMS_HDWD2( special, s_lun )

        implicit none

        integer special
        integer s_lun

*       Output the special warning for XTE PCA when used as input
*       Modified 2005 September to output instrument specific messages
*         for many missions
*       Minor mode 2005 November for Suzaku
*       Updated 2006 April for Swift/BAT
*       Minor update for XMM EPIC-pn change to PATTERN==0, 2007 September
*       Updated 2010 Jun not to output funny characters when there is
*         no special message to output
*       Updated 2011 Feb to include Integral JEM-X, which now has
*         operational units.
*       Last updated in 2011 Apr to include NuSTAR and ASTRO-H SXS

*       Current "SPECIAL"s (mirrored from PMS_SPECL)
*         101 & 102:	ASCA SIS & GIS
*	  201 - 215:	XTE (PCA, old HEXTE, ASM; XTE with LLDs)
*         304:          HEAO-1 A4
*         401 - 404:    SAX (LECS, MECS, HPGSPC and PDS)
*         501 - 531:    CHANDRA (ACIS pileup)
*         601, 612-614:	XMM (MOS and PN pile-up)
*         701-721:      ASTRO-E/Astro-E2/Suzaku
*         801 & 821:    Integral
*         901-913:      Swift
*         1001:         NuSTAR
*         1101-1104:    ASTRO-H SXS

        character*78 pw_strng
        integer lpw
        integer LENTRIM

        pw_strng = ' '
        if( special .eq. 101 ) then
*         ASCA SIS
          pw_strng = '%!% Integration over the entire chip '
     &                      // '(not just in the source region) assumed'
        else if (special .eq. 102 ) then
*         ASCA GIS
          pw_strng = '%!% Integration over full GIS FOV '
     &                      // '(not just in the source region) assumed'

        else if( special .eq. 201 ) then
*         This hard-wired code means we're dealing with XTE PCA
          pw_strng = '%!% Count rate is assumed to be per PCU'
        else if( special .ge. 210 .and. special .le. 215 ) then
*         This hard-wired code means we're dealing with XTE HEXTE
          pw_strng = '%!% Net count rate per cluster assumed'

        else if ( special .eq. 402 ) then
*         This hard-wired code means we're dealing with SAX MECS
          pw_strng = '%!% Count rate assumed to be for 2 MECS'

        else if( special .eq. 501 .or. special .eq. 502 .or.
     &           special .eq. 507 .or. special .eq. 511 .or.
     &           special .eq. 527 .or. special .eq. 531 ) then
*         501 for ACIS-I, 502 for ACIS-S-BI; others are HETG/LETG order0
          pw_strng = '%!% No pile-up correction will be applied'
        else if( special .eq. 504 .or. special .eq. 505 .or.
     &           special .eq. 506 .or. special .eq. 510 .or.
     &           special .eq. 520 ) then
*         Hard-wired code 504: we're dealing with Chandra HETG-ACIS-S HEG1
*         Hard-wired code 505: we're dealing with Chandra HETG-ACIS-S HEG1MEG1
*         Hard-wired code 506: we're dealing with Chandra HETG-ACIS-S MEG1
*         Hard-wired code 510: we're dealing with Chandra LETG-ACIS-S LETG1
*         Hard-wired code 520: we're dealing with Chandra LETG-HRC-S LETG1
          pw_strng = '%!% Count rate assumed to be sum of +/-1st orders'

        else if( special .eq. 601 ) then
*         Hard-wired code 601 means we're dealing with XMM EPIC MOS
*         Pile-up calculations are assumed to be independent of the
*           choice of filters (Open, Med, Thin or Thick)
          pw_strng = '%!% Pile-up corrected PATTERN=0-12 rate ' //
     &                                     'in 5 arcmin region assumed'
        else if( special .ge. 612 .and. special .le. 614 ) then
*         Hard-wired codes 612-614 mean we're dealing with XMM EPIC-PN
*         612: thin, 613: medium, 614: thick
          pw_strng = '%!% Pile-up corrected PATTERN=0-4 rate ' //
     &                                     'in 5 arcmin region assumed'
*         Will probably need Suzaku special, but turn it off for now

        else if( special .eq. 711 ) then
*         Hard-wired code 711 means we're dealing with Suzaku XIS FI
          pw_strng = '%!% Count rate per 1 FI unit assumed'

        else if( special .eq. 721 ) then
*         Hard-wired code 721 means we're dealing with Suzaku HXD PIN
          pw_strng =
     &     '%!% Observation at "XIS Nominal" pointing position assumed'

        else if( special .eq. 821 ) then
*         Hard-wired code 821 means we're dealing with Integral JEM-X
          pw_strng = '%!% Count rate per unit assumed ' //
     &              '(Note that 2 units are operational since 2010 Oct)'

        else if ( special .eq. 901 ) then
*         Hard-wired code 901 means we are dealing with Swift BAT
          pw_strng = '%!% Source count rate (channels 3-78) '
     &                       // 'per fully illuminated detector assumed'

        else if ( special .eq. 911 ) then
*         Hard-wired code 911: Swift XRT pc (photon counting mode)
          pw_strng =
     &'%!% Total (not just source region) Grade 0-12 count rate assumed'
        else if ( special .eq. 912 ) then
*         Hard-wired code 911: Swift XRT wt (window timing mode)
          pw_strng =
     & '%!% Total (not just source region) Grade 0-2 count rate assumed'
        else if ( special .eq. 913 ) then
*         Hard-wired code 911: Swift XRT pd (photodiode mode)
          pw_strng =
     & '%!% Total (not just source region) Grade 0-5 count rate assumed'

        else if ( special .eq. 1001 ) then
*         Hard-wired code 1001: NuSTAR
          pw_strng = '%!% Dead-time corrected 2-module rate for ' //
     &                               '50% PSF extraction region assumed'

        else if ( special .ge. 1101 .and. special .le. 1104 ) then
*         Hard-wired code for ASTRO-H SXS (1101: open, 1102: Be25 filter;
*                                          1103: Be50 filter; 1104: OBF)
          pw_strng = '%!% Count rate assumed to be for ' //
     &                    'a point source excatly centered on the array'

        end if

        lpw = LENTRIM( pw_strng )
        if( lpw .gt. 1 ) then
          call PWRITE( pw_strng( : lpw ) )
        end if

        end
