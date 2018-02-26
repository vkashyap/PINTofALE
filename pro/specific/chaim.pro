function chaim,ra,dec,ChGEOM=ChGEOM,ChOCAT=ChOCAT,sobs=sobs,$
	verbose=verbose, _extra=e
;+
;function	chaim
;	find and return Chandra ObsIDs that contain observations
;	of list of RA,Dec
;
;syntax
;	chmatch=chaim(ra,dec,ChGEOM=ChGEOM,ChOCAT=ChOCAT,sobs=sobs,$
;	verbose=verbose, caldb=caldb,version=version)
;
;parameters
;	RA	[INPUT; required] RA_2000.0 in [deg]
;	Dec	[INPUT; required] Dec_2000.0 in [deg]
;		* size of DEC must match that of RA
;
;keywords
;	ChGEOM	[I/O] Chandra geometry file
;		* read in from CALDB location if not present on input
;	ChOCAT	[I/O] name of file containing the Chandra Observation CATalog,
;		or output of RDB structure containing the same -- contains list
;		of Chandra ObsIDs and and attendant information such as RA_NOM,
;		Dec_NOM, ROLL_NOM, and which chips are on, etc.
;		* if valid structure, will use it as is; otherwise:
;		-- if filename, will read from that file and return
;		   the output structure in this keyword
;		* if not specified, will look to read from file 'axafocat.rdb'
;	sobs	[INPUT] an integer list of ObsIDs to consider as a
;		subset of ChOCAT
;	verbose	[INPUT] controls chatter
;		* if .GE.5, makes a plot displaying the locations
;		  of the matches on the detector
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		CALDB [rd_chandra_geom]
;		VERSION [rd_chandra_geom]
;
;warning
;	spacecraft roll is approximated in Euclidean space.  this will
;	result in some deviations from true, especially for those
;	observations involving the grating arrays and aligned N-S.
;	however, this will be an issue only for the extreme chips edges
;	and is unlikely to have a practical effect.
;
;subroutines
;	RD_CHANDRA_GEOM
;	RDB
;
;history
;	vinay kashyap (Jul2005)
;	bug correction: roll off by -90 deg (VK; Aug2005)
;-

;	usage
ok='ok' & np=n_params() & nr=n_elements(ra) & nd=n_elements(dec)
sra=size(ra) & nsra=n_elements(sra)
sdec=size(dec) & nsdec=n_elements(sdec)
if np lt 2 then ok='Insufficient parameters' else $
 if nr eq 0 then ok='RA: undefined' else $
  if nd eq 0 then ok='Dec: undefined' else $
   if nr ne nd then ok='RA and DEC are incompatible' else $
    if sra[nsra-2] eq 7 then ok='RA must be numeric' else $
     if sdec[nsdec-2] eq 7 then ok='DEC must be numeric' else $
      if max(ra) gt 360 then ok='RA must be in [degrees]' else $
       if min(ra) lt 0 then ok='RA is negative?' else $
       if max(abs(dec)) gt 90 then ok='DEC must be in [degrees]'
if ok ne 'ok' then begin
  print,'Usage: chmatch=chaim(ra,dec,ChGEOM=ChGEOM,ChOCAT=ChOCAT,sobs=sobs,$'
  print,'       verbose=verbose, caldb=caldb,version=version)'
  print,'  find and return Chandra ObsIDs that have observations of RA,DEC'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	read in Chandra geometry if necessary
if n_elements(ChGEOM) eq 0 then ok='ChGEOM not defined'
ng=n_tags(ChGEOM) & if ng eq 0 then ok='ChGEOM not a structure'
if ng gt 0 then begin
  gnam=tag_names(ChGEOM)
  ;	check here for the existence of all relevant extensions
  ;	after this, correct structure will be assumed to exist
  inam=0L
  for i=0L,ng-1L do begin
    if strpos(gnam[i],'INSTRUMENTS',0) ge 0 then inam=inam+1L
    if strpos(gnam[i],'GRATINGS',0) ge 0 then inam=inam+1L
    if strpos(gnam[i],'INSTGEOM1',0) ge 0 then inam=inam+1L
    if strpos(gnam[i],'INSTGEOM2',0) ge 0 then inam=inam+1L
    if strpos(gnam[i],'INSTGEOM3',0) ge 0 then inam=inam+1L
    if strpos(gnam[i],'INSTGEOM4',0) ge 0 then inam=inam+1L
    if strpos(gnam[i],'INSTGEOM5',0) ge 0 then inam=inam+1L
  endfor
  if inam lt 7 then ok='ChGEOM does not have correct format'
endif
if ok ne 'ok' then begin
  if n_elements(ChGEOM) ne 0 then message,ok,/informational
  ChGEOM=rd_chandra_geom(_extra=e)
  if n_tags(ChGEOM) eq 0 then begin
    message,'ERROR: ChGEOM cannot be found',/informational
    return,-1L
  endif
endif

;	read in Chandra ObsID list if necessary
nc=n_elements(ChOCAT) & mc=n_tags(ChOCAT)
if nc eq 0 then ChOCAT='axafocat.rdb'
if mc eq 0 then begin		;(ChOCAT not a structure
  szc=size(ChOCAT) & nszc=n_elements(szc)
  if szc[nszc-2] eq 7 then begin
    fil=findfile(ChOCAT[0],count=nfil)
    if nfil eq 0 then begin
      message,'ERROR: '+ChOCAT[0]+': file not found',/informational
      return,-1L
    endif
    tmp=rdb(ChOCAT[0])
    ChOCAT=tmp		;overwrite!
  endif
endif				;MC=0)
mc=n_tags(ChOCAT) & mnam=tag_names(ChOCAT) & inam=0L
for i=0L,mc-1L do begin
  if strpos(mnam[i],'OBSID',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'INSTRUMENT',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'GRATING',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'SI_MODE',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'RA',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'DEC',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'SOE_ROLL',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'Y_DET_OFFSET',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'Z_DET_OFFSET',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDI0_ON',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDI1_ON',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDI2_ON',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDI3_ON',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDS0_ON',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDS1_ON',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDS2_ON',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDS3_ON',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDS4_ON',0) ge 0 then inam=inam+1L
  if strpos(mnam[i],'CCDS5_ON',0) ge 0 then inam=inam+1L
endfor
if inam lt 19 then begin
  message,'ERROR: ChOCAT is not in the right format',/informational
  return,-1L
endif

obsid=long(ChOCAT.OBSID)
instrument=ChOCAT.INSTRUMENT
grating=ChOCAT.GRATING
si_mode=ChOCAT.SI_MODE
status=ChOCAT.STATUS
chra=ChOCAT.RA
  oo=where(strtrim(chra,2) eq 'NULL',moo)
  if moo gt 0 then chra[oo]='0' & chra=float(chra)
chdec=ChOCAT.DEC
  oo=where(strtrim(chdec,2) eq 'NULL',moo)
  if moo gt 0 then chdec[oo]='0' & chdec=float(chdec)
chroll=ChOCAT.SOE_ROLL
  oo=where(strtrim(chroll,2) eq 'NULL',moo)
  if moo gt 0 then chroll[oo]='0' & chroll=float(chroll)
y_offset=ChOCAT.Y_DET_OFFSET
  ;oo=where(strtrim(y_offset,2) eq 'NULL',moo)
  ;if moo gt 0 then y_offset[oo]='0' & y_offset=float(y_offset)
z_offset=ChOCAT.Z_DET_OFFSET
  ;oo=where(strtrim(z_offset,2) eq 'NULL',moo)
  ;if moo gt 0 then z_offset[oo]='0' & z_offset=float(z_offset)
trans_offset=ChOCAT.TRANS_OFFSET
  ;oo=where(strtrim(trans_offset,2) eq 'NULL',moo)
  ;if moo gt 0 then trans_offset[oo]='0' & trans_offset=float(trans_offset)
subarray=ChOCAT.SUBARRAY
subarray_row=ChOCAT.SUBARRAY_START_ROW
  oo=where(strtrim(subarray_row,2) eq 'NULL',moo)
  if moo gt 0 then subarray_row[oo]='1' & subarray_row=fix(subarray_row)
ccdi0_on=ChOCAT.CCDI0_ON & ccdi1_on=ChOCAT.CCDI1_ON
ccdi2_on=ChOCAT.CCDI2_ON & ccdi3_on=ChOCAT.CCDI3_ON
ccds0_on=ChOCAT.CCDS0_ON & ccds1_on=ChOCAT.CCDS1_ON
ccds2_on=ChOCAT.CCDS2_ON & ccds3_on=ChOCAT.CCDS3_ON
ccds4_on=ChOCAT.CCDS4_ON & ccds5_on=ChOCAT.CCDS5_ON
approved_exposure_time=ChOCAT.APPROVED_EXPOSURE_TIME
rem_exp_time=ChOCAT.REM_EXP_TIME
soe_st_sched_date=ChOCAT.SOE_ST_SCHED_DATE
ncat=n_elements(obsid)
icat=lonarr(ncat)-1L	;will only look at ObsIDs where ICAT.GE.0
rcat=lindgen(ncat)+1	;this keeps track of repeated observations

;	cue repeat observations, those with the exact same
;	RA, DEC, ROLL, INSTRUMENT, GRATING, SI_MODE, OFFSETS,
;	SUBARRAY, and CCDI?_ON
;for i=0L,ncat-2L do begin
;  if 100*long(i/100) eq i then kilroy
;  if instrument[i] eq 'ACIS-I' or instrument[i] eq 'ACIS-S' then begin
;    oo=where(rcat[i+1L:*] gt 0 and $
;	instrument[i+1L:*] eq instrument[i] and $
;	grating[i+1L:*] eq grating[i] and $
;	si_mode[i+1L:*] eq si_mode[i] and $
;	abs(chra[i+1L:*]-chra[i]) lt 1e-4 and $
;	abs(chdec[i+1L:*]-chdec[i]) lt 1e-4 and $
;	abs(chroll[i+1L:*]-chroll[i]) lt 0.1 and $
;	y_offset[i+1L:*] eq y_offset[i] and $
;	z_offset[i+1L:*] eq z_offset[i] and $
;	trans_offset[i+1L:*] eq trans_offset[i] and $
;	subarray[i+1L:*] eq subarray[i] and $
;	subarray_row[i+1L:*] eq subarray_row[i] and $
;	ccdi0_on[i+1L:*] eq ccdi0_on[i] and $
;	ccdi1_on[i+1L:*] eq ccdi1_on[i] and $
;	ccdi2_on[i+1L:*] eq ccdi2_on[i] and $
;	ccdi3_on[i+1L:*] eq ccdi3_on[i] and $
;	ccds0_on[i+1L:*] eq ccds0_on[i] and $
;	ccds1_on[i+1L:*] eq ccds1_on[i] and $
;	ccds2_on[i+1L:*] eq ccds2_on[i] and $
;	ccds3_on[i+1L:*] eq ccds3_on[i] and $
;	ccds4_on[i+1L:*] eq ccds4_on[i] and $
;	ccds5_on[i+1L:*] eq ccds5_on[i],$
;	moo)
;    if moo gt 0 then rcat[oo+i+1L]=-i-1L
;  endif
;endfor

;	check keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
nobs=n_elements(sobs)
if nobs eq 0 then $				;... look at all ObsIDs
	icat=lindgen(ncat) else begin		;(look only at a subset
  for i=0L,nobs-1L do begin
    oo=where(long(obsid) eq sobs[i],moo)
    if moo gt 0 then icat[oo]=oo
  endfor
endelse						;NOBS)

;	output
chmatch=strarr(nr)

;	nominal aimpoint locations [mm]
acismmpix = 0.023957031	;[mm/pix]
acisdegmm = 0.492/acismmpix/3600. ;[deg/mm]=[arcsec/pix]*[pix/mm]*[deg/arcsec]
hrcmmpix = 0.00642938	;[mm/pix]
hrcdegmm = 0.13175/3600./hrcmmpix ;[deg/mm]=[arcsec/pix]*[deg/arcsec]*[pix/mm]
nom_acisi_y = 962.*acismmpix & nom_acisi_z = 964.*acismmpix	;[mm] on I3
nom_aciss_y = 252.*acismmpix & nom_aciss_z = 510.*acismmpix	;[mm] on S3
nom_hrci_y = 0. & nom_hrci_z = 0.				;[mm] on I
nom_hrcs_y = -4. & nom_hrcs_z = 0.				;[mm] on S2

;	go through each row of ChOCAT, construct the view, and
;	see how many of the RA,Dec overlap the FOV
ok=where(icat ge 0,mok)
for i=0L,mok-1L do begin
  j=ok[i] & ramin=0. & ramax=0. & decmin=0. & decmax=0.
  tht=chroll[j]*!pi/180. - !pi/2.
  yoff=0. & zoff=0.
  if instrument[j] eq 'ACIS-S' and grating[j] eq 'LETG' then $
	yoff=(-0.33/60.)/acisdegmm	;[mm] = [deg]/[deg/mm]

  if strpos(instrument[j],'ACIS',0) ge 0 then begin		;(ACIS

    NCHIP=10
    if y_offset[j] ne 'NULL' then yoff=float(y_offset[j])/60./acisdegmm
    if z_offset[j] ne 'NULL' then zoff=float(z_offset[j])/60./acisdegmm

    if instrument[j] eq 'ACIS-I' then begin
      ll=chgeom.INSTGEOM1.LL & lr=chgeom.INSTGEOM1.LR
      ur=chgeom.INSTGEOM1.UR & ul=chgeom.INSTGEOM1.UL
      ori_y=nom_acisi_z-ll[1,3] & ori_z=ll[2,3]+nom_acisi_y
      ori_y=ori_y+yoff & ori_z=ori_z+zoff
      ll[1,*]=ll[1,*]-ori_y & ll[2,*]=ll[2,*]-ori_z
      lr[1,*]=lr[1,*]-ori_y & lr[2,*]=lr[2,*]-ori_z
      ur[1,*]=ur[1,*]-ori_y & ur[2,*]=ur[2,*]-ori_z
      ul[1,*]=ul[1,*]-ori_y & ul[2,*]=ul[2,*]-ori_z
    endif
    if instrument[j] eq 'ACIS-S' then begin
      ll=chgeom.INSTGEOM1.LL & lr=chgeom.INSTGEOM1.LR
      ur=chgeom.INSTGEOM1.UR & ul=chgeom.INSTGEOM1.UL
      ori_y=ll[1,7]+nom_aciss_y & ori_z=ll[2,7]+nom_aciss_z
      ori_y=ori_y+yoff & ori_z=ori_z+zoff
      ll[1,*]=ll[1,*]-ori_y & ll[2,*]=ll[2,*]-ori_z
      lr[1,*]=lr[1,*]-ori_y & lr[2,*]=lr[2,*]-ori_z
      ur[1,*]=ur[1,*]-ori_y & ur[2,*]=ur[2,*]-ori_z
      ul[1,*]=ul[1,*]-ori_y & ul[2,*]=ul[2,*]-ori_z
    endif

    k=0
    yreg_i0=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_i0=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_i0=yreg_i0*sin(tht)+zreg_i0*cos(tht) + chra[j]
    dec_reg_i0=-yreg_i0*cos(tht)+zreg_i0*sin(tht) + chdec[j]
    k=1
    yreg_i1=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_i1=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_i1=yreg_i1*sin(tht)+zreg_i1*cos(tht) + chra[j]
    dec_reg_i1=-yreg_i1*cos(tht)+zreg_i1*sin(tht) + chdec[j]
    k=2
    yreg_i2=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_i2=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_i2=yreg_i2*sin(tht)+zreg_i2*cos(tht) + chra[j]
    dec_reg_i2=-yreg_i2*cos(tht)+zreg_i2*sin(tht) + chdec[j]
    k=3
    yreg_i3=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_i3=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_i3=yreg_i3*sin(tht)+zreg_i3*cos(tht) + chra[j]
    dec_reg_i3=-yreg_i3*cos(tht)+zreg_i3*sin(tht) + chdec[j]
    k=4
    yreg_s0=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_s0=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_s0=yreg_s0*sin(tht)+zreg_s0*cos(tht) + chra[j]
    dec_reg_s0=-yreg_s0*cos(tht)+zreg_s0*sin(tht) + chdec[j]
    k=5
    yreg_s1=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_s1=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_s1=yreg_s1*sin(tht)+zreg_s1*cos(tht) + chra[j]
    dec_reg_s1=-yreg_s1*cos(tht)+zreg_s1*sin(tht) + chdec[j]
    k=6
    yreg_s2=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_s2=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_s2=yreg_s2*sin(tht)+zreg_s2*cos(tht) + chra[j]
    dec_reg_s2=-yreg_s2*cos(tht)+zreg_s2*sin(tht) + chdec[j]
    k=7
    yreg_s3=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_s3=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_s3=yreg_s3*sin(tht)+zreg_s3*cos(tht) + chra[j]
    dec_reg_s3=-yreg_s3*cos(tht)+zreg_s3*sin(tht) + chdec[j]
    k=8
    yreg_s4=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_s4=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_s4=yreg_s4*sin(tht)+zreg_s4*cos(tht) + chra[j]
    dec_reg_s4=-yreg_s4*cos(tht)+zreg_s4*sin(tht) + chdec[j]
    k=9
    yreg_s5=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*acisdegmm
    zreg_s5=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*acisdegmm
    ra_reg_s5=yreg_s5*sin(tht)+zreg_s5*cos(tht) + chra[j]
    dec_reg_s5=-yreg_s5*cos(tht)+zreg_s5*sin(tht) + chdec[j]
    ;
    ramin=min([ra_reg_i0,ra_reg_i1,ra_reg_i2,ra_reg_i3,$
      ra_reg_s0,ra_reg_s1,ra_reg_s2,ra_reg_s3,ra_reg_s4,ra_reg_s5],$
      max=ramax)
    decmin=min([dec_reg_i0,dec_reg_i1,dec_reg_i2,dec_reg_i3,$
      dec_reg_s0,dec_reg_s1,dec_reg_s2,dec_reg_s3,dec_reg_s4,dec_reg_s5],$
      max=decmax)

  endif							;ACIS)

  if strpos(instrument[j],'HRC',0) ge 0 then begin	;(HRC

    if y_offset[j] ne 'NULL' then yoff=float(y_offset[j])/60./hrcdegmm
    if z_offset[j] ne 'NULL' then zoff=float(z_offset[j])/60./hrcdegmm

    if instrument[j] eq 'HRC-I' then begin
      NCHIP=1
      ll=chgeom.INSTGEOM2.LL & lr=chgeom.INSTGEOM2.LR
      ur=chgeom.INSTGEOM2.UR & ul=chgeom.INSTGEOM2.UL
      ori_y=nom_hrci_z-ll[1] & ori_z=ll[2]+nom_hrci_y
      ori_y=ori_y+yoff & ori_z=ori_z+zoff
      ll[1]=ll[1]-ori_y & ll[2]=ll[2]-ori_z
      lr[1]=lr[1]-ori_y & lr[2]=lr[2]-ori_z
      ur[1]=ur[1]-ori_y & ur[2]=ur[2]-ori_z
      ul[1]=ul[1]-ori_y & ul[2]=ul[2]-ori_z

      yreg_i0=[ll[1],lr[1],ur[1],ul[1],ll[1]]*hrcdegmm
      zreg_i0=[ll[2],lr[2],ur[2],ul[2],ll[2]]*hrcdegmm
      ra_reg_i0=yreg_i0*sin(tht)+zreg_i0*cos(tht) + chra[j]
      dec_reg_i0=-yreg_i0*cos(tht)+zreg_i0*sin(tht) + chdec[j]

      ramin=min(ra_reg_i0,max=ramax)
      decmin=min(dec_reg_i0,max=decmax)
    endif
    if instrument[j] eq 'HRC-S' then begin
      NCHIP=3
      ll=chgeom.INSTGEOM3.LL & lr=chgeom.INSTGEOM3.LR
      ur=chgeom.INSTGEOM3.UR & ul=chgeom.INSTGEOM3.UL
      ori_y=ll[1,1]+nom_hrcs_y & ori_z=ll[2,1]+nom_hrcs_z
      ori_y=ori_y+yoff & ori_z=ori_z+zoff
      ll[1,*]=ll[1,*]-ori_y & ll[2,*]=ll[2,*]-ori_z
      lr[1,*]=lr[1,*]-ori_y & lr[2,*]=lr[2,*]-ori_z
      ur[1,*]=ur[1,*]-ori_y & ur[2,*]=ur[2,*]-ori_z
      ul[1,*]=ul[1,*]-ori_y & ul[2,*]=ul[2,*]-ori_z

      k=0
      yreg_i1=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*hrcdegmm
      zreg_i1=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*hrcdegmm
      ra_reg_i1=yreg_i1*sin(tht)+zreg_i1*cos(tht) + chra[j]
      dec_reg_i1=-yreg_i1*cos(tht)+zreg_i1*sin(tht) + chdec[j]
      k=1
      yreg_i2=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*hrcdegmm
      zreg_i2=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*hrcdegmm
      ra_reg_i2=yreg_i2*sin(tht)+zreg_i2*cos(tht) + chra[j]
      dec_reg_i2=-yreg_i2*cos(tht)+zreg_i2*sin(tht) + chdec[j]
      k=2
      yreg_i3=[ll[1,k],lr[1,k],ur[1,k],ul[1,k],ll[1,k]]*hrcdegmm
      zreg_i3=[ll[2,k],lr[2,k],ur[2,k],ul[2,k],ll[2,k]]*hrcdegmm
      ra_reg_i3=yreg_i3*sin(tht)+zreg_i3*cos(tht) + chra[j]
      dec_reg_i3=-yreg_i3*cos(tht)+zreg_i3*sin(tht) + chdec[j]

      ramin=min([ra_reg_i1,ra_reg_i2,ra_reg_i3],max=ramax)
      decmin=min([dec_reg_i1,dec_reg_i2,dec_reg_i3],max=decmax)
    endif
  endif							;HRC)

    om=where(ra ge ramin and ra le ramax and $
	dec ge decmin and dec le decmax,mom)
    if mom gt 0 then begin			;(any possible matches?
      if vv ge 5 then begin
	cc=strtrim(obsid[j],2)+' : '+instrument[j]
	if grating[j] ne 'NONE' then cc=cc+'/'+grating[j]
        plot,[0],xr=[ramax,ramin],yr=[decmin,decmax],title=cc
	oplot,ra[om],dec[om],psym=4,col=3
	if nchip eq 10 then begin
          oplot,ra_reg_i0,dec_reg_i0,col=2 & oplot,ra_reg_i1,dec_reg_i1,col=2
          oplot,ra_reg_i2,dec_reg_i2,col=2 & oplot,ra_reg_i3,dec_reg_i3,col=2
          oplot,ra_reg_s0,dec_reg_s0,col=2 & oplot,ra_reg_s1,dec_reg_s1,col=2
          oplot,ra_reg_s2,dec_reg_s2,col=2 & oplot,ra_reg_s3,dec_reg_s3,col=2
          oplot,ra_reg_s4,dec_reg_s4,col=2 & oplot,ra_reg_s5,dec_reg_s5,col=2
	endif
	if nchip eq 1 then oplot,ra_reg_i0,dec_reg_i0,col=2
	if nchip eq 3 then begin
          oplot,ra_reg_i1,dec_reg_i1,col=2
          oplot,ra_reg_i2,dec_reg_i2,col=2
	  oplot,ra_reg_i3,dec_reg_i3,col=2
	endif
        plots,0,0,psym=1
      endif
      for ichip=0,nchip-1 do begin		;{check for each CCD
	ccd_on='Y'
	if NCHIP eq 10 then begin
	  case ichip of
	    0: begin & ccd_on=ccdi0_on[j] & x=ra_reg_i0 & y=dec_reg_i0 & end
	    1: begin & ccd_on=ccdi1_on[j] & x=ra_reg_i1 & y=dec_reg_i1 & end
	    2: begin & ccd_on=ccdi2_on[j] & x=ra_reg_i2 & y=dec_reg_i2 & end
	    3: begin & ccd_on=ccdi3_on[j] & x=ra_reg_i3 & y=dec_reg_i3 & end
	    4: begin & ccd_on=ccds0_on[j] & x=ra_reg_s0 & y=dec_reg_s0 & end
	    5: begin & ccd_on=ccds1_on[j] & x=ra_reg_s1 & y=dec_reg_s1 & end
	    6: begin & ccd_on=ccds2_on[j] & x=ra_reg_s2 & y=dec_reg_s2 & end
	    7: begin & ccd_on=ccds3_on[j] & x=ra_reg_s3 & y=dec_reg_s3 & end
	    8: begin & ccd_on=ccds4_on[j] & x=ra_reg_s4 & y=dec_reg_s4 & end
	    9: begin & ccd_on=ccds5_on[j] & x=ra_reg_s5 & y=dec_reg_s5 & end
	  endcase
	endif
	if NCHIP eq 1 then begin
	  x=ra_reg_i0 & y=dec_reg_i0
	endif
	if NCHIP eq 3 then begin
	  case ichip of
	    0: begin
	      if si_mode[j] eq 'SCENTER' or si_mode[j] eq 'SPWING' then $
		ccd_on='N'
	      x=ra_reg_i1 & y=dec_reg_i1
	    end
	    1: begin
	      if si_mode[j] eq 'SNWING' or si_mode[j] eq 'SPWING' then $
		ccd_on='N'
	      x=ra_reg_i2 & y=dec_reg_i2
	    end
	    2: begin
	      if si_mode[j] eq 'SCENTER' or si_mode[j] eq 'SNWING' then $
		ccd_on='N'
	      x=ra_reg_i3 & y=dec_reg_i3
	    end
	  endcase
	endif
	;print,ccd_on,ichip
        if ccd_on eq 'Y' then begin		;(if CCD is on
	  if vv ge 5 then oplot,x,y,thick=3
	  n=n_elements(x)
	  theta=fltarr(mom)
	  for ir=0,n-2 do begin			;{calculate angles
	    x1p=x[ir]-ra[om] & y1p=y[ir]-dec[om]
	    x2p=x[ir+1]-ra[om] & y2p=y[ir+1]-dec[om]
	    dpp=x1p*x2p+y1p*y2p & cpp=x1p*y2p-x2p*y1p
	    theta=theta+atan(cpp,dpp)
	  endfor				;IR=0,N-2}
	  oin=where(abs(theta) gt 2.*!pi-0.01,moin)
	  if moin gt 0 then begin		;(there are matches
	    chmatch[om[oin]]=chmatch[om[oin]]+strtrim(obsid[j],2)+' '
	    if vv ge 5 then oplot,ra[om[oin]],dec[om[oin]],psym=2
	    if vv ge 6 then xyouts,ra[om[oin]],dec[om[oin]],$
		'  '+strtrim(1+om[oin],2)
	    if vv ge 10 then wait,vv/10.
	  endif					;MOIN>0)
        endif					;CCD_ON)
      endfor					;ICHIP=0,NCHIP-1}
    endif					;MOM)

endfor

return,chmatch
end
