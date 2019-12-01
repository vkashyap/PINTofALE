pro wrt_ogip_rmf,rmf,rmfile,hdrfil=hdrfil,$
    matext=matext,ebext=ebext,telescop=telescop,instrume=instrume,$
    filter=filter,detnam=detnam,grating=grating,chantype=chantype,$
    detchans=detchans,hduclass=hduclass,hduclas1=hduclas1,$
    hduclas2_rm=hduclas2_rm,hduclass2_eb=hduclas2_eb,hduclas3=hduclas3,$
    hduvers=hduvers,tlmin4=tlmin4,numgrp=numgrp,numelt=numelt,$
    lo_thres=lo_thres,ccls=ccls,ccnm=ccnm,cdtp=cdtp,cvsd=cvsd,$
    cvst=cvst,cdes=cdes,rmfversn=rmfversn,hduvers1=hduvers1,$
    hduvers2=hduvers2,origin=origin,creator=creator,revision=revision,$
    mission=mission, _extra=e
;+
;procedure	wrt_ogip_rmf
;	write out a response matrix into an OGIP-compatible FITS file
;
;syntax
;	wrt_ogip_rmf,rmf,rmfile,hdrfil=hdrfil,matext=matext,ebext=ebext,$
;	telescop=telescop,instrume=instrume,filter=filter,detnam=detnam,$
;	grating=grating,chantype=chantype,detchans=detchans,$
;	hduclass=hduclass,hduclas1=hduclas1,hduclas2_rm=hduclas2_rm,$
;	hduclas2_eb=hduclas2_eb,hduclas3=hduclas3,hduvers=hduvers,$
;	tlmin4=tlmin4,numgrp=numgrp,numelt=numelt,lo_thres=lo_thres,$
;	/ccls,/ccnm,/cdtp,cvsd=cvsd,cvst=cvst,cdes=cdes,$
;	rmfversn=rmfversn,hduvers1=hduvers1,hduvers2=hduvers2,$
;	origin=origin,creator=creator,revision=revision,mission=mission
;
;parameters
;	rmf	[INPUT; required] response matrix structure, in the
;		same format as is read in by, e.g., RD_OGIP_RMF()
;	rmfile	[INPUT; required] name of output file
;
;keywords
;	hdrfil	[INPUT] name of existing OGIP-style response matrix file
;		from which to copy the headers as far as possible
;
;		the following are header keywords recommended by OGIP,
;		along with a brief description of what they mean.  the
;		default values used by this program are also listed.
;		see http://legacy.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html
;
;	matext		= 'MATRIX' 	/MATRIX or SPECRESP MATRIX, name of extension
;	ebext		= 'EBOUNDS'	/name of extension
;	telescop	= 'TELESCOP'	/telescope used
;	instrume	= 'INSTRUME'	/instrument used
;	filter		= 'NONE'	/the instrument filter in use (if any)
;	detnam		= 'DETNAM'	/detector used
;	grating		= 'GRATING'	/grating used (if any)
;	chantype	= 'PHA' 	/PHA or PI, uncorrected or corrected for gain
;	detchans	= n_elements(rmf.EMN)	/number of raw detector channels
;	hduclass	= 'OGIP'	/file format is OGIP standard
;	hduclas1	= 'RESPONSE'	/extension contains response data
;	hduclas2_rm	= 'RSP_MATRIX'	/extension contains a response matrix
;	hduclas2_eb	= 'EBOUNDS'	/extension contains a response matrix
;	hduvers		= '1.3.0'	/version of the file format
;	tlmin4		= rmf.FIRSTCHAN	/first channel in the response
;	numgrp		= total(rmf.N_GRP)	/number of channel subsets
;	numelt		= total(rmf.N_CHAN)	/number of response elements
;	lo_thres	= 0.000001	/lower threshold used to construct the matrix
;	hduclas3	= 'DETECTOR'	/REDIST (RMF) or DETECTOR (with QE, no ARF) or FULL (everything)
;	rmfversn	= '1992a'	/obsolete
;	hduvers1	= '1.1.0'	/obsolete
;	hduvers2	= '1.2.0'	/obsolete
;	ccls		= 'CPF'		/the OGIP-class of this calibration file
;	ccnm		= 'MATRIX' or 'EBOUNDS'	/CIF codename for this type of calibration dataset
;	cdtp		= 'DATA'	/OGIP code for the form of the contents of the file
;	cvsd		= YYYY-MM-DD	/UTC date of beginning applicability
;	cvst		= HH:MM:SS	/UTC time of beginning applicability
;	cdes		= 'RSP'		/brief descriptive summary of this dataset
;	origin		= whoami@host	/source of FITS file
;	creator		= 'PINTofALE'	/program creating this file
;	revision	= 1
;	mission		= 'MISSION'	/mission
;	
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	* requires the IDL-Astro library
;	* if not a UNIX system, keyword ORIGIN must be explicitly set
;
;side-effects
;	writes out a FITS file to disk
;
;history
;	vinay kashyap (Nov2001)
;	force variable length arrays even when N_GRP=1 (VK; Jul02)
;	if input RMF structure has a different number of rows than the header template,
;	  NAXIS2 is updated accordingly (VK; Jul19)
;-

;	usage
ok='ok' & np=n_params()
nr=n_elements(rmf) & nf=n_elements(rmfile) & nrr=n_tags(rmf)
if np lt 2 then ok='Insufficient parameters' else $
 if nr eq 0 then ok='RMF is undefined' else $
  if nf eq 0 then ok='RMFILE is undefined' else $
   if nrr eq 0 then ok='RMF must be a structure'
if ok ne 'ok' then begin
  print,'Usage: wrt_ogip_rmf,rmf,rmfile,hdrfil=hdrfil,matext=matext,$'
  print,'       ebext=ebext,telescop=telescop,instrume=instrume,$'
  print,'       filter=filter,detnam=detnam,grating=grating,$'
  print,'       chantype=chantype,detchans=detchans,hduclass=hduclass,$'
  print,'       hduclas1=hduclas1,hduclas2_rm=hduclas2_rm,$'
  print,'       hduclas2_eb=hduclas2_eb,hduclas3=hduclas3,hduvers=hduvers,$'
  print,'       tlmin4=tlmin4,numgrp=numgrp,numelt=numelt,lo_thres=lo_thres,$'
  print,'       /ccls,/ccnm,/cdtp,cvsd=cvsd,cvst=cvst,cdes=cdes,$'
  print,'       rmfversn=rmfversn,hduvers1=hduvers1,hduvers2=hduvers2,$
  print,'       origin=origin,creator=creator,revision=revision,mission=mission'
  print,'  write out a response matrix into an OGIP-compatible FITS file'
  if np ne 0 then message,ok,/info
  return
endif

;	verify RMF
t=tag_names(rmf) & k=0
for i=0,nrr-1 do begin
  if t[i] eq 'ELO' then k=k+1
  if t[i] eq 'EHI' then k=k+1
  if t[i] eq 'EMN' then k=k+1
  if t[i] eq 'EMX' then k=k+1
  if t[i] eq 'N_GRP' then k=k+1
  if t[i] eq 'F_CHAN' then k=k+1
  if t[i] eq 'N_CHAN' then k=k+1
  if t[i] eq 'MATRIX' then k=k+1
  if t[i] eq 'FIRSTCHAN' then k=k+1
endfor
if k lt 9 then begin
  message,'RMF not a standard response matrix structure; see RD_OGIP_RMF()',/info
  return
endif

;	make a HISTORY card out of current call
hist='wrt_ogip_rmf,rmf,'+rmfile[0]
if keyword_set(hdrfil) then hist=hist+',hdrfil='+strtrim(hdrfil[0],2)
if keyword_set(matext) then hist=hist+',matext='+strtrim(matext[0],2)
if keyword_set(ebext) then hist=hist+',ebext='+strtrim(ebext[0],2)
if keyword_set(telescop) then hist=hist+',telescop='+strtrim(telescop[0],2)
if keyword_set(instrume) then hist=hist+',instrume='+strtrim(instrume[0],2)
if keyword_set(filter) then hist=hist+',filter='+strtrim(filter[0],2)
if keyword_set(detnam) then hist=hist+',detnam='+strtrim(detnam[0],2)
if keyword_set(grating) then hist=hist+',grating='+strtrim(grating[0],2)
if keyword_set(chantype) then hist=hist+',chantype='+strtrim(chantype[0],2)
if keyword_set(detchans) then hist=hist+',detchans='+strtrim(detchans[0],2)
if keyword_set(hduclass) then hist=hist+',hduclass='+strtrim(hduclass[0],2)
if keyword_set(hduclas1) then hist=hist+',hduclas1='+strtrim(hduclas1[0],2)
if keyword_set(hduclas2_rm) then hist=hist+',hduclas2_rm='+strtrim(hduclas2_rm[0],2)
if keyword_set(hduclas2_eb) then hist=hist+',hduclas2_eb='+strtrim(hduclas2_eb[0],2)
if keyword_set(hduvers) then hist=hist+',hduvers='+strtrim(hduvers[0],2)
if n_elements(tlmin4) gt 0 then hist=hist+',tlmin4='+strtrim(tlmin4[0],2)
if keyword_set(numgrp) then hist=hist+',numgrp='+strtrim(numgrp[0],2)
if keyword_set(numelt) then hist=hist+',numelt='+strtrim(numelt[0],2)
if keyword_set(lo_thres) then hist=hist+',lo_thres='+$
	strtrim(string(lo_thres[0],'(f10.8)'),2)
if keyword_set(hduclas3) then hist=hist+',hduclas3='+strtrim(hduclas3[0],2)
if keyword_set(rmfversn) then hist=hist+',rmfversn='+strtrim(rmfversn[0],2)
if keyword_set(hduvers1) then hist=hist+',hduvers1='+strtrim(hduvers1[0],2)
if keyword_set(hduvers2) then hist=hist+',hduvers2='+strtrim(hduvers2[0],2)
if keyword_set(ccls) then hist=hist+',/ccls'
if keyword_set(ccnm) then hist=hist+',/ccnm'
if keyword_set(cdtp) then hist=hist+',/cdtp'
if keyword_set(cvsd) then hist=hist+',cvsd='+strtrim(cvsd[0],2)
if keyword_set(cvst) then hist=hist+',cvst='+strtrim(cvst[0],2)
if keyword_set(cdes) then hist=hist+',cdes='+strtrim(cdes[0],2)
if keyword_set(origin) then hist=hist+',origin='+strtrim(origin[0],2)
if keyword_set(creator) then hist=hist+',creator='+strtrim(creator[0],2)
if keyword_set(revision) then hist=hist+',revision='+strtrim(revision[0],2)
if keyword_set(mission) then hist=hist+',mission='+strtrim(mission[0],2)
;	and make sure it will fit into the 80 char limit..
ncol=72
;ncol=60
nhist=long(strlen(hist)/ncol) & histo=strarr(nhist+1L)
for i=0L,nhist do histo[i]=strmid(hist,ncol*i,ncol)

;	initialize some variables
n1=n_elements(rmf.ELO) & n2=n_elements(rmf.EMN)
if not keyword_set(origin) then begin
  spawn,'whoami ; hostname',jnk & origin=jnk[0]+'@'+jnk[1]
endif
if not keyword_set(creator) then creator='wrt_ogip_rmf.pro'
szm=size(rmf.matrix) & nummat=szm[1] & numrow=szm[2]

;	figure out the header
if keyword_set(hdrfil) then begin
  jnk=mrdfits(hdrfil[0],0,hdr0)
  jnk=mrdfits(hdrfil[0],1,hdr1)
  jnk=mrdfits(hdrfil[0],2,hdr2)
  if keyword_set(hdr1) then x1=strtrim(sxpar(hdr1,'EXTNAME'),2) else x1=''
  if keyword_set(hdr2) then x2=strtrim(sxpar(hdr2,'EXTNAME'),2) else x2=''
  if x1 eq 'EBOUNDS' then begin
    hdrsm=hdr2 & hdreb=hdr1
  endif else begin
    hdrsm=hdr1 & hdreb=hdr2
  endelse
endif
if not keyword_set(hdr0) then begin
  fxhmake,hdr0,/extend,/date,/INITIALIZE
  fxaddpar,hdr0,'LONGSTRN','OGIP 1.0',$
	'The HEASARC Long String Convention may be used'
  fxaddpar,hdr0,'ORIGIN',origin
endif
sxaddhist,histo,hdr0
if not keyword_set(hdrsm) then fxbhmake,hdrsm,n1,'SPECRESP MATRIX',/date,/INITIALIZE
;
if not keyword_set(hdreb) then fxbhmake,hdreb,n2,'EBOUNDS',/date,/INITIALIZE

;	force appropriate changes to header
tmp=sxpar(hdrsm,'EXTNAME') & c='(MATRIX or SPECRESP MATRIX)'
if not keyword_set(tmp) then begin
  if keyword_set(matext) then fxaddpar,hdrsm,'EXTNAME',matext,c else fxaddpar,hdrsm,'EXTNAME','SPECRESP MATRIX',c
endif
tmp=sxpar(hdrsm,'NAXIS1') & c='Number of rows'
fxaddpar,hdrsm,'NAXIS2',rmf.NNRG,c
tmp=sxpar(hdreb,'EXTNAME') & c='EBOUNDS'
if not keyword_set(tmp) then begin
  if keyword_set(ebext) then fxaddpar,hdreb,'EXTNAME',ebext,c else fxaddpar,hdreb,'EXTNAME','EBOUNDS',c
endif
tmp=sxpar(hdrsm,'TELESCOP') & c='Telescope used'
if not keyword_set(tmp) then begin
  if keyword_set(telescop) then begin
    fxaddpar,hdr0,'TELESCOP',telescop,c
    fxaddpar,hdrsm,'TELESCOP',telescop,c & fxaddpar,hdreb,'TELESCOP',telescop,c
  endif else begin
    fxaddpar,hdr0,'TELESCOP','CHANDRA',c
    fxaddpar,hdrsm,'TELESCOP','CHANDRA',c & fxaddpar,hdreb,'TELESCOP','CHANDRA',c
  endelse
endif
tmp=sxpar(hdrsm,'INSTRUME') & c='Instrument used'
if not keyword_set(tmp) then begin
  if keyword_set(instrume) then begin
    fxaddpar,hdr0,'INSTRUME',instrume,c
    fxaddpar,hdrsm,'INSTRUME',instrume,c & fxaddpar,hdreb,'INSTRUME',instrume,c
  endif else begin
    fxaddpar,hdr0,'INSTRUME','HRC',c
    fxaddpar,hdrsm,'INSTRUME','HRC',c & fxaddpar,hdreb,'INSTRUME','HRC',c
  endelse
endif
tmp=sxpar(hdrsm,'FILTER') & c='the instrument filter in use (if any)'
if not keyword_set(tmp) then begin
  if keyword_set(filter) then begin
    fxaddpar,hdr0,'FILTER',filter,c
    fxaddpar,hdrsm,'FILTER',filter,c & fxaddpar,hdreb,'FILTER',filter,c
  endif else begin
    fxaddpar,hdr0,'FILTER','NONE',c
    fxaddpar,hdrsm,'FILTER','NONE',c & fxaddpar,hdreb,'FILTER','NONE',c
  endelse
endif
tmp=sxpar(hdrsm,'DETNAM') & c='detector used'
if not keyword_set(tmp) then begin
  if keyword_set(DETNAM) then begin
    fxaddpar,hdr0,'DETNAM',DETNAM,c
    fxaddpar,hdrsm,'DETNAM',DETNAM,c & fxaddpar,hdreb,'DETNAM',DETNAM,c
  endif else begin
    fxaddpar,hdr0,'DETNAM','HRC-S',c
    fxaddpar,hdrsm,'DETNAM','HRC-S',c & fxaddpar,hdreb,'DETNAM','HRC-S',c
  endelse
endif
tmp=sxpar(hdrsm,'GRATING') & c='grating used, if any'
if not keyword_set(tmp) then begin
  if keyword_set(GRATING) then begin
    fxaddpar,hdr0,'GRATING',GRATING,c
    fxaddpar,hdrsm,'GRATING',GRATING,c & fxaddpar,hdreb,'GRATING',GRATING,c
  endif else begin
    fxaddpar,hdr0,'GRATING','LETG',c
    fxaddpar,hdrsm,'GRATING','LETG',c & fxaddpar,hdreb,'GRATING','LETG',c
  endelse
endif
tmp=sxpar(hdrsm,'CHANTYPE') & c='(PHA or PI) uncorrected or corrected for gain'
if not keyword_set(tmp) then begin
  if keyword_set(chantype) then begin
    fxaddpar,hdrsm,'CHANTYPE',chantype,c & fxaddpar,hdreb,'CHANTYPE',chantype,c
  endif else begin
    fxaddpar,hdrsm,'CHANTYPE','PI',c & fxaddpar,hdreb,'CHANTYPE','PI',c
  endelse
endif
tmp=sxpar(hdrsm,'DETCHANS') & c='raw detector channels'
if not keyword_set(tmp) then begin
  if keyword_set(detchans) then begin
    fxaddpar,hdrsm,'DETCHANS',detchans,c
  endif else begin
    fxaddpar,hdrsm,'DETCHANS',n_elements(rmf.EMN),c
  endelse
endif
tmp=sxpar(hdrsm,'HDUCLASS') & c='file format is OGIP standard'
if not keyword_set(tmp) then begin
  if keyword_set(hduclass) then begin
    fxaddpar,hdrsm,'HDUCLASS',hduclass,c & fxaddpar,hdreb,'HDUCLASS',hduclass,c
  endif else begin
    fxaddpar,hdrsm,'HDUCLASS','OGIP',c & fxaddpar,hdreb,'HDUCLASS','OGIP',c
  endelse
endif
tmp=sxpar(hdrsm,'HDUCLAS1') & c='extension contains response data'
if not keyword_set(tmp) then begin
  if keyword_set(HDUCLAS1) then begin
    fxaddpar,hdrsm,'HDUCLAS1',HDUCLAS1,c & fxaddpar,hdreb,'HDUCLAS1',HDUCLAS1,c
  endif else begin
    fxaddpar,hdrsm,'HDUCLAS1','RESPONSE',c & fxaddpar,hdreb,'HDUCLAS1','RESPONSE',c
  endelse
endif
tmp=sxpar(hdrsm,'HDUCLAS2') & c='extension contains a response matrix'
if not keyword_set(tmp) then begin
  if keyword_set(HDUCLAS2_RM) then begin
    fxaddpar,hdrsm,'HDUCLAS2',HDUCLAS2_RM,c
  endif else begin
    fxaddpar,hdrsm,'HDUCLAS2','RSP_MATRIX',c
  endelse
endif
tmp=sxpar(hdreb,'HDUCLAS2') & c='extension contains a response matrix'
if not keyword_set(tmp) then begin
  if keyword_set(HDUCLAS2_EB) then begin
    fxaddpar,hdreb,'HDUCLAS2',HDUCLAS2_EB,c
  endif else begin
    fxaddpar,hdreb,'HDUCLAS2','EBOUNDS',c
  endelse
endif
tmp=sxpar(hdrsm,'HDUVERS') & c='version of the file format'
if not keyword_set(tmp) then begin
  if keyword_set(HDUVERS) then begin
    fxaddpar,hdrsm,'HDUVERS',HDUVERS,c
  endif else begin
    fxaddpar,hdrsm,'HDUVERS','1.3.0',c
  endelse
endif
tmp=sxpar(hdreb,'HDUVERS') & c='version of the file format'
if not keyword_set(tmp) then begin
  if keyword_set(HDUVERS) then begin
    fxaddpar,hdreb,'HDUVERS',HDUVERS,c
  endif else begin
    fxaddpar,hdreb,'HDUVERS','1.2.0',c
  endelse
endif
tmp=sxpar(hdrsm,'TLMIN4') & c='first channel in the response'
if not keyword_set(tmp) then begin
  if n_elements(TLMIN4) gt 0 then begin
    fxaddpar,hdrsm,'TLMIN4',TLMIN4[0],c & fxaddpar,hdreb,'TLMIN4',TLMIN4[0],c
  endif else begin
    fxaddpar,hdrsm,'TLMIN4',rmf.FIRSTCHAN,c & fxaddpar,hdreb,'TLMIN4',rmf.FIRSTCHAN,c
  endelse
endif
tmp=sxpar(hdrsm,'NUMGRP') & c='number of channel subsets'
if not keyword_set(tmp) then begin
  if keyword_set(NUMGRP) then begin
    fxaddpar,hdrsm,'NUMGRP',NUMGRP,c
  endif else begin
    fxaddpar,hdrsm,'NUMGRP',long(total(rmf.N_GRP)),c
  endelse
endif
tmp=sxpar(hdrsm,'NUMELT') & c='number of response elements'
if not keyword_set(tmp) then begin
  if keyword_set(NUMELT) then begin
    fxaddpar,hdrsm,'NUMELT',NUMELT,c
  endif else begin
    fxaddpar,hdrsm,'NUMELT',long(total(rmf.N_CHAN)),c
  endelse
endif
tmp=sxpar(hdrsm,'LO_THRES') & c='lower threshold used to construct the matrix'
if not keyword_set(tmp) then begin
  if keyword_set(LO_THRES) then begin
    fxaddpar,hdrsm,'LO_THRES',LO_THRES,c,format='f10.8'
  endif else begin
    fxaddpar,hdrsm,'LO_THRES',0.000001,c,format='f10.8'
  endelse
endif
tmp=sxpar(hdrsm,'HDUCLAS3') & c='REDIST (RMF) or DETECTOR (w.QE,no ARF) or FULL (everything)'
if not keyword_set(tmp) then begin
  if keyword_set(HDUCLAS3) then begin
    fxaddpar,hdrsm,'HDUCLAS3',HDUCLAS3,c
  endif else begin
    fxaddpar,hdrsm,'HDUCLAS3','DETECTOR',c
  endelse
endif
tmp=sxpar(hdrsm,'RMFVERSN') & c='obsolete'
if not keyword_set(tmp) then begin
  if keyword_set(RMFVERSN) then begin
    fxaddpar,hdrsm,'RMFVERSN',RMFVERSN,c & fxaddpar,hdreb,'RMFVERSN',RMFVERSN,c
  endif else begin
    fxaddpar,hdrsm,'RMFVERSN','1992a',c & fxaddpar,hdreb,'RMFVERSN','1992a',c
  endelse
endif
tmp=sxpar(hdrsm,'HDUVERS1') & c='obsolete'
if not keyword_set(tmp) then begin
  if keyword_set(HDUVERS1) then begin
    fxaddpar,hdrsm,'HDUVERS1',HDUVERS1,c & fxaddpar,hdreb,'HDUVERS1',HDUVERS1,c
  endif else begin
    fxaddpar,hdrsm,'HDUVERS1','1.1.0',c & fxaddpar,hdreb,'HDUVERS1','1.1.0',c
  endelse
endif
tmp=sxpar(hdrsm,'HDUVERS2') & c='obsolete'
if not keyword_set(tmp) then begin
  if keyword_set(HDUVERS2) then begin
    fxaddpar,hdrsm,'HDUVERS2',HDUVERS2,c & fxaddpar,hdreb,'HDUVERS2',HDUVERS2,c
  endif else begin
    fxaddpar,hdrsm,'HDUVERS2','1.2.0',c & fxaddpar,hdreb,'HDUVERS2','1.2.0',c
  endelse
endif
tmp=sxpar(hdrsm,'CCLS0001') & c='OGIP-class of this calibration file'
if not keyword_set(tmp) then begin
  if keyword_set(CCLS) then begin
    fxaddpar,hdrsm,'CCLS0001','CPF',c & fxaddpar,hdreb,'CCLS0001','CPF',c
  endif
endif
tmp=sxpar(hdrsm,'CCNM0001') & c='CIF codename for this type of calibration dataset'
if not keyword_set(tmp) then begin
  if keyword_set(CCNM) then begin
    fxaddpar,hdrsm,'CCNM0001','MATRIX',c & fxaddpar,hdreb,'CCNM0001','EBOUNDS',c
  endif
endif
tmp=sxpar(hdrsm,'CDTP0001') & c='OGIP code for the form of the contents of the file'
if not keyword_set(tmp) then begin
  if keyword_set(CDTP) then begin
    fxaddpar,hdrsm,'CDTP0001','DATA',c & fxaddpar,hdreb,'CDTP0001','DATA',c
  endif
endif
tmp=sxpar(hdrsm,'CVSD0001') & c='UTC date of beginning applicability'
caldat,systime(/julian),mon,day,yr,hr,min,sec
if not keyword_set(tmp) then begin
  if keyword_set(CVSD) then begin
    fxaddpar,hdrsm,'CVSD0001',CVSD,c & fxaddpar,hdreb,'CVSD0001',CVSD,c
  endif else begin
    yyyymmdd=string(yr,'(i4)')+'-'+string(mon,'(i2)')+'-'+string(day,'(i2)')
    fxaddpar,hdrsm,'CVSD0001',yyyymmdd,c & fxaddpar,hdreb,'CVSD0001',yyyymmdd,c
  endelse
endif
tmp=sxpar(hdrsm,'CVST0001') & c='UTC time of beginning applicability'
if not keyword_set(tmp) then begin
  if keyword_set(CVST) then begin
    fxaddpar,hdrsm,'CVST0001',CVST,c & fxaddpar,hdreb,'CVST0001',CVST,c
  endif else begin
    hhmmss=string(hr,'(i2)')+':'+string(min,'(i2)')+':'+string(sec,'(i2)')
    fxaddpar,hdrsm,'CVST0001',hhmmss,c & fxaddpar,hdreb,'CVST0001',hhmmss,c
  endelse
endif
tmp=sxpar(hdrsm,'CDES0001') & c='brief descriptive summary of this dataset'
if not keyword_set(tmp) then begin
  if keyword_set(CDES) then begin
    fxaddpar,hdrsm,'CDES0001',CDES,c & fxaddpar,hdreb,'CDES0001',CDES,c
  endif else begin
    fxaddpar,hdrsm,'CDES0001','RSP',c & fxaddpar,hdreb,'CDES0001','RSP',c
  endelse
endif
tmp=sxpar(hdrsm,'ORIGIN') & c='source of FITS file'
if not keyword_set(tmp) then begin
  fxaddpar,hdrsm,'ORIGIN',ORIGIN,c & fxaddpar,hdreb,'ORIGIN',ORIGIN,c
endif
tmp=sxpar(hdrsm,'CREATOR') & c='program creating this file'
if not keyword_set(tmp) then begin
  fxaddpar,hdrsm,'CREATOR',CREATOR,c & fxaddpar,hdreb,'CREATOR',CREATOR,c
endif
tmp=sxpar(hdrsm,'REVISION') & c='Revision'
if not keyword_set(tmp) then begin
  if keyword_set(REVISION) then begin
    fxaddpar,hdrsm,'REVISION',REVISION,c & fxaddpar,hdreb,'REVISION',REVISION,c
  endif else begin
    fxaddpar,hdrsm,'REVISION',1,c & fxaddpar,hdreb,'REVISION',1,c
  endelse
endif
tmp=sxpar(hdrsm,'MISSION') & c='Mission'
if not keyword_set(tmp) then begin
  if keyword_set(MISSION) then begin
    fxaddpar,hdrsm,'MISSION',MISSION,c & fxaddpar,hdreb,'MISSION',MISSION,c
  endif else begin
    fxaddpar,hdrsm,'MISSION','AXAF',c & fxaddpar,hdreb,'MISSION','AXAF',c
  endelse
endif

;	define the columns in MATRIX extension
elo=float(rmf.ELO) & ehi=float(rmf.EHI)
n_grp=fix(rmf.N_GRP) & f_chan=fix(rmf.F_CHAN) & n_chan=fix(rmf.N_CHAN)
matrix=float(rmf.MATRIX)
szf=size(f_chan)
tmp=sxpar(hdrsm,'TTYPE1')
if keyword_set(tmp) then begin
  if tmp ne 'ENERG_LO' then fxaddpar,hdrsm,'TTYPE1','ENERG_LO','Energy bin lower bound'
  fxaddpar,hdrsm,'TUNIT1','keV'
endif else $
  fxbaddcol,i1,hdrsm,elo[0],'ENERG_LO','Energy bin lower bound',TUNIT='keV'
fxaddpar,hdrsm,'TFORM1','1E','format of field'
fxaddpar,hdrsm,'TLMIN1',min(rmf.ELO)
fxaddpar,hdrsm,'TLMAX1',max(rmf.ELO)
tmp=sxpar(hdrsm,'TTYPE2')
if keyword_set(tmp) then begin
  if tmp ne 'ENERG_HI' then fxaddpar,hdrsm,'TTYPE2','ENERG_HI','Energy bin upper bound'
  fxaddpar,hdrsm,'TUNIT2','keV'
endif else $
  fxbaddcol,i2,hdrsm,ehi[0],'ENERG_HI','Energy bin upper bound',TUNIT='keV'
fxaddpar,hdrsm,'TFORM2','1E','format of field'
fxaddpar,hdrsm,'TLMIN2',min(rmf.EHI)
fxaddpar,hdrsm,'TLMAX2',max(rmf.EHI)
tmp=sxpar(hdrsm,'TTYPE3')
if keyword_set(tmp) then begin
  if tmp ne 'N_GRP' then fxaddpar,hdrsm,'TTYPE3','N_GRP','number of channel subsets'
endif else $
  fxbaddcol,i3,hdrsm,n_grp[0],'N_GRP','number of channel subsets'
fxaddpar,hdrsm,'TFORM3','I','format of field'
fxaddpar,hdrsm,'TLMIN3',min([0,(rmf.N_GRP)[*]])
fxaddpar,hdrsm,'TLMAX3',max(rmf.N_GRP)
tmp=sxpar(hdrsm,'TTYPE4')
if keyword_set(tmp) then begin
  if tmp ne 'F_CHAN' then fxaddpar,hdrsm,'TTYPE4','F_CHAN','channel number of first channel'
endif else $
  if szf[0] eq 1 then fxbaddcol,i4,hdrsm,F_CHAN[0],$
	'F_CHAN','channel number of first channel',/variable else $
	fxbaddcol,i4,hdrsm,reform(F_CHAN[*,0]),$
	'F_CHAN','channel number of first channel',/variable
fxaddpar,hdrsm,'TFORM4','PI('+strtrim(max(rmf.N_GRP),2)+')','format of field; variable length array'
;fxaddpar,hdrsm,'TLMIN4',min(rmf.F_CHAN)
fxaddpar,hdrsm,'TLMAX4',max(rmf.F_CHAN)
tmp=sxpar(hdrsm,'TTYPE5')
if keyword_set(tmp) then begin
  if tmp ne 'N_CHAN' then fxaddpar,hdrsm,'TTYPE5','N_CHAN','number of channels in each group'
endif else $
  if szf[0] eq 1 then fxbaddcol,i5,hdrsm,N_CHAN[0],$
	'N_CHAN','number of channels in each group',/variable else $
	fxbaddcol,i5,hdrsm,reform(N_CHAN[*,0]),$
	'N_CHAN','number of channels in each group',/variable
fxaddpar,hdrsm,'TFORM5','PI('+strtrim(max(rmf.N_GRP),2)+')','format of field; variable length array'
fxaddpar,hdrsm,'TLMIN5',min([0,(rmf.N_CHAN)[*]])
fxaddpar,hdrsm,'TLMAX5',max(rmf.N_CHAN)
tmp=sxpar(hdrsm,'TTYPE6')
if keyword_set(tmp) then begin
  if tmp ne 'MATRIX' then fxaddpar,hdrsm,'TTYPE6','MATRIX','response values'
endif else $
  fxbaddcol,i6,hdrsm,reform(MATRIX[*,0]),'MATRIX','response values',/variable
fxaddpar,hdrsm,'TFORM6','PE('+strtrim(nummat,2)+')','format of field; variable length array'
fxaddpar,hdrsm,'TLMIN6',min(rmf.MATRIX)
fxaddpar,hdrsm,'TLMAX6',max(rmf.MATRIX)

;	define the columns in EBOUNDS extension
emn=float(rmf.EMN) & emx=float(rmf.EMX)
tmp=sxpar(hdreb,'TTYPE1')
if keyword_set(tmp) then begin
  if tmp ne 'CHANNEL' then fxaddpar,hdreb,'TTYPE1','CHANNEL','raw channel number'
endif else $
  fxbaddcol,j1,hdreb,1,'CHANNEL','raw channel number'
fxaddpar,hdreb,'TFORM1','I','format of field'
fxaddpar,hdreb,'TLMIN1',1
fxaddpar,hdreb,'TLMAX1',n2
tmp=sxpar(hdreb,'TTYPE2')
if keyword_set(tmp) then begin
  if tmp ne 'E_MIN' then fxaddpar,hdreb,'TTYPE2','E_MIN','detector channel bin lower bound'
  fxaddpar,hdreb,'TUNIT1','keV'
endif else $
  fxbaddcol,j2,hdreb,EMN[0],'E_MIN','detector channel bin lower bound',TUNIT='keV'
fxaddpar,hdreb,'TFORM2','1E','format of field'
fxaddpar,hdreb,'TLMIN2',min(rmf.EMN)
fxaddpar,hdreb,'TLMAX2',max(rmf.EMN)
tmp=sxpar(hdreb,'TTYPE3')
if keyword_set(tmp) then begin
  if tmp ne 'E_MAX' then fxaddpar,hdreb,'TTYPE3','E_MAX','detector channel bin upper bound'
  fxaddpar,hdreb,'TUNIT2','keV'
endif else $
  fxbaddcol,j3,hdreb,EMX[0],'E_MAX','detector channel bin upper bound',TUNIT='keV'
fxaddpar,hdreb,'TFORM3','1E','format of field'
fxaddpar,hdreb,'TLMIN3',min(rmf.EMX)
fxaddpar,hdreb,'TLMAX3',max(rmf.EMX)

;	write the primary header
message,'Writing to file:'+rmfile[0],/info
fxwrite,rmfile[0],hdr0

;	write the MATRIX extension
fxbcreate,um,rmfile[0],hdrsm
for i=0L,numrow-1L do begin
  fxbwrite,um,elo[i],1,i+1L
  fxbwrite,um,ehi[i],2,i+1L
  fxbwrite,um,n_grp[i],3,i+1L
  if szf[0] eq 1 then fxbwrite,um,f_chan[i],4,i+1L else $
    fxbwrite,um,reform(f_chan[*,i]),4,i+1L
  if szf[0] eq 1 then fxbwrite,um,n_chan[i],5,i+1L else $
    fxbwrite,um,reform(n_chan[*,i]),5,i+1L
  fxbwrite,um,reform(matrix[*,i]),6,i+1L
endfor
fxbfinish,um

;	write the EBOUNDS extension
fxbcreate,ue,rmfile[0],hdreb
for i=0L,n_elements(emn)-1L do begin
  fxbwrite,ue,fix(i+1),1,i+1L
  fxbwrite,ue,emn[i],2,i+1L
  fxbwrite,ue,emx[i],3,i+1L
endfor
fxbfinish,ue

return
end
