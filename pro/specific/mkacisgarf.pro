function mkacisgarf,offset,wgrid,order=order,arm=arm,caldb=caldb,$
	hrmaea=hrmaea,acisqe=acisqe,greff=greff,contam=contam,$
	tstart=tstart,verbose=verbose, _extra=e
;+
;function	mkacisgarf
;	returns effective area as a function of wavelength for specified
;	CXC grating and order on the ACIS-S
;
;syntax
;	arf=mkacisgarf(offset,wgrid,order=order,arm=arm,caldb=caldb,$
;	hrmaea=hrmaea,acisqe=acisqe,greff=greff,contam=contam,$
;	tstart=tstart,verbose=verbose)
;
;parameters
;	offset	[INPUT; required] offset of source position from the
;		nominal aim point on S3 [arcmin]
;		* +ve is towards the center of S3
;		* someday, after I figure out how to do it, I'll change
;		  the required inputs to (SRC_X,SRC_Y), and possibly ROLL.
;		  meanwhile, this parameter goes straight to ACISS_GAPS.
;	wgrid	[OUTPUT] wavelengths at which effective area is computed
;
;keywords
;	order	[INPUT] diffraction order.  if not set, assumes +ve 1st order.
;	arm	[INPUT] string specifying which grating arm to be used.
;		* 'MEG', 'HEG', or 'LEG', in decreasing precedence
;		* if illegal or undecipherable, NONE is assumed
;		* Note that MEG intercepts rays from the two outer HRMA shells
;		  (1 and 3, denoted 1100) and the HEG intercepts rays from the
;		  two inner shells (4 and 6, denoted 0011)
;	caldb	[INPUT] root $CALDB directory
;		* default: /soft/ciao/CALDB/
;		* assumes that the following tree exists:
;		- CALDBhrma=CALDB/data/chandra/default
;		- CALDBacis=CALDB/data/chandra/acis/
;		- CALDBhetg=CALDB/data/chandra/default/greff/
;		- CALDBletg=CALDB/data/chandra/default/greff/
;	hrmaea	[INPUT] full path name of file containing HRMA effective areas
;		* must be a FITS file with extensions AXAF_EFFEA1..5, one
;		  for each mirror shell (viz. header keyword SHELL), and
;		  each extension having the columns
;		  	ENERG_LO,ENERG_HI,EFFAREA
;		* if not specified, assumed to be
;			CALDBhrma/effarea/hrmaD1996-12-20axeffaN0008.fits
;	acisqe	[INPUT] full path name of file containing ACIS QEs
;		* must be a FITS file with an extension for each chip
;		  (AXAF_QE1..10), and each extension having columns
;		  	ENERGY,QE
;		  (latter assumed to be product of detector QE and
;		  optical blocking filter)
;		* if not specified, assumed to be
;			CALDBacis/qe/acisD1997-04-17qeN0006.fits
;	greff	[INPUT] full path name of file containing grating efficiencies
;		* must be a FITS file with an extension (AXAF_GREFF1..4),
;		  for each mirror shell (viz., header keyword SHELL), and
;		  each extension having the columns
;		  	ENERGY, EFF[--..-1,0,+1,..++]
;		* if not specified, assumed to be
;		  ARM='LEG':	CALDBletg/letgD1996-11-01greffpr001N0005.fits
;		  else: 	CALDBhetg/hetgD1996-11-01greffpr001N0006.fits
;	contam	[INPUT] full path name of file containing ACIS contamination
;		coefficients
;		* must be FITS file with 3 extensions, one each for C-K, O-K, and F-K
;		  edges, with columns
;		  	FEATURE,N_ENERGY,ENERGY,ALPHA,N_TIME,TIME,BETA
;		* if not specified, assumed to be
;		  	CALDBacis/contam/acisD1997-04-17contamN0006.fits
;		* NOTE: this file does not exist in caldb yet, so for the
;		  nonce, the CONTAM file is ignored unless explicitly specified
;		* NOTE: the contamination is a multiplicative factor to the ARF,
;		  exp(SUM{-alpha(E)*beta}), where beta is a time dependent factor
;	tstart	[INPUT] the time of observation to be used to determine BETA for
;		contamination corrections, in units of spacecraft time in [s]
;		since launch.
;		* if not specified, the latest time encoded in CONTAM is assumed
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	does not handle the 0th order
;	does not handle ACIS-I CCDs being used as the spectroscopic array
;
;subroutines
;	ACISS_GAPS
;
;history
;	modified from ACISGARF (VK; MMIVjan)
;	updated default files (VK; MMVIIsep)
;-

;	usage
ok='ok' & np=n_params() & noff=n_elements(offset)
if np eq 0 then ok='Insufficient parameters' else $
 if noff eq 0 then ok='OFFSET is undefined' else $
  if noff gt 1 then ok='cannot handle array of OFFSETs'
if ok ne 'ok' then begin
  print,'Usage: effar=mkacisgarf(offset,wvlar,order=order,arm=arm,caldb=caldb,$'
  print,'  hrmaea=hrmaea,acisqe=acisqe,greff=greff,contam=contam,tstart=tstart,$
  print,'  verbose=verbose)'
  print,'  return effective areas for ACIS/HEG, ACIS/MEG, or ACIS/LEG'
  if np ne 0 then message,ok,/informational
  return,-1L
endif
srcpos=offset(0)

;	keywords

;	verbosity
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	order
gord=1 & if n_elements(order) gt 0 then gord=fix(order[0])

;	grating arm
garm='NONE' & sza=size(arm) & nsza=n_elements(sza)
if sza(nsza-2) eq 7 then begin
  if strpos(strlowcase(strtrim(arm[0],2)),'leg',0) ge 0 then garm='LEG'
  if strpos(strlowcase(strtrim(arm[0],2)),'heg',0) ge 0 then garm='HEG'
  if strpos(strlowcase(strtrim(arm[0],2)),'meg',0) ge 0 then garm='MEG'
endif

;	$CALDB
dbcal=getenv('CALDB')
if not keyword_set(dbcal) then dbcal='/soft/ciao/CALDB/'
if keyword_set(caldb) then begin
  if size(caldb,/type) eq 7 then dbcal=strtrim(caldb[0],2)
endif
CALDBhrma=dbcal+'/data/chandra/default/'
CALDBacis=dbcal+'/data/chandra/acis/'
CALDBhetg=dbcal+'/data/chandra/default/greff/'
CALDBletg=dbcal+'/data/chandra/default/greff/'

;	mirror effective area
hrmafil=CALDBhrma+'/axeffa/hrmaD1996-12-20axeffaN0008.fits'
if keyword_set(hrmaea) then begin
  if size(hrmaea,/type) eq 7 then hrmafil=strtrim(hrmaea[0],2)
endif
if vv gt 1 then $
  message,'reading HRMA effective areas from: '+hrmafil,/informational

;	ACIS QE
acisqefil=CALDBacis+'/qe/acisD1997-04-17qeN0006.fits'
if keyword_set(acisqe) then begin
  if size(acisqe,/type) eq 7 then acisqefil=strtrim(acisqe[0],2)
endif
if vv gt 1 then $
  message,'reading ACIS QE from: '+acisqefil,/informational

;	grating efficiencies
if garm eq 'LEG' then $
  greffil=CALDBletg+'/letgD1996-11-01greffpr001N0005.fits' else $
  greffil=CALDBhetg+'/hetgD1996-11-01greffpr001N0006.fits'
if keyword_set(greff) then begin
  if size(greff,/type) eq 7 then greffil=strtrim(greff[0],2)
endif
if vv gt 1 then $
  message,'reading grating efficiencies from: '+greffil,/informational

;	contamination
;contamfil='/data/draft_caldb/caldbmgr/AO-6_effarea/2297/primary/acisD1999-08-13contamN0002.fits'
;contamfil='/data/draft_caldb/CALDB/data/chandra/acis/bcf/contam/acisD1999-08-13contamN0003.fits'
contamfil=CALDBacis+'/contam/acisD1997-04-17contamN0006.fits'
if keyword_set(contam) then begin
  if size(contam,/type) eq 7 then contamfil=strtrim(contam[0],2)
endif
if vv gt 1 then $
  message,'reading contamination coefficients from: '+contamfil,/informational

;{	initialize

legmin=3.0 & legmax=180.0	;LEG 0th order wavelength range
legbin=0.0075			;LEG wavelength grid spacing
megmin=1.0 & megmax=42.0	;MEG 0th order wavelength range
megbin=0.005			;MEG wavelength grid spacing
hegmin=0.5 & hegmax=21.0	;HEG 0th order wavelength range
hegbin=0.0025			;HEG wavelength grid spacing
if keyword_set(gord) then begin
  legmin=legmin/abs(gord) & legmax=legmax/abs(gord) & legbin=legbin/abs(gord)
  megmin=megmin/abs(gord) & megmax=megmax/abs(gord) & megbin=megbin/abs(gord)
  hegmin=hegmin/abs(gord) & hegmax=hegmax/abs(gord) & hegbin=hegbin/abs(gord)
endif

;	reference wavelength grid
nleg=(legmax-legmin)/legbin & wleg=findgen(nleg)*legbin+legmin
nmeg=(megmax-megmin)/megbin & wmeg=findgen(nmeg)*megbin+megmin
nheg=(hegmax-hegmin)/hegbin & wheg=findgen(nheg)*hegbin+hegmin

;	chip gap locations
aciss_gaps,hgap,mgap,lgap,onchip=onchip,offset=srcpos,order=gord
;input OFFSET [arcmin], outputs ONCHIP [S-array#], ?GAP [Ang]

case garm of
  'LEG': begin
    ngrid=nleg & wgrid=wleg & gaps=lgap
  end
  'MEG': begin
    ngrid=nmeg & wgrid=wmeg & gaps=mgap
  end
  'HEG': begin
    ngrid=nheg & wgrid=wheg & gaps=hgap
  end
  else: begin
    ;	this should NEVER toggle
    message,'Watchewtakinaboot, Willis?',/info
    return,-1L
  end
endcase

;	initialzed}

;	file existence checks
ok='ok'
tmp=findfile(hrmafil,count=nhrma)
tmp=findfile(acisqefil,count=nacis)
tmp=findfile(greffil,count=ngref)
if keyword_set(contam) then begin	;TEMPORARY
tmp=findfile(contamfil,count=ncontam)
endif					;TEMPORARY
if nhrma eq 0 then ok=hrmafil+': HRMA effective area file does not exist'
if nacis eq 0 then ok=acisqefil+': ACIS QE file does not exist'
if ngref eq 0 then ok=greffil+': grating efficiency file does not exist'
if keyword_set(contam) then begin	;TEMPORARY
if ncontam eq 0 then ok=contamfil+': contamination coefficients file does not exist'
endif					;TEMPORARY
if nhrma gt 1 then ok=hrmafil+': too many files'
if nacis gt 1 then ok=acisqefil+': too many files'
if ngref gt 1 then ok=greffil+': too many files'
if keyword_set(contam) then begin	;TEMPORARY
if ncontam gt 1 then ok=contamfil+': too many files'
endif					;TEMPORARY
if ok ne 'ok' then begin
  message,ok,/informational & return,0*wgrid
endif

;	read in CALDB files
;	HRMA EA
for i=1,5 do begin
  tmp=mrdfits(hrmafil,i,hdr)
  shell=strtrim(sxpar(hdr,'SHELL'),2)
  case shell of
    '1000': hrma1000=tmp
    '0100': hrma0100=tmp
    '0010': hrma0010=tmp
    '0001': hrma0001=tmp
    '1111': ;ignore this case
    else: begin
      message,hrmafil+': could not read HRMA EA',/informational
      if vv gt 100 then stop,'HALTing.  type .CON to continue'
      return,0*wgrid
    end
  endcase
endfor
;	ACIS QEs
for i=5,10 do begin
  tmp=mrdfits(acisqefil,i,hdr)
  detnam=strtrim(sxpar(hdr,'DETNAM'),2)
  case detnam of
    'ACIS-0': acisqe01=tmp
    'ACIS-1': acisqe02=tmp
    'ACIS-2': acisqe03=tmp
    'ACIS-3': acisqe04=tmp
    'ACIS-4': acisqe05=tmp
    'ACIS-5': acisqe06=tmp
    'ACIS-6': acisqe07=tmp
    'ACIS-7': acisqe08=tmp
    'ACIS-8': acisqe09=tmp
    'ACIS-9': acisqe10=tmp
    else: begin
      message,acisqefil+': could not read ACIS QE',/informational
      if vv gt 100 then stop,'HALTing.  type .CON to continue'
      return,0*wgrid
    endelse
  endcase
endfor
;	grating efficiencies
for i=1,4 do begin
  tmp=mrdfits(greffil,i,hdr)
  shell=strtrim(sxpar(hdr,'SHELL'),2)
  case shell of
    '1000': gref1000=tmp
    '0100': gref0100=tmp
    '0010': gref0010=tmp
    '0001': gref0001=tmp
    else: begin
      message,greffil+': could not read grating efficiecies',/informational
      if vv gt 100 then stop,'HALTing.  type .CON to continue'
      return,0*wgrid
    end
  endcase
endfor
;	contamination coefficients
if keyword_set(contam) then begin	;TEMPORARY
cntmdat=mrdfits(contamfil,1,hdr)
oCK=where(cntmdat.(0) eq 'CK',moCK)
oOK=where(cntmdat.(0) eq 'OK',moOK)
oFK=where(cntmdat.(0) eq 'FK',moFK)
oFF=where(cntmdat.(0) eq 'FF',moFF)
if mocK eq 0 or moOK eq 0 or moFK eq 0 then begin
  message,contamfil+': not in the right format',/informational
  return,0*wgrid
endif
nrow=n_elements(cntmdat.n_energy) & start_time=0.D
start_time=max(cntmdat.TIME) & if keyword_set(tstart) then start_time=double(tstart[0])
if start_time gt max(cntmdat.TIME) then begin
  message,'CONTAM file does not reach specified TSTART; setting to maximum',/informational
  start_time=max(cntmdat.TIME)
endif
if start_time lt min(cntmdat.TIME) then begin
  message,'CONTAM file does not reach specified TSTART; setting to minimum',/informational
  start_time=min(cntmdat.TIME)
endif
for i=0,nrow-1 do begin
  clambda=12.3985/((cntmdat.ENERGY)[0L:(cntmdat.N_ENERGY)[i]-1L,i])
  if i eq 0 then lambda=clambda else lambda=[lambda,clambda]
endfor
os=sort(lambda) & cexp=0.*lambda & lambda=lambda[os]
for i=0,nrow-1 do begin
  clambda=12.3985/((cntmdat.ENERGY)[0L:(cntmdat.N_ENERGY)[i]-1L,i])
  calpha=(cntmdat.ALPHA)[0L:(cntmdat.N_ENERGY)[i]-1L,i]
  ctime=(cntmdat.TIME)[0L:(cntmdat.N_TIME)[i]-1L,i]
  cbeta=(cntmdat.BETA)[0L:(cntmdat.N_TIME)[i]-1L,i]
  xbeta=interpol(cbeta,ctime,start_time) > 0
  cexp=cexp+interpol(xbeta[0]*calpha,clambda,lambda)
  if vv gt 400 then stop,i
endfor
cexp=cexp>0
cfactor=exp(-cexp)
endif					;TEMPORARY

;	interpolate everything onto the same wavelength grid
for i=1,4 do begin
  ;HRMA: ENERG_LO,ENERG_HI,EFFAREA
  ;GREFF: ENERGY, EFF[--..-1,0,+1,..++]
  case i of
    1: begin & hrma=hrma1000 & gref=gref1000 & end
    2: begin & hrma=hrma0100 & gref=gref0100 & end
    3: begin & hrma=hrma0010 & gref=gref0010 & end
    4: begin & hrma=hrma0001 & gref=gref0001 & end
  endcase
  egrid=0.5*(hrma.energ_lo+hrma.energ_hi)
  h_ea=interpol(hrma.effarea,12.3985/egrid,wgrid)
  iord=(size(gref.EFF))[1]/2 & oord=lindgen(2L*iord+1L)-iord
  if abs(gord) gt iord then begin
    message,'Grating efficiencies only present up to order '+$
    	strtrim(iord,2),/informational
    return,0*wgrid
  endif
  kord=gord+iord
  gr=interpol(gref.eff[kord,*],12.3985/gref.energy,wgrid)
  case i of
    1: begin & hea1000=h_ea & gr1000=gr & end
    2: begin & hea0100=h_ea & gr0100=gr & end
    3: begin & hea0010=h_ea & gr0010=gr & end
    4: begin & hea0001=h_ea & gr0001=gr & end
  endcase
endfor
for i=5,10 do begin
  ;QE: ENERGY,QE
  case i of
    1: aqe01=interpol(acisqe01.QE,12.3985/acisqe01.energy,wgrid)
    2: aqe02=interpol(acisqe02.QE,12.3985/acisqe02.energy,wgrid)
    3: aqe03=interpol(acisqe03.QE,12.3985/acisqe03.energy,wgrid)
    4: aqe04=interpol(acisqe04.QE,12.3985/acisqe04.energy,wgrid)
    5: aqe05=interpol(acisqe05.QE,12.3985/acisqe05.energy,wgrid)
    6: aqe06=interpol(acisqe06.QE,12.3985/acisqe06.energy,wgrid)
    7: aqe07=interpol(acisqe07.QE,12.3985/acisqe07.energy,wgrid)
    8: aqe08=interpol(acisqe08.QE,12.3985/acisqe08.energy,wgrid)
    9: aqe09=interpol(acisqe09.QE,12.3985/acisqe09.energy,wgrid)
    10: aqe10=interpol(acisqe10.QE,12.3985/acisqe10.energy,wgrid)
  endcase
endfor
contamfac=0.*wgrid+1.
if keyword_set(contam) then begin	;TEMPORARY
oo=where(wgrid gt min(lambda) and wgrid lt max(lambda),moo)
if moo gt 0 then contamfac[oo]=interpol(cfactor,lambda,wgrid[oo]) > 0
oo=where(wgrid ge max(lambda),moo)
if moo gt 0 then begin
  ok=where(lambda gt 44,mok)
  if mok gt 0 then begin
    ab=linfit(lambda[ok],alog(cfactor[ok]))
    contamfac[oo]=exp(ab[0]+ab[1]*wgrid[oo])
  endif else contamfac[oo]=0.
endif
endif					;TEMPORARY

if vv gt 200 then stop,'HALTing, type .CON to continue'

;	compute effective area for each chip over specified wavelength grid
;	EA_CHIP = CHIP_QE * Sum_SHELL (HRMA_EA*GRATING_EFFICIENCIES)

case garm of
  'MEG': sumshell=hea1000*gr1000+hea0100*gr0100
  'HEG': sumshell=hea0010*gr0010+hea0001*gr0001
  else: sumshell=hea1000*gr1000+hea0100*gr0100+hea0010*gr0010+hea0001*gr0001
endcase


;	figure out which chips we have to deal with and multiply the
;	appropriate chip QE with HRMA effar and grating efficiency
garf=0.*wgrid
if gord lt 0 and onchip lt 0 then return,garf	;spectrum off of array
if gord gt 0 and onchip gt 5 then return,garf	;spectrum off of array
if gord eq 0 then begin
  case onchip of
    0: garf=aqe05*sumshell*contamfac
    1: garf=aqe06*sumshell*contamfac
    2: garf=aqe07*sumshell*contamfac
    3: garf=aqe08*sumshell*contamfac
    4: garf=aqe09*sumshell*contamfac
    5: garf=aqe10*sumshell*contamfac
    else: begin
      help,onchip
      message,'0th order not on any known ACIS-S chip',/informational
      garf=0*wgrid
    end
  endcase
  return,garf
endif
for i=5,10 do begin
  j=i-5	;the S-array chip#
  w0=max(gaps)+10. & w1=w0
  if gord lt 0 and j le onchip then begin
    w1=gaps[2*j] & w0=0. & if onchip gt j then w0=gaps[2*j+1]
  endif
  if gord gt 0 and j ge onchip then begin
    w1=gaps[2*j+1] & w0=0. & if onchip lt j then w0=gaps[2*j]
  endif
  ow=where(wgrid ge w0 and wgrid le w1,mow)
  if mow gt 0 then begin
    case i of
      ;1: garf[ow]=aqe01[ow]*sumshell[ow]*contamfac[ow]
      ;2: garf[ow]=aqe02[ow]*sumshell[ow]*contamfac[ow]
      ;3: garf[ow]=aqe03[ow]*sumshell[ow]*contamfac[ow]
      ;4: garf[ow]=aqe04[ow]*sumshell[ow]*contamfac[ow]
      5: garf[ow]=aqe05[ow]*sumshell[ow]*contamfac[ow]
      6: garf[ow]=aqe06[ow]*sumshell[ow]*contamfac[ow]
      7: garf[ow]=aqe07[ow]*sumshell[ow]*contamfac[ow]
      8: garf[ow]=aqe08[ow]*sumshell[ow]*contamfac[ow]
      9: garf[ow]=aqe09[ow]*sumshell[ow]*contamfac[ow]
      10: garf[ow]=aqe10[ow]*sumshell[ow]*contamfac[ow]
      else: begin
	message,'ACIS-S'+strtrim(j,2)+$
	': cannot understand which CCD we are on',/informational
	return,garf
      end
    endcase
  endif
  if vv gt 300 then stop,i,j,w0,w1,'	HALTing; type .CON to continue'
endfor

if vv gt 200 then stop,'HALTing, type .CON to continue'

return,garf
end
