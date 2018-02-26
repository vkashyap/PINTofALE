function mkhrciarf,wgrid,order=order,arm=arm,caldb=caldb,$
	hrmaea=hrmaea,qefil=qefil,greff=greff,verbose=verbose,$
	_extra=e
;+
;function	mkhrciarf
;	returns effective area as a function of wavelength for specified
;	CXC grating and order on the HRC-I
;
;syntax
;	arf=mkhrciarf(wgrid,order=order,arm=arm,caldb=caldb,$
;	hrmaea=hrmaea,qefil=qefil,greff=greff,verbose=verbose)
;
;parameters
;	wgrid	[OUTPUT] wavelengths at which effective area is computed
;
;keywords
;	order	[INPUT] diffraction order.  if not set, assumes +ve 1st order.
;	arm	[INPUT] string specifying which grating arm to be used.
;		* 'NONE', 'MEG', 'HEG', or 'LEG', in decreasing precedence
;		* if illegal or undecipherable, NONE is assumed
;		* Note that MEG intercepts rays from the two outer HRMA shells
;		  (1 and 3, denoted 1100) and the HEG intercepts rays from the
;		  two inner shells (4 and 6, denoted 0011)
;	caldb	[INPUT] root $CALDB directory
;		* default: /soft/ciao/CALDB/
;		* assumes that the following tree exists:
;		- CALDBhrma=CALDB/data/chandra/default/
;		- CALDBhrc=CALDB/data/chandra/hrc/
;		- CALDBhetg=CALDB/data/chandra/default/greff/
;		- CALDBletg=CALDB/data/chandra/default/greff/
;	hrmaea	[INPUT] full path name of file containing HRMA effective areas
;		* must be a FITS file with extensions AXAF_EFFEA1..5, one
;		  for each mirror shell (viz. header keyword SHELL), and
;		  each extension having the columns
;		  	ENERG_LO,ENERG_HI,EFFAREA
;		* if not specified, assumed to be
;			CALDBhrma/axeffa/hrmaD1996-12-20axeffaN0008.fits
;	qefil	[INPUT] full path name of file containing HRC-I QEs
;		* must be a FITS file with an extension for each chip
;		  (AXAF_QE1..4; but currently only AXAF_QE1), and each
;		  extension having columns ENERGY,QE
;		  (latter assumed to be product of detector QE and
;		  optical blocking filter)
;		* if not specified, assumed to be
;			CALDBhrc/qe/hrciD1999-07-22qeN0007.fits
;	greff	[INPUT] full path name of file containing grating efficiencies
;		* must be a FITS file with an extension (AXAF_GREFF1..4),
;		  for each mirror shell (viz., header keyword SHELL), and
;		  each extension having the columns
;		  	ENERGY, EFF[--..-1,0,+1,..++]
;		* if not specified, assumed to be
;		  ARM='NONE':	ignored
;		  ARM='LEG':	CALDBletg/letgD1996-11-01greffpr001N0005.fits
;		  else: 	CALDBhetg/hetgD1996-11-01greffpr001N0006.fits
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;
;history
;	modified from MKACISGARF (VK; MMVjan)
;	brought defaults up to date (VK; MMIXnov)
;-

;	usage
ok='ok' & np=n_params()
if np eq 0 then ok='Insufficient parameters'
if ok ne 'ok' then begin
  print,'Usage: effar=mkhrciarf(wvlar,order=order,arm=arm,caldb=caldb,$'
  print,'  hrmaea=hrmaea,qefil=qefil,greff=greff,verbose=verbose)'
  print,'  return effective areas for HRC-I/NONE or HRC-I/[HML]EG'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

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
CALDBhrc=dbcal+'/data/chandra/hrc/'
CALDBhetg=dbcal+'/data/chandra/default/greff/'
CALDBletg=dbcal+'/data/chandra/default/greff/'

;	mirror effective area
hrmafil=CALDBhrma+'/axeffa/hrmaD1996-12-20axeffaN0008.fits'
if keyword_set(hrmaea) then begin
  if size(hrmaea,/type) eq 7 then hrmafil=strtrim(hrmaea[0],2)
endif
if vv gt 1 then $
  message,'reading HRMA effective areas from: '+hrmafil,/informational

;	HRC-I QE
hrciqefil=CALDBhrc+'/qe/hrciD1999-07-22qeN0007.fits'
if keyword_set(qefil) then begin
  if size(qefil,/type) eq 7 then hrciqefil=strtrim(qefil[0],2)
endif
if vv gt 1 then $
  message,'reading HRC-I QE from: '+hrciqefil,/informational

;	grating efficiencies
if garm eq 'NONE' then greffil='NONE' else $
  if garm eq 'LEG' then $
  greffil=CALDBletg+'/letgD1996-11-01greffpr001N0005.fits' else $
  greffil=CALDBhetg+'/hetgD1996-11-01greffpr001N0006.fits'
if keyword_set(greff) then begin
  if size(greff,/type) eq 7 then greffil=strtrim(greff[0],2)
endif
if vv gt 1 then $
  message,'reading grating efficiencies from: '+greffil,/informational

;{	initialize

legmin=3.0 & legmax=180.0	;LEG 0th order wavelength range
legbin=0.0075			;LEG wavelength grid spacing
megmin=1.0 & megmax=42.0	;MEG 0th order wavelength range
megbin=0.005			;MEG wavelength grid spacing
hegmin=0.5 & hegmax=21.0	;HEG 0th order wavelength range
hegbin=0.0025			;HEG wavelength grid spacing
allmin=legmin < megmin < hegmin
allmax=hegmax > megmax > legmax
allbin=legbin < megbin < hegbin
if keyword_set(gord) then begin
  legmin=legmin/abs(gord) & legmax=legmax/abs(gord) & legbin=legbin/abs(gord)
  megmin=megmin/abs(gord) & megmax=megmax/abs(gord) & megbin=megbin/abs(gord)
  hegmin=hegmin/abs(gord) & hegmax=hegmax/abs(gord) & hegbin=hegbin/abs(gord)
endif

;	reference wavelength grid
nall=(allmax-allmin)/allbin & wall=findgen(nall)*allbin+allmin
nleg=(legmax-legmin)/legbin & wleg=findgen(nleg)*legbin+legmin
nmeg=(megmax-megmin)/megbin & wmeg=findgen(nmeg)*megbin+megmin
nheg=(hegmax-hegmin)/hegbin & wheg=findgen(nheg)*hegbin+hegmin

;;	chip gap locations
;aciss_gaps,hgap,mgap,lgap,onchip=onchip,offset=srcpos,order=gord
;;input OFFSET [arcmin], outputs ONCHIP [S-array#], ?GAP [Ang]

case garm of
  'NONE': begin
    ngrid=nall & wgrid=wall ;& gaps=lall
  end
  'LEG': begin
    ngrid=nleg & wgrid=wleg ;& gaps=lgap
  end
  'MEG': begin
    ngrid=nmeg & wgrid=wmeg ;& gaps=mgap
  end
  'HEG': begin
    ngrid=nheg & wgrid=wheg ;& gaps=hgap
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
tmp=findfile(hrciqefil,count=nhrci)
if greffil ne 'NONE' then tmp=findfile(greffil,count=ngref) else ngref=0
if nhrma eq 0 then ok=hrmafil+': HRMA effective area file does not exist'
if nhrci eq 0 then ok=hrciqefil+': ACIS QE file does not exist'
if ngref eq 0 and greffil ne 'NONE' then $
	ok=greffil+': grating efficiency file does not exist'
if nhrma gt 1 then ok=hrmafil+': too many files'
if nhrci gt 1 then ok=hrciqefil+': too many files'
if ngref gt 1 then ok=greffil+': too many files'
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
;	HRC-I QE
detnam=strtrim(sxpar(hdr,'DETNAM'),2)
for i=1,1 do begin
  tmp=mrdfits(hrciqefil,i,hdr)
  detnam=strtrim(sxpar(hdr,'DETNAM'),2)
  case detnam of
    'HRC-I': hrciqe=tmp
    ;'HRC-S1': begin
    ;  hrcs102qe=create_struct('QE',reform((tmp.QE)[*,0]),$
    ;	'ENERGY',reform((tmp.ENERGY)[*,0]))
    ;  hrcs101qe=create_struct('QE',reform((tmp.QE)[*,1]),$
    ;	'ENERGY',reform((tmp.ENERGY)[*,1]))
    ;end
    ;'HRC-S2': begin
    ;  hrcs204qe=create_struct('QE',reform((tmp.QE)[*,0]),$
    ;	'ENERGY',reform((tmp.ENERGY)[*,0]))
    ;  hrcs201qe=create_struct('QE',reform((tmp.QE)[*,1]),$
    ;	'ENERGY',reform((tmp.ENERGY)[*,1]))
    ;  hrcs203qe=create_struct('QE',reform((tmp.QE)[*,2]),$
    ;	'ENERGY',reform((tmp.ENERGY)[*,2]))
    ;  hrcs202qe=create_struct('QE',reform((tmp.QE)[*,3]),$
    ;	'ENERGY',reform((tmp.ENERGY)[*,3]))
    ;end
    ;'HRC-S3': begin
    ;  hrcs302qe=create_struct('QE',reform((tmp.QE)[*,0]),$
    ;	'ENERGY',reform((tmp.ENERGY)[*,0]))
    ;  hrcs301qe=create_struct('QE',reform((tmp.QE)[*,1]),$
    ;	'ENERGY',reform((tmp.ENERGY)[*,1]))
    ;end
    else: begin
      message,hrciqefil+': could not read HRC-I QE',/informational
      if vv gt 100 then stop,'HALTing.  type .CON to continue'
      return,0*wgrid
    endelse
  endcase
endfor
;	grating efficiencies
if greffil ne 'NONE' then begin
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
endif

;	interpolate everything onto the same wavelength grid
for i=1,4 do begin
  ;HRMA: ENERG_LO,ENERG_HI,EFFAREA
  ;GREFF: ENERGY, EFF[--..-1,0,+1,..++]
  case i of
    1: begin & hrma=hrma1000 & if n_tags(gref1000) gt 0 then gref=gref1000 & end
    2: begin & hrma=hrma0100 & if n_tags(gref1000) gt 0 then gref=gref0100 & end
    3: begin & hrma=hrma0010 & if n_tags(gref1000) gt 0 then gref=gref0010 & end
    4: begin & hrma=hrma0001 & if n_tags(gref1000) gt 0 then gref=gref0001 & end
  endcase
  egrid=0.5*(hrma.energ_lo+hrma.energ_hi)
  h_ea=interpol(hrma.effarea,12.3985/egrid,wgrid)
  if n_tags(gref) gt 0 then begin
    iord=(size(gref.EFF))[1]/2 & oord=lindgen(2L*iord+1L)-iord
    if abs(gord) gt iord then begin
      message,'Grating efficiencies only present up to order '+$
    	strtrim(iord,2),/informational
      return,0*wgrid
    endif
    kord=gord+iord
    gr=interpol(gref.eff[kord,*],12.3985/gref.energy,wgrid)
  endif else gr=0.*wgrid+1.0
  case i of
    1: begin & hea1000=h_ea & gr1000=gr & end
    2: begin & hea0100=h_ea & gr0100=gr & end
    3: begin & hea0010=h_ea & gr0010=gr & end
    4: begin & hea0001=h_ea & gr0001=gr & end
  endcase
endfor
for i=1,1 do begin
  ;QE: ENERGY,QE
  case i of
    1: hiqe=interpol(hrciqe.QE,12.3985/hrciqe.energy,wgrid)
    ;2: h102qe=interpol(hrcs102qe.QE,12.3985/hrcs102qe.energy,wgrid)
    ;3: h101qe=interpol(hrcs101qe.QE,12.3985/hrcs101qe.energy,wgrid)
    ;4: h204qe=interpol(hrcs204qe.QE,12.3985/hrcs204qe.energy,wgrid)
    ;5: h201qe=interpol(hrcs201qe.QE,12.3985/hrcs201qe.energy,wgrid)
    ;6: h203qe=interpol(hrcs203qe.QE,12.3985/hrcs203qe.energy,wgrid)
    ;7: h202qe=interpol(hrcs202qe.QE,12.3985/hrcs202qe.energy,wgrid)
    ;8: h302qe=interpol(hrcs302qe.QE,12.3985/hrcs302qe.energy,wgrid)
    ;9: h301qe=interpol(hrcs301qe.QE,12.3985/hrcs301qe.energy,wgrid)
  endcase
endfor

if vv gt 200 then stop,'HALTing, type .CON to continue'

;	compute effective area for each chip over specified wavelength grid
;	EA_CHIP = CHIP_QE * Sum_SHELL (HRMA_EA*GRATING_EFFICIENCIES)


case garm of
  'MEG': sumshell=hea1000*gr1000+hea0100*gr0100
  'HEG': sumshell=hea0010*gr0010+hea0001*gr0001
  else: sumshell=hea1000*gr1000+hea0100*gr0100+hea0010*gr0010+hea0001*gr0001
endcase
garf=hiqe*sumshell > 0

if vv gt 200 then stop,'HALTing, type .CON to continue'

return,garf
end
