function chandra_psfsize,skyx,skyy,x0=x0,y0=y0,eefrac=eefrac,energy=energy,$
	instrum=instrum,psftbl=psftbl,defoc=defoc,verbose=verbose,$
	_extra=e
;+
;function	chandra_psfsize
;	read in the Chandra PSF size table and return the radius [arcsec]
;	that encloses a specified energy of the PSF 
;
;syntax
;	psfsize=chandra_psfsize(skyx,skyy,eefrac=eefrac,energy=energy,$
;	instrum=instrum,psftbl=psftbl,verbose=verbose)
;
;parameters
;	skyx	[INPUT; required] the x pixel locations of the points at
;		which to compute the size of the PSF
;		* or off-axis angles in [arcmin] if SKYY is illegal
;	skyy	[INPUT] the y pixel locations corresponding to SKYX
;		* sizes of SKYX and SKYY _must_ match.  if they don't --
;		  -- if N(SKYY)=1, SKYY[0] is expanded out to N(SKYX)
;		  -- if N(SKYY)=0 or N(SKYY)>N(SKYX), then SKYY is ignored
;		     and SKYX is assumed to be the offaxis angle in [arcmin]
;
;keywords
;	x0	[INPUT] x-location of the aimpoint [pixel]
;	y0	[INPUT] y-location of the aimpoint [pixel]
;		* if not given, X0 and Y0 are presumed to be at the
;		  nominal Chandra aimpoints
;		* they are ignored if SKYY is not given
;	eefrac	[INPUT; default=0.9] the fraction of energy enclosed
;	energy	[INPUT; default=1.5 keV] the energy at which to determine
;		the PSF size
;		* note: there is no interpolation -- the entries closest
;		  to the specified value will be chosen.  any interpolation
;		  can always be done post facto by the user.
;	instrum	[INPUT; default='HRC-I'] the detector
;	psftbl	[INPUT; default='/soft/ciao/data/psfsize20010416.fits']
;		full path name to the wavdetect-compatible FITS file that
;		contains the PSF size information
;	defoc	[INPUT; default=0] the defocus value [mm] at which to
;		determine the PSF size
;		* note that the default Chandra PSFTBL only lists DEFOC=0
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Apr2006)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(skyx) & ny=n_elements(skyy)
if np eq 0 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='SKYX is undefined'
if ok ne 'ok' then begin
  print,'Usage: psfsize=chandra_psfsize(skyx,skyy,eefrac=eefrac,energy=energy,$'
  print,'       instrum=instrum,psftbl=psftbl,defoc=defoc,verbose=verbose)'
  print,'  read in Chandra PSF-size table and return the appropriate values'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
frac=0.9
if keyword_set(eefrac) then begin
  frac=abs(eefrac[0])
  if frac ge 1 and frac lt 100 then frac=frac/100.
  if frac ge 100 then frac=1.D - 1./frac
endif
;
nrg=1.5 & if keyword_set(energy) then nrg=float(energy[0])
;
focusdep=0. & if keyword_set(defoc) then focusdep=float(defoc[0])
;
psftab='/soft/ciao/data/psfsize20010416.fits'
if keyword_set(psftbl) then begin
  szp=size(psftbl) & nszp=n_elements(szp)
  if szp[nszp-2L] eq 7 then psftab=psftbl[0] else message,$
	'PSFTBL cannot be understood; using default',/informational
endif
;
det='HRC-I'
if keyword_set(instrum) then begin
  cc=strupcase(instrum)
  if strpos(cc,'HRC-I') ge 0 and strpos(cc,'L') lt 0 then det='HRC-I'
  if strpos(cc,'HRC-I') ge 0 and strpos(cc,'L') ge 0 then det='HRC-I/LEG'
  ;
  if strpos(cc,'HRC-S') ge 0 and $
     (strpos(cc,'L') lt 0 and strpos(cc,'M') lt 0 and strpos(cc,'HEG') lt 0) $
     then det='HRC-S'
  if strpos(cc,'HRC-S') ge 0 and strpos(cc,'L') ge 0 then det='HRC-S/LEG'
  if strpos(cc,'HRC-S') ge 0 and strpos(cc,'M') ge 0 then det='HRC-S/MEG'
  if strpos(cc,'HRC-S') ge 0 and strpos(cc,'HEG') ge 0 then det='HRC-S/HEG'
  ;
  if strpos(cc,'ACIS-I') ge 0 then det='ACIS-I'
  ;
  if strpos(cc,'ACIS-S') ge 0 and $
     (strpos(cc,'L') lt 0 and strpos(cc,'M') lt 0 and strpos(cc,'H') lt 0) $
     then det='ACIS-S'
  if strpos(cc,'ACIS-S') ge 0 and strpos(cc,'L') ge 0 then det='ACIS-S/LEG'
  if strpos(cc,'ACIS-S') ge 0 and strpos(cc,'M') ge 0 then det='ACIS-S/MEG'
  if strpos(cc,'ACIS-S') ge 0 and strpos(cc,'H') ge 0 then det='ACIS-S/HEG'
  ;
  if strpos(cc,'HRMA') ge 0 and $
     (strpos(cc,'L') lt 0 and strpos(cc,'M') lt 0 and strpos(cc,'HEG') lt 0) $
     then det='HRMA'
  if strpos(cc,'HRMA') ge 0 and strpos(cc,'L') ge 0 then det='HRMA/LEG'
  if strpos(cc,'HRMA') ge 0 and strpos(cc,'M') ge 0 then det='HRMA/MEG'
  if strpos(cc,'HRMA') ge 0 and strpos(cc,'HEG') ge 0 then det='HRMA/HEG'
endif
if strpos(det,'ACIS') ge 0 then pixsiz=0.4920 else $  ;arcsec -- for ACIS
 if strpos(det,'HRC') ge 0 then pixsiz=0.13175 else $ ;arcsec -- for HRC
  pixsiz=60.	;input assumed to be in [arcmin]
;
xnom=0. & ynom=0.
if strpos(det,'HRC-I') ge 0 then begin
  xnom=16384. & ynom=16384.
endif
if strpos(det,'HRC-S') ge 0 then begin
  xnom=32768. & ynom=32768.
endif
if strpos(det,'ACIS-I') ge 0 then begin
  xnom=4096. & ynom=4096.
endif
if strpos(det,'ACIS-S') ge 0 then begin
  xnom=4096. & ynom=4096.
endif
;
xx=[skyx[*]] & yy=0*xx+ynom
if ny eq 1 then yy[*]=skyy[0]
if ny eq nx then yy=[skyy[*]]
offx=sqrt((xx-xnom)^2+(yy-ynom)^2)*pixsiz/60.	;convert [pixel] to [arcmin]
if nx ne ny then offx=xx	;input assumed to be in [arcmin]

;	read in PSFTBL
idxtbl=mrdfits(psftab,1,hdridx)
table=idxtbl.table
case det of
  'HRC-I': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'HRC-I' and $
	strtrim(idxtbl.GRATING,2) eq 'NONE',modet)
  'HRC-S': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'HRC-S' and $
	strtrim(idxtbl.GRATING,2) eq 'NONE',modet)
  'ACIS-I': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'ACIS-I' and $
	strtrim(idxtbl.GRATING,2) eq 'NONE',modet)
  'ACIS-S': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'ACIS-S' and $
	strtrim(idxtbl.GRATING,2) eq 'NONE',modet)
  'HRC-S/LEG': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'HRC-S' and $
	strtrim(idxtbl.GRATING,2) eq 'LEG',modet)
  'HRC-S/HEG': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'HRC-S' and $
	strtrim(idxtbl.GRATING,2) eq 'HEG',modet)
  'HRC-S/MEG': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'HRC-S' and $
	strtrim(idxtbl.GRATING,2) eq 'MEG',modet)
  'ACIS-S/LEG': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'ACIS-S' and $
	strtrim(idxtbl.GRATING,2) eq 'LEG',modet)
  'ACIS-S/HEG': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'ACIS-S' and $
	strtrim(idxtbl.GRATING,2) eq 'HEG',modet)
  'ACIS-S/MEG': odet=where(strtrim(idxtbl.DETECTOR,2) eq 'ACIS-S' and $
	strtrim(idxtbl.GRATING,2) eq 'MEG',modet)
  else: begin
    message,det+': not implemented',/informational
    if vv gt 1000 then stop,'HALTing; type .CON to continue'
    return,-1L
  end
endcase
if modet gt 0 then extdet=strtrim(table(odet[0]),2) else begin
  message,det+': not found in '+psftab,/informational
  if vv gt 1000 then stop,'HALTing; type .CON to continue'
  return,-1L
endelse

;	find the appropriate linking extension
kext=2L & skipext=0L & go_on=1
while go_on do begin
  tblk=mrdfits(psftab,kext,hdrk)
  if not keyword_set(tblk) then begin
    message,psftab+': does not contain extension '+extdet,/informational
    if vv gt 1000 then stop,'HALTing; type .CON to continue'
    return,-1L
  endif
  if strtrim(sxpar(hdrk,'EXTNAME'),2) eq extdet then go_on=0 else begin
    skipext=skipext+n_elements(tblk.(0))
  endelse
  kext=kext+1L
endwhile
defocus=tblk.DEFOCUS
keV=tblk.ENERGY
theta=tblk.THETA
table3=strtrim(tblk.TABLE,2)

;	these are the extensions that will contain the tables
;	of interest
ok=where(abs(defocus-focusdep) eq min(abs(defocus-focusdep)) and $
	abs(keV-nrg) eq min(abs(keV-nrg)),mok)
if mok eq 0 then begin
  message,psftab+': has nothing of interest',/informational
  if vv gt 1000 then stop,'HALTing; type .CON to continue'
  return,-1L
endif
tht=theta[ok] & tab3=table3[ok]
rr=fltarr(mok) & ii=intarr(mok)

;	read the EEfrac(radius) grid
kext=kext+skipext & go_on=1
while go_on do begin
  tbl=mrdfits(psftab,kext,hdr)
  if not keyword_set(tbl) then begin
    message,psftab+': is incomplete',/informational
    if vv gt 1000 then stop,'HALTing; type .CON to continue'
    return,-1L
  endif
  extnam=strtrim(sxpar(hdr,'EXTNAME'),2)
  if vv gt 5 then print,extnam,kext
  oj=where(extnam eq tab3,moj)
  if moj eq 1 then begin
    ii[oj[0]]=1
    rr[oj[0]]=interpol(tbl.RADIUS,tbl.FRACTION,frac)>0
  endif
  kext=kext+1L
  if total(ii) eq mok then go_on=0
endwhile

;	output
psfsize=interpol(rr,tht,offx)>0

return,psfsize
end
