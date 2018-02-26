function fold_aria,x,spec,arstr=arstr,effar=effar,wvlar=wvlar,lines=lines,$
	unfold=unfold, _extra=e
;+
;function	fold_aria
;	returns spectrum with effective area folded in
;
;syntax
;	aspec=fold_aria(x,spec,arstr=arstr,effar=effar,wvlar=wvlar,$
;	/lines,/unfold,/perbin)
;
;parameters
;	x	[INPUT; required] wavelengths at which spectrum is given
;	spec	[INPUT] spectrum, SPEC[X]
;		* size must match that of WVL, unless
;		* if scalar, expanded out to match size of X, or
;		* if missing, assumed to be unity
;
;keywords
;	arstr	[INPUT] structure containing effective areas and wavelengths
;		* see ARIA.PRO for definition and construction
;		* if not given, or is not a structure, looks to EFFAR and WVLAR
;	effar	[INPUT] if ARSTR is not defined, take effective areas from
;		this keyword
;		* if missing or undefined, assumed to be unity
;	wvlar	[INPUT] wavelengths at which EFFAR is defined
;		* size MUST match that of EFFAR
;	lines	[INPUT] if set, assumes that SPEC[X] are individual lines
;		and not a spectrum.
;	unfold	[INPUT] if set, *divides* the input spectrum by the
;		effective area
;		* WARNING: where the area is zero, the output is set to NaN
;	perbin	[INPUT] if set, assumes that units of SPEC are [.../bin]
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	* requires IDL5+
;	* requires subroutines
;	  -- ARIA
;	  -- REBINW
;	* X, ARSTR.WVL, and WVLAR must all have the same units
;	* if grating orders > 1 are present in ARSTR
;	  -- X must be a wavelength scale
;	  -- if LINES and/or UNFOLD are set, output will be garbage
;
;history
;	vinay kashyap (JunMM)
;-

;	usage
ok='ok' & nx=n_elements(x) & ns=n_elements(spec) & np=n_params()
if np eq 0 then ok='' else $
 if nx eq 0 then ok='missing wavelengths' else $
  if (ns gt 1 and ((ns ne nx) or (ns ne (nx-1L)))) then ok='mismatched inputs'
if ok ne 'ok' then begin
  print,'Usage: ay=fold_spec(x,y,arstr=arstr,effar=effar,wvlar=wvlar,$'
  print,'       /lines,/unfold,/perbin)'
  print,'  return spectrum with effective area folded in'
  message,ok,/info
  return,-1L
endif

;	check inputs
xx=x & nxx=nx & y=fltarr(nx)+1.
if ns eq 0 then spec=y else $		;SPEC is undefined
 if ns eq 1 then y[*]=spec[0] else $	;SPEC is a scalar
  if ns eq nx then y=spec else $	;SPEC=SPEC[X]
   if ns eq nx-1L then begin
     message,'Inputs assumed to be bin-boundaries and spectrum',/info
     if keyword_set(lines) then begin
       message,'but keyword LINES has been set.',/info
       message,'this cannot be.  obviously something is seriously wrong.',/info
       message,'quitting without further ado',/info
       return,spec
     endif
     xx=0.5*(x[1:*]+x) & nxx=nx-1L & y=spec
   endif

;	check input effective areas
na=n_tags(arstr) & nea=n_elements(effar) & nwa=n_elements(wvlar)
areff=fltarr(nxx)+1. & arwvl=xx
if na eq 0 then begin	;(ARSTR not a structure
  if nea ne 0 and nea eq nwa then begin
    areff=[effar[*]] & arwvl=[wvlar[*]]
  endif
  areae=aria(areff,arwvl,1)
endif else begin	;ARSTR)(but is it a valid structure?
  arnam=tag_names(arstr.(0))
  if arnam[0] ne 'AREA' then begin
    message,'input area structure in unknown format.  ignoring.',/info
    areae=aria(areff,arwvl,1)
  endif else areae=arstr
endelse			;ARSTR)
naa=n_tags(areae)

;	UNFOLD
ok='ok'
if keyword_set(unfold) then begin
  ;	make an "effective" effective area
  areff=fltarr(nxx)
  for i=0L,naa-1L do begin
    a=areae.(i).AREA & w=areae.(i).WVL & o=areae.(i).ORDER
    ww=o[0]*w & amax=max(a) & ow=sort(ww) & ww=ww[ow] & a=a[ow]
    aa=(interpol(a,ww,xx) > 0) < amax
    areff=areff+aa
    if abs(o[0]) gt 1 then ok='orders > 1 are present.  CAVEAT VSVATOR!'
  endfor
  ;	divide SPEC by "effective" effective area
  z=0*y
  oo=where(areff gt 0,moo) & if moo gt 0 then z[oo]=y[oo]/areff[oo]
  oo=where(areff eq 0,moo) & if moo gt 0 then z[oo]=!VALUES.F_NAN
  if ok ne 'ok' then message,ok,/info
  return,z
endif

;	FOLD the LINES
ok='ok'
if keyword_set(lines) then begin
  ;	make a 2D matrix to hold all the folded flux values
  z=fltarr(nxx,naa)
  for i=0L,naa-1L do begin
    a=areae.(i).AREA & w=areae.(i).WVL & o=areae.(i).ORDER & amax=max(a)
    ow=sort(w) & w=w[ow] & a=a[ow]
    aa=(interpol(a,w,xx) > 0) < amax
    z[*,i] = y[*]*aa[i]
    if abs(o[0]) gt 1 then ok='orders > 1 are present.  CAVEAT VSVATOR!'
  endfor
  if ok ne 'ok' then message,ok,/info
  return,z
endif

;	sort, if necessary (this is needed because of call to REBINW)
if xx[1] lt xx[0] then begin
  unsort=1	;make sure the spectrum does gets unsorted later
  ox=sort(xx) & xx=xx[ox] & y=y[ox] & ii=lindgen(nxx) & ii=ii[ox]
endif

;	FOLD the SPECTRUM
z=fltarr(nxx)
for i=0L,naa-1L do begin
  a=areae.(i).AREA & w=areae.(i).WVL & o=areae.(i).ORDER & amax=max(a)
  ow=sort(w) & w=w[ow] & a=a[ow]
  aa=(interpol(a,w,xx) > 0) < amax
  tmp=y*aa & ww=o[0]*x
  if nx eq nxx then begin
    dw=ww[1:*]-ww & wx=[ww-0.5*dw,ww[nx-1L]+0.5*dw[nx-2L]]
  endif else wx=ww
  tmp2=rebinw(tmp,wx,xx,perbin=perbin)
  z=z+tmp2
endfor

;	unsort, if necessary
if keyword_set(unsort) then begin
  oi=sort(ii) & z=z[oi]
endif

return,z
end
