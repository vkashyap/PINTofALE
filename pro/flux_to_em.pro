function flux_to_em,emis,flux=flux,logT=logT,wvl=wvl,Z=Z,NH=NH,defEM=defEM,$
	noph=noph,thresh=thresh, _extra=e
;+
;function	flux_to_em
;	return emission measures (-1 where undeterminable) at various
;	temperatures consistent with given line fluxes and line emissivities.
;
;	for each temperature, assume a delta-function emission measure,
;	compute line fluxes seen through some instrument, account for
;	interstellar absorption, and scale to match the observed flux.
;
;syntax
;	em=flux_to_em(emis,flux=flux,logT=logT,wvl=wvl,Z=Z,NH=NH,$
;	defEM=defEM,noph=noph,thresh=thresh, abund=abund,/temp,/ikev,$
;	effar=effar,wvlar=wvlar,fh2=fh2,he1=he1,heII=heII,/fano)
;
;parameters
;	emis	[INPUT; required] 2D array of line emissivities, EMIS(LOGT,WVL)
;		* if 1-D, assumed to be EMIS(LOGT)
;
;keywords
;	flux	[INPUT] observed fluxes
;		* if not given, assumed to be 1 [(ph|erg)/s/cm^2]
;		* size is forced to match 2nd dimension of EMIS
;		  if <, filled with 1s
;		  if >, extra elements are ignored
;		* NOTE: the author realizes that there may be multiple IDs
;		  for the single observed line, and directs the inquirer
;		  to the function ID_TO_EMIS, which returns the appropriately
;		  concatenated emissivity table.
;	logT	[INPUT] log_10(Temperature [K]) at which EMIS is defined
;		* if not given, assumed to go from 4 to 8 in even steps
;		* interpolated if size does not match 1st dimension of EMIS
;	wvl	[INPUT] wavelengths at which EMIS is defined
;		* if not given or if size does not match 2nd dimension of EMIS,
;		  (a) interstellar absorption correction is not made
;		  (b) all fluxes >gotta be< [ergs/s/cm^2]
;	Z	[INPUT] atomic numbers
;		* if not given, assumed to be 1s (i.e., H. i.e., no
;		  abundance corrections are made)
;		* size is forced to match 2nd dimension of EMIS
;		  if <, filled with 1s
;		  if >, extra elements are ignored
;	NH	[INPUT] H column density [cm^-2]
;		* ignored if WVL is incorrect
;	defEM	[INPUT; default=1e14] default EM to use prior to scaling [cm^-5]
;	noph	[INPUT] if set, FLUX is assumed to be in [ergs/s/cm^2], and
;		no attempts are made to convert things to [ph/...]
;		* forcibly set if WVL is incorrect
;	thresh	[INPUT; default=1e-2] dynamic range in output curves,
;		relative to maximum in each curve
;	_extra	[INPUT ONLY] use this to pass defined keywords to subroutines
;		LINEFLX [ABUND,TEMP,IKEV,EFFAR,WVLAR]
;		ISMTAU [FH2,HE1,HEII,FANO]
;
;restrictions
;	* EFFAR and WVLAR must be correctly defined, or else the output
;	  will be garbage
;	* requires subroutines
;	  LINEFLX
;	      WHEE
;	  ISMTAU
;
;history
;	vinay kashyap (Oct98)
;	bug correction re 0 EMIS (VK; Nov98)
;	bug correction re 0 NH (VK; AugMM)
;-

;	usage
sze=size(emis) & nsze=n_elements(sze) & n=sze(nsze-2)
if n eq 0 then begin
  print,'Usage: em=flux_to_em(emis,flux=flux,logT=logT,wvl=wvl,Z=Z,NH=NH,$'
  print,'       defEM=defEM,noph=noph,thresh=thresh;'
  print,'       LINEFLX:abund,temp,ikev,effar,wvlar; ISMTAU:fh2,he1,heII,fano)'
  return,-1L
endif

;	consistency checks
edim=sze(0)		;dimension of EMIS
if edim eq 0 then begin & mt=1L & mw=1L & endif		;scalar
if edim eq 1 then begin & mt=sze(1) & mw=1L & endif	;EMIS=EMIS(LOGT)
if edim eq 2 then begin & mt=sze(1) & mw=sze(2) & endif	;EMIS=EMIS(LOGT,WVL)
if edim gt 2 then begin
  message,'emissivity array is too complex to understand',/info
  return,0*emis
endif
nf=n_elements(flux) & nt=n_elements(logT) & nw=n_elements(wvl) & nz=n_elements(Z)
;
fx=fltarr(mw)+1.				;fluxes
if nf gt 0 and nf lt mw then fx(0:nf-1)=flux
if nf ge mw then fx=flux(0:mw-1)
;
tlog=interpol([4.,8.],mt)			;temperatures
if nt gt 0 and nt ne mt then tlog=interpol(logt,mt)
if nt eq mt then tlog=logt
;
if keyword_set(NH) then setabs=1 else setabs=0	;include ISM absorption?
;
ww=findgen(mw) 					;wavelengths
if nw eq 0 or nw ne mw then begin
  setnoph=1 & setabs=0
endif else ww=wvl
;
zz=intarr(mw)+1					;atomic numbers
if nz gt 0 and nz lt mw then zz(0:nf-1)=Z
if nz ge mw then zz=Z(0:mw-1)
;
if keyword_set(noph) then setnoph=1		;everything in [ergs/...]
;
if not keyword_set(defEM) then dem=1d14 else dem=double(defEM)
						;default EM [cm^-5]
if not keyword_set(thresh) then thresh=1e-2	;dynamic range

;	get optical depths
tau=fltarr(mw)
if keyword_set(setabs) then begin
  tau=ismtau(ww,NH=NH,fH2=0,/Fano,_extra=e) < 69.
endif

;	get the predicted fluxes
ct=dblarr(nt,nw)
for iw=0,nw-1 do begin
  if edim eq 0 then line=[emis]
  if edim eq 1 then line=[emis(*)]
  if edim eq 2 then line=reform(emis(*,iw))
  for it=0,nt-1 do begin
    tmp=lineflx(line(it),tlog(it),ww(iw),zz(iw),noph=setnoph,dem=dem, _extra=e)
    ct(it,iw)=tmp*exp(-tau(iw))
  endfor
endfor

;	compare to observed and scale EM
em=dblarr(nt,nw)-1
for iw=0,nw-1 do begin
  if edim eq 0 then begin
    pct=[ct] & line=[emis]
  endif
  if edim eq 1 then begin
    pct=[ct(*)] & line=[emis(*)]
  endif
  if edim eq 2 then begin
    pct=reform(ct(*,iw)) & line=reform(emis(*,iw))
  endif
  oo=where(pct ge thresh*max(pct),moo)
  if max(pct) le 0 then moo=0L
  if moo gt 0 then begin
    frac=dem*fx(iw)/pct(oo)
    em(oo,iw)=frac(*)
  endif
endfor

return,em
end
