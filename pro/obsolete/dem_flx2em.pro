function dem_flx2em,fx,wvl,Z,ion,eden=eden, _extra=e
;+
;function	dem_flx2em
;	returns an array of Emission Measures [cm^-5] (or [cm^-3] if
;	keywords EFFAR and WVLAR are set) computed at various densities.
;
;	read in the line emissivities for the matched line at specified
;	density, fold in ion balance, assume delta-function EMs at each
;	temperature, compute line fluxes, account for ISM absorption,
;	scale to obtain observed flux, determine Pottasch point for
;	each match, repeat for a different electron density.
;
;	if more than one line in the database is matched to the observed
;	line, then a separate set of EM(T) is calculated for each match
;	(the relative values of which are determined by the relative values
;	of the ion-balanced, abundance included, emissivities) and are
;	then averaged over the matches at each temperature.
;
;	returns -1 where undefinable
;
;parameters
;	fx	[INPUT; required] measured flux in spectral line [counts/s]
;	wvl	[INPUT; required] wavelength(s) of matching spectral line [Ang]
;	Z	[INPUT; required] atomic numbers of matching line(s)
;	ion	[INPUT; required] ionic states of matching line(s)
;
;keywords
;	eden	[INPUT] electron densities [cm^-3] at which EMs are to
;		be computed
;		* default: [5e7, 1e8, 5e8, ... 1e14]
;	_extra	[INPUT] allows setting of defined keywords to
;		FLX2EM [ABUND, NH, SKIPRD, DEFEM]
;		RD_LINE [DBDIR, ALLAH]
;		FOLD_IONEQ [CHIFIL [EQFILE, CHIDIR]]
;		LINEFLX [EFFAR, WVLAR]
;		ISMTAU [FH2, HE1, HEII, FANO]
;		POTTASCH [LEVEL]
;
;restrictions
;	* requires LINEID to have been run first
;	* requires constant-density line emissivity database to be present
;	* requires subroutines FLX2EM, RD_LINE, FOLD_IONEQ, GETABUND,
;	  LINEFLX, ISMTAU, POTTASCH, RD_IONEQ, READ_IONEQ, SYMB2ZION,
;	  LAT2ARAB, KILROY, WHEE
;
;history
;	vinay kashyap (Jun97)
;-

message,'OBSOLETE!',/informational

;	check input for problems
ok=''			;assume everything is OK
nf=n_elements(fx) & nw=n_elements(wvl) & nz=n_elements(Z) & ni=n_elements(ion)
if nf eq 0 then ok='Missing flux' else $
  if nw eq 0 then ok='Missing wavelength' else $
    if nz eq 0 then ok='Missing atomic number' else $
      if ni eq 0 then ok='Missing ionic state' else $
        if (nw ne nz) or (nw ne ni) or (nz ne ni) then ok='mismatched inputs'

;	usage
if ok ne '' then begin
  message,ok,/info
  print,''
  print,'Usage: em=dem_flx2em(fx,wvl,Z,ion,lfx,eden=eden)'
  print,'  compute Pottasch-corrected EMs at various densities'
  print,'  also accepts defined keywords ABUND,NH,SKIPRD,DEFEM (FLX2EM);'
  print,'  DBDIR,ALLAH (RD_LINE); CHIFIL,EQFILE,CHIDIR (FOLD_IONEQ);'
  print,'  EFFAR,WVLAR (LINEFLX); FH2,HE1,HEII,FANO (ISMTAU); LEVEL (POTTASCH)'
  return,-1L
endif

;	at which densities?
if not keyword_set(eden) then $
  eden=[5e7,1e8,5e8,1e9,5e9,1e10,5e10,1e11,5e11,1e12,5e12,1e13,5e13,1e14]
nden=n_elements(eden)
emden=dblarr(nden)-1.		;initialize output array

;	call FLX2EM
for i=0,nden-1 do begin
  n_e=eden(i)			;electron density [cm^-3]

  lfx=0
  em=flx2em(fx,wvl,Z,ion,lfx,logT=logT,n_e=n_e, _extra=e)	;get EM(T)
  pot=pottasch(lfx,logT,Z=Z,ion=ion, _extra=e)	;Pottasch correction factors

  emm=dblarr(nw)				;corrected EMs
  for iw=0,nw-1 do begin
    oo=where(em(*,iw) gt 0,moo)
    if moo gt 0 then emm(iw)=min(em(oo,iw))/pot(iw) else emm(iw)=-1.
  endfor

  oo=where(emm gt 0,moo)
  if moo gt 0 then emden(i)=total(emm(oo))/moo		;averaged EM
endfor

return,emden
end
