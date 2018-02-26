function flx2em,fx,wvl,z,ion,lfx,logT=logT,abund=abund,NH=NH,$
	skiprd=skiprd,defEM=defEM, _extra=e
;+
;function	flx2em
;	generate a set of emission measures at various temperatures
;	that are consistent with the given flux and given line ID(s)
;
;	read in the line emissivities for the matched line, fold in
;	ion balance.  For each temperature, assume a delta-function
;	EM, compute line fluxes seen through some instrument, account
;	for interstellar absorption, and scale this EM to obtain the
;	observed flux.
;
;	if more than 1 line is matched, then a separate set of EMs
;	is generated for each match, the relative values of which
;	is determined by the relative values of the ion-balanced,
;	abundance included emissivities.
;
;	returns -1 where undefinable
;
;parameters
;	fx	[INPUT; required] measured flux in spectral line [counts/s]
;	wvl	[INPUT; required] wavelength(s) of matching spectral line [Ang]
;	Z	[INPUT; required] atomic numbers of matching line(s)
;	ion	[INPUT; required] ionic states of matching line(s)
;	lfx	[I/O] if set on input as a 2D array (N_LOGT,N_WVL), call to
;		RD_LINE and FOLD_IONEQ is skipped.
;		on output contains all the line emissivities of interest.
;		* ignored if LOGT is not set or if size of 2nd dimension
;		  does not match LOGT
;		* ignored if keyword SKIPRD is not set
;
;keywords
;	logT	[OUTPUT] temperatures at which EM is computed -- this
;		is read in along with line emissivities.
;	abund	[INPUT] abundances (default: Anders & Grevesse)
;	NH	[INPUT] H column density [cm^-2]
;	skiprd	[INPUT] if set AND LFX is properly defined, skips call
;		to RD_LINE
;	defEM	[INPUT; default=1e14] default EM to use prior to scaling [cm^-5]
;
;	_extra	[INPUT] allows setting defined keywords to
;		RD_LINE [LOGP, N_E, DBDIR, ALLAH]
;		FOLD_IONEQ [CHIFIL [EQFILE, CHIDIR]]
;		LINEFLX [EFFAR, WVLAR]
;		ISMTAU [FH2,HE1,HEII,FANO]
;
;restrictions
;	* requires LINEID to have been run first
;	* EFFAR and WVLAR must be correctly defined, or else the output
;	  will be garbage
;	* requires subroutines
;	  RD_LINE [KILROY, SYMB2ZION [LAT2ARAB]]
;	  FOLD_IONEQ [WHEE, GETLOCMAX, RD_IONEQ [READ_IONEQ]]
;	  GETABUND
;	  LINEFLX [WHEE]
;	  ISMTAU
;
;history
;	vinay kashyap (Feb97)
;	added keyword DEFEM; restructured EM loop (VK; Apr97)
;	sum up emissivities if many lines at same wavelength (VK; Jul97)
;	added message saying routine is obsolete (VK; Oct98)
;-

message,'This routine is obsolete.  Use ID_TO_EMIS and FLUX_TO_EM instead',/info

;	usage
np=n_params()
if np lt 4 then begin
  print,'Usage: em=flx2em(fx,wvl,Z,ion,lfx,logT=logT,abund=abund,NH=NH,$
  print,'  defEM=defEM,/skiprd)'
  print,'  return emission measures at various temperatures consistent with'
  print,'  observed line flux'
  print,'also accepts defined keywords:'
  print,'  LOGP,DBDIR,ALLAH (RD_LINE); EFFAR,WVLAR (LINEFLX)'
  print,'  CHIFIL,EQFILE,CHIDIR (FOLD_IONEQ) and FH2,HE1,HEII,FANO (ISMTAU)'
  return,-1L
endif

;	inputs
nf=n_elements(fx) & nw=n_elements(wvl) & nz=n_elements(z) & nion=n_elements(ion)
szf=size(lfx) & mt=szf(1) & mw=szf(2) & nt=n_elements(logt)
nab=n_elements(abund)

;	catch trivial errors
c1='y' & c2='y'
if nf eq 0 then c1='Observed flux not given'
if nf gt 1 then c1='Cannot handle more than one flux value at a time'
if nw eq 0 then c1='Matching wavelengths not given'
if nz gt 1 and nz ne nw then c1='Atomic numbers not correctly specified'
if nion gt 1 and nion ne nw then c1='Ionic states not correctly specified'
if strmid(c1,0,1) ne 'y' then begin
  message,c1,/info & return,-1L				;input errors
endif
if szf(0) ne 0 then begin
  if mw ne nw then c2='Line emissivities incorrect; will be overwritten'
  if nt eq 0 then c2='Temperatures must be supplied; will read in '+$
	'both T and line emissivities'
  if nt ge 1 and nt ne mt then c2='Temperature v/s line emissivity mismatch;'+$
	' emissivities will be overwritten'
  if nt gt 1 and szf(0) gt 2 then c2='T v/s line emissivity mismatch;'+$
	' emissivities will be overwritten'
  if strmid(c2,0,1) ne 'y' then message,c2,/info	;non-fatal mistakes
endif
if keyword_set(skiprd) and strmid(c2,0,1) ne 'y' then begin	;oy vey
  print,'type q to return instanter, z to halt precipitously'
  wait,1 & c1=get_kbrd(0)
  if c1 eq 'q' then return,-1L
  if c1 eq 'z' then stop
endif
if strmid(c2,0,1) ne 'y' then skiprd=0
if nab ne 30 then abund=getabund('anders & grevesse')

;	initialize, store, etc.
ff=[fx] & ww=[wvl] & zz=intarr(nw) & jon=zz
if nz eq 1 then zz(*)=z(0) else zz=z
if nion eq 1 then jon(*)=ion(0) else jon=ion
rom=['I','II','III','IV','V','VI','VII','VIII','IX','X']
rom=[rom,'X'+rom,'XX'+rom,'XXX'+rom]	;roman numerals from 1-40
atom=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',$
	'Ni','Cu','Zn']			;elements from 1-30
if not keyword_set(defEM) then dem=1d14 else dem=double(defEM)
					;default EM [cm^-5]

;	read in the line emissivities and fold in ion balance
if not keyword_set(skiprd) then begin
  lfx=[-1.]
  for iw=0,nw-1 do begin
    elem=atom(zz(iw)-1)+rom(jon(iw)-1) & wr=abs(ww(iw))
    tmp=rd_line(elem,wrange=wr,logT=logT,Z=tZ,_extra=e)	;emissivities
    if n_elements(tZ) gt 1 then begin
      tmp2=reform(tmp(*,0))
      for i=1,n_elements(tZ)-1 do tmp2(*,0)=tmp2(*,0)+tmp(*,i)
      tmp=tmp2
    endif
    lfx=[lfx,tmp] & nt=n_elements(logT) & mt=nt
  endfor
  lfx=lfx(1:*) & mf=n_elements(lfx)
  if mf ne nw*nt then begin
    message,'HALT: Some IDs appear to be missing'
  endif else lfx=reform(lfx,nt,nw)
  ;
  lfx=fold_ioneq(lfx,zz,jon,logt=logt, _extra=e)	;ion balance
  ;
  ;smooth result of ion balance to eliminate stepped look
  for i=0,nw-1 do begin
    tmp=reform(lfx(*,i)) & if nt gt 41 then tmp=smooth(tmp,3) & lfx(*,i)=tmp
  endfor
endif

;	get optical depths
tau=0*ww
if keyword_set(NH) then tau=ismtau(ww,NH=NH,fH2=0,/Fano,_extra=e) < 69.

;	now to get counts...
em=dblarr(nt,nw)-1 & ct=dblarr(nt,nw) & fia=fltarr(nt,nw)
for it=0,nt-1 do begin				;{at each temperature
  fi=reform(lfx(it,*))			;emissivities w. ion balance
  fia(it,*)=fi*abund(zz-1)		;and abundances
  tfia=total(fia(it,*)) & if tfia le 0 then tfia=1.
  fia(it,*)=fia(it,*)/tfia		;relative contribution of each match
  for iw=0,nw-1 do begin			;{at each matching wvl
    if fi(iw) gt 0 then $
      ct(it,iw)=lineflx(fi(iw),logt(it),ww(iw),zz(iw),$
	abund=abund,dem=dem,_extra=e)			;[ph/s]
  endfor					;iw=0,nw-1}
endfor						;it=0,nt-1}

;	... and emission measures
for iw=0,nw-1 do begin
  oo=where(ct(*,iw) gt 1e-2*max(ct(*,iw)),moo)
  if moo gt 0 then begin
    tmp=reform(dem*fx(0)*fia(oo,iw)/ct(oo,iw)/exp(-tau(iw)))
	;default_EM*(observed_counts*fraction_in_line)/(model_counts)/absorp.
    em(oo,iw)=tmp(*)		;is this because of an IDL bug?
				;em(oo,iw)=tmp(oo),tmp(oo,iw),etc don't work!
  endif
endfor

return,em
end
