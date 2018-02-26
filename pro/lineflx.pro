function lineflx,line,logT,wvl,Z,DEM=DEM,temp=temp,abund=abund,noph=noph,$
	nhne=nhne,effar=effar,wvlar=wvlar,ikev=ikev,noabund=noabund,$
	regrid=regrid, _extra=e
;+
;function	lineflx
;	computes observeable counts [ct/s] or flux [ergs/s] from spectral
;	line for specified DEM by simple trapezoidal integration over
;	temperatures.
;	if effective area is not specified, output will be [ph/s/cm^2]
;	or [ergs/s/cm^2]
;
;	input line emissivities are in [1e-23 ergs cm^3/s].  these include
;	ionization balance, but not (necessarily) abundances.  the line
;	emissivities are multiplied by the differential emission measure
;	and element abundance to get fluxes, then divided by the line energy
;	to get photon fluxes, and multiplied by the supplied effective areas
;	to get count rates.  no effort is made to determine line profiles
;	or to include spectral responses (beyond the effective areas).
;
;syntax
;	fx=lineflx(line,logT,wvl,Z,DEM=DEM,/temp,abund=abund,/noph,$
;	nhne=nhne,effar=effar,wvlar=wvlar,/ikev,/noabund)
;
;parameters
;	line	[INPUT; required] array of line cooling emissivities
;		in units of 1e-23 erg cm^3/s
;		* if 2D array, LINE==LINE(logT,WVL)
;
;		* WARNING: will return garbage if given 1D array LINE(WVL)
;		  (or even 2D array LINE(1,WVL))
;		  use a for-loop to handle such a case
;
;		* WARNING: will get converted to 1-element vector if input
;		  as scalar
;	logT	[INPUT; required] array of log10(Temperature [K]) at
;		which emissivities are given.
;		* array size MUST match that of LINE
;		* if not regularly gridded, remember to set REGRID!
;	wvl	[INPUT; required] wavelength(s) [Angstrom] of lines at
;		which LINE is given
;		* array size MUST match LINE
;	Z	[INPUT] atomic number of element that generates each WVL
;		* if Z has less elements than WVL, Z is ignored
;		* if Z is to be ignored, all lines are assumed to be from
;		  same element, and abundance values are ignored
;
;keywords
;	DEM	[INPUT] Differential Emission Measure at each T [cm^-5/logK]
;		* default is a constant=1e14!!
;		* if array size does not match that of LOGT, then DEM is
;		  assumed to cover the same range and is linear interpolated
;		* if defined at only 1 temperature, assumed to be EM [cm^-5]
;	temp	[INPUT] if set, assumes that logT is actually in T [K], and
;		that DEM is given in units of [cm^-5/K]
;		* ignored if logT is defined at only one bin
;	abund	[INPUT] abundances relative to H (abund(0)=1)
;		* abund(Z-1) contains the abundance for element Z
;		* if array size is smaller than the largest present Z,
;		  the last element is replicated to fill the gap
;		* default: Anders & Grevesse (1989)
;	noph	[INPUT] if set, does not do the conversion to [ph] from [ergs]
;	nhne	[INPUT] the ratio of N(H)/n_e needed to complete the
;		intensity calculation
;		* if not set, assumes 0.83, which is the limiting value
;		  for high-temperature plasma with approximately cosmic
;		  abundances
;		* if array and size does not match LOGT, then uses only
;		  the first element
;		* if Z is not defined, or is 1, then it is **IGNORED**
;		* see discussion in the documentation of POPSOL()
;	effar	[INPUT] effective area [cm^2]
;		* if not set, assumed to be 1 cm^2
;	wvlar	[INPUT] wavelengths at which effective area is defined.
;		* if not set, assumed to be WVL
;		* array size MUST match that of EFFAR
;	ikev	[INPUT] if set, assumes that WVL and WVLAR are in [keV],
;		NOT [Angstrom]
;	noabund	[INPUT] if set, ABUND values are ignored in the calculations,
;		same as though they (or Z) had been set to 1
;		* but NHNE gets applied
;	regrid	[INPUT] if set, assumes that LOGT is irregularly sampled
;		and regrids it (and LINE) to a regular grid
;		* NOT IMPLEMENTED YET
;	_extra	[INPUT] junk -- here only to prevent crashing the program
;
;subroutines
;	GETABUND
;	WHEE
;
;usage examples
;	* f=lineflx(rd_line(logP=20,wvl=w,logT=T,/allah),T,w,dem=10.D^(5*t))
;	* fl=rd_line(wr=16.776,wvl=w,logT=T,/allah) & t=10.^(t)
;	  f=lineflx(fl,T,w,dem=10.D^(5*alog10(t)),/temp)
;
;history
;	vinay kashyap (Nov96)
;	added keywords EFFAR and WVLAR (VK; Jan96)
;	added call to WHEE, removed call to KILROY,
;	  ensured DEM stretches to fit logT, added case of EM masquerading
;	  as DEM, removed the [/sr] from the output (VK; Feb97)
;	changed integration style; restructured to ignore effective area
;	  if not supplied; added keyword REGRID (VK; Apr97)
;	corrected ln v/s log10 bug (VK; Jun97)
;	added keyword NOPH; bug correction -- dlogT should in in base 10
;	  (VK; Jan98)
;	changed keyword KEV to IKEV (VK; DecMM)
;	EFFAR and WVLAR ignored if set to 0 (VK; JanMMI)
;	added keyword NHNE (VK; Jun02)
;	added keyword NOABUND (VK; Apr03)
;	changed bad input check from LINE(0)=-1 to total(LINE)=-N(LINE)
;	  (VK; Feb06)
;	changed back integration style to simple from trapezoidal (VK; Mar06)
;	changed check for whether WVL are given (VK; Jun07)
;	bug fix: was including 0.83 factor even when continuua and broad-band
;	  emissivities were being used [thx P.Testa]; now anything with Z=1
;	  will exclude the NHNE correction -- for H lines, fold the NHNE factors
;	  into the emissivities manually (VK; Dec10)
;-

;	check input
;is the cooling function supplied?
if n_elements(line) eq 0 then begin
  print,'Usage: fx=lineflx(line,logT,wvl,Z,DEM=DEM,/temp,abund=abund,$'
  print,'       /noph,nhne=nhne,effar=effar,wvlar=wvlar,/ikev,/noabund)'
  print,'  compute observed flux ([ph/s/cm^2] or [ph/s]) from spectral line'
  return,0.
endif
szc=size(line) & nt=szc(1) & nw=szc(2)
if szc(0) eq 0 and szc(1) ne 0 then begin
  ;LINE is input as a scalar -- convert to vector
  line=[line] & szc=size(line) & nt=szc(1) & nw=1L
endif
if szc(0) eq 1 then nw=1L

;was the cooling function properly read?
;if line(0) eq -1. then begin	-- this is needlessly restrictive
if total(line) eq -1L*n_elements(line) then begin	;(means filled with -1s
  message,'	No lines input; returning',/info & return,line
endif						;LINE==-1)

;are the temperatures given?
if n_elements(logT) ne nt then begin		;size mismatch
  message,'Cooling Function and Temperature arrays do not match',/info
  return,line
endif
tt=[logT]

;are the wavelengths given?
ang=[0.] & if keyword_set(wvl) then ang=[wvl]
if szc(0) ne 1 then begin
  if n_elements(wvl) ne nw then begin		;size mismatch
    message,'Cooling Function and wavelength arrays do not match',/info
    return,line
  endif
  ang=[wvl]
endif
ang=abs(ang)					;the CHIANTI case

;if wavelengths are given, so also should atomic numbers be
atom=intarr(nw)+1				;default is H
;if ang(0) ne 0 then begin			;wavelengths are given
if total(ang) ne 0 then begin			;wavelengths are given
  nz=n_elements(z)
  if nz gt 0 and nz le nw then atom(0)=z
  if nz gt nw then atom(*)=z(0:nw-1)
endif
oz=where(atom eq 1,moz)

;is the DEM given?
diffem=dblarr(nt)+1d14					;default
if keyword_set(dem) then begin
  ndem=n_elements(dem)
  if ndem eq 1 then diffem(*)=dem		;scalar
  if ndem eq nt then diffem=[dem]		;what the doc ordered
  if ndem gt 1 and ndem ne nt then begin	;interpolate
    message,'Interpolating DEM onto a default grid',/info
    t0=logT(0) & t1=logT(nt-1) & ti=findgen(ndem)*(t1-t0)/(ndem-1.)+t0
    tmp=alog10(dem>1e-20)
    diffem=interpol(tmp,ti,logT) & diffem=10.^(diffem)
    oo=where(diffem lt 1e-19) & if oo(0) ne -1 then diffem(oo)=0.
  endif
endif

;are the temperatures in T or logT?
if keyword_set(temp) then begin
  if nt gt 1 then diffem=diffem*tt		;[cm^-5/K]->[cm^-5/logK]
  ;(if nt eq 1, then DEM is actually EM)
  tt=alog10(abs(logT))				;convert to log10(T)
endif

;is the N(H)/N(e) fraction given?
nh_ne=fltarr(nt)+0.83 & nhe=n_elements(nhne)
if nhe gt 0 then begin
  if nhe eq 1 or nhe ne nt then begin
    if nhne(0) le 0 then nh_ne(*)=1. else nh_ne(*)=nhne(0)
  endif else nh_ne=float(nhne)>0.
endif
if nz eq 0 then nh_ne(*)=1.	;ignore if Z is not given
if moz ne 0 then nh_ne(oz)=1.	;or is "H" because this is how
				;continuum and solar calcs are done

;check if temperature gridding is uniform
if nt gt 1 then begin
  dlogt=tt(1:*)-tt & mdt=total(dlogt)/(nt-1.)
  sdt=sqrt(total((dlogt-mdt)^2)/(nt-1.))
  if sdt/mdt gt 0.1 and not keyword_set(regrid) then $
	message,'Temperature gridding is not regular, but ignoring!',/info
endif

;regrid if required
if keyword_set(regrid) then begin
  message,'	sorry, REGRID not implemented yet!',/info
endif

;convert dlog10 to dl
;if nt gt 1 then mdt=alog(10.^(tt(1)-tt(0)))	;THIS IS WRONG!!
if nt gt 1 then mdt=alog(10.)*mdt

;are the abundances given?
if n_elements(abund) ne 30 then abund=getabund('anders & grevesse')

;are the effective areas and wavelengths given?
nar=n_elements(effar) & nwv=n_elements(wvlar)
if nar eq 1 and not keyword_set(effar) then nar=0
if nwv eq 1 and not keyword_set(wvlar) then nwv=0
if nar eq 0 and nwv eq 0 then area=0*ang+1. else begin
  if nar ne 0 then areff=[effar] else areff=[1.,1.]
  if nwv ne 0 then arwvl=[wvlar] else arwvl=[min(abs(ang)),max(abs(ang))]
  if nar ne nwv then begin
    message,'Effective area v/s wavelength mismatch; doing nothing',/info
    area=0*ang+1.
  endif else begin
    ow=sort(arwvl) & arwvl=arwvl(ow) & areff=areff(ow)	;sort
    area=interpol(areff,arwvl,abs(ang))>0	;effective areas at WVLs
  endelse
endelse

;	integrate
flx=dblarr(nw) & avdem=diffem*1d-23
for iw=0L,nw-1L do begin
  if keyword_set(noabund) then zab=1. else $
	zab=(abund([atom(iw)-1]))(0)			;abundance
  if nt eq 1 then begin			;(EM, not DEM
    cfnc=[nh_ne(0)*line(0)]
    flx(iw)=zab*avdem*cfnc
	;[ergs/s/cm^2]=[1e23 cm^-5]*[1e-23 ergs cm^3/s]
  endif else begin			;)(DEM, not EM
    cfnc=reform(nh_ne(*)*line(*,iw))
    tmp=zab*avdem*cfnc
    ;tmp([0,nt-1])=0.5*tmp([0,nt-1])
    ;;	the 0.5 weighting at the ends is due to trapezoidal integration
    ;which has been commented out by VK because it is just very annoying
    ;to have those 2x dips at the temperature grid extremes (Mar2006)
    flx(iw)=mdt*total(tmp,/nan)
	;[ergs/s/cm^2]=[1e23 cm^-5/logK]*[logK]*[1e-23 ergs cm^3/s]
	;NOTE:	divide by (4.D*!dpi) to include [/sr] -- not done here
    ;;;;flx(iw)=int_tabulated(tt,zab*avdem*cfnc)
    ;;;;	uses 5-pt Newton-Cotes, but slows LINEFLX to a
    ;;;;	crawl even on a SPARC Ultra
  endelse				;NT.ne.1)
endfor

;	convert from ergs/s/cm^2 to ph/s
if not keyword_set(noph) then begin
  h=6.626176e-27 & c=2.9979e10		;[ergs s] & [cm s^-1]
  nrg=h*c*1e8/abs(ang)			;[ergs]
  if keyword_set(ikev) then nrg=1.6021892e-9*abs(ang)	;[ergs]
  flx=flx/nrg				;[ph/s/cm^2]
endif
flx=flx*area				;[(erg|ph)/s]

return,flx
end
