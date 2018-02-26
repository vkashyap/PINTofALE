function discriminatemp,emis,ferr,EM0=EM0,outflx=outflx,$
	wvl=wvl,Z=Z,logT=logT,abund=abund,normark=normark,$
	outeflx=outeflx, verbose=verbose, _extra=e
;+
;function	discriminatemp
;	an implementation of Mark Weber's scheme for determining
;	the temperature discrimination ability of a given set of
;	emissivity functions.  the output is a square matrix
;	of the same size as the temperature grid, showing the value
;	of a statistic that can be used to infer whether two
;	corresponding temperatures can be distinguished.
;
;	the technique is simple -- compute fluxes at two temperatures
;	for the same emission measure, and for an assumed uncertainty,
;	ask whether the two sets of fluxes differ in a statistically
;	significant manner.  if they do, the two temperatures are
;	distinguishable, and not if they do not.
;
;syntax
;	dtemp=disriminatemp(emis,ferr,EM0=EM0,outflx=outflx,$
;	wvl=wvl,Z=Z,logT=logT,abund=abund,normark=normark,$
;	outeflx=outeflx,verbose=verbose,$
;	/noph,effar=effar,wvlar=wvlar,/ikev,metals=metals,fipbias=fipbias)
;
;parameters
;	emis	[INPUT; required] array of emissivities
;		* EMIS is a 2D array with size (LOGT, WVL)
;		* the units are assumed to be [1e-23 ergs cm^3/s]
;		  for anything else (e.g., for solar work), put the
;		  difference into EM0 and be sure to set the keyword
;		  NOPH=1, because it gets pushed into LINEFLX()
;	ferr	[INPUT] the assumed uncertainties on the lines
;		* if not given, taken to be 0.1 of the computed
;		  fluxes for all lines
;		* if scalar, then
;		    if >0 and <1, taken to be that fraction of the
;		    computed fluxes for all lines
;		    if >1 and <100, taken to be that percentage of
;		    the computed fluxes for all lines
;		    if >100, the reciprocal is taken to be the
;		    fraction of the computed fluxes for all lines
;		    if <0, the abs value is taken to be a constant
;		    absolute uncertainty for all lines
;		* if vector, each individual value is used as given,
;		  with the same condition on each as for the scalar
;		  case
;		* e.g., if there are 4 lines, and FERR=[-10,0.3,5,0],
;		  after the 4 fluxes are computed, the corresponding
;		  errors are set to [10.,0.3*F[1],0.05*F[2],0.1*F[3]]
;		* the best way to set reasonable errors is to run
;		  this program once, extract the computed fluxes
;		  using keyword OUTFLX, figure out the appropriate
;		  errors, and then feed it back into another run
;
;keywords
;	EM0	[INPUT] the emission measure to use
;		* default is 1d14 cm^-5, fwiw
;	outflx	[OUTPUT] the fluxes calculated for each line for each
;		temperature, is an array of size (LOGT,WVL), same as EMIS
;	outeflx	[OUTPUT] the flux errors calculated for each line for each
;		temperature, is an array of size (LOGT,WVL), same as EMIS
;	wvl	[INPUT] line wavelengths for which EMIS is given
;		* used only if keyword NOPH is _not_ set (i.e.,
;		  in the conversion of [erg/...] to [ph/...])
;		  and if size matches the 2nd dimension of EMIS
;	Z	[INPUT] atomic numbers of elements contributing
;		to EMIS
;		* used only if size matches 2nd dimension of EMIS,
;		  and is used to multiply EMIS with the appropriate
;		  abundance
;		* if not given, assumed to be 1 (i.e., H), which
;		  is essentially a way to ignore ABUND
;	logT	[INPUT] log_10(Temperature [K]) at which EMIS are given
;		* unused
;	abund	[I/O] element abundances
;		* if size smaller than 30, calls GETABUND() and resets
;		  to use Anders & Grevesse
;	normark	[INPUT] a normalization factor, usually set to the
;		maximum number of filters or lines that can be used,
;		and whose squared reciprocal divided into the output --
;		Mark Weber uses this to "normalize" the results across
;		different filter choices
;		* hard minimum value is 1, also the default
;		* if this is set, the diagonal elements of the output,
;		  which are identically zero, are reset to the mean of
;		  the nearest neighbors
;		* if vector, only the first element is used and all the
;		  extra elements are ignored
;	verbose	[INPUT] controls chatter
;
;	_extra	[INPUT ONLY] pass defined variables to subroutines
;		LINEFLX : /NOPH, /IKEV, EFFAR, WVLAR
;		GETABUND : METALS, FIPBIAS
;
;description
;	Weber et al., 2008, SPD?
;
;history
;	vinay kashyap (2009mar)
;	bigfix: wasn't dealing with vector FERR gracefully (VK; 2010aug)
;-

;	usage
ok='ok' & np=n_params() & nem=n_elements(emis) & sze=size(emis)
if np eq 0 then ok='Insufficient parameters' else $
 if nem eq 0 then ok='EMIS is undefined' else $
  if sze[0] gt 3 then ok='EMIS cannot be more than 2D'
if ok ne 'ok' then begin
  print,'Usage: dtemp=discriminatemp(emis,ferr,EM0=EM0,outflx=outflx,$'
  print,'       wvl=wvl,Z=Z,logT=logT,abund=abund,normark=normark,$'
  print,'       outeflx=outeflx,verbose=verbose,$'
  print,'       /noph,effar=effar,wvlar=wvlar,/ikev,metals=metals,$'
  print,'       fipbias=fipbias)
  print,'  determines the temperature discrimination ability of a chosen'
  print,'  set of emissivities'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
nT=sze[1] & nL=sze[2]
;
nerr=n_elements(ferr) & fsig=fltarr(nL)+0.1
if nerr gt 0 then begin
  if nerr eq 1 then fsig[*]=ferr[0] else $
    fsig[0L:(nerr<nL)-1L]=ferr[0L:(nerr<nL)-1L]
endif
for k=0L,nL-1L do begin
  ff=fsig[k]
  ;if ff gt 0 and ff lt 1 then fsig[k]=fsig[k]
  if ff ge 1 and ff lt 100 then fsig[k]=fsig[k]/100.
  if ff ge 100 then fsig[k]=1./fsig[k]
endfor
;
EMdef=1d14 & if keyword_set(EM0) then EMdef=EM0[0]
;
wave=fltarr(nL)+1. & if n_elements(wvl) eq nL then wave=wvl[*]
;
ZZ=intarr(nL)+1 & if n_elements(Z) eq nL then ZZ=Z[*]
;
if n_elements(abund) lt 30 then abund=getabund('anders & grevesse', _extra=e)
;
marknorm=1.0 & if n_elements(normark) gt 0 then marknorm=abs(normark[0])>1.

;	output
dtemp=bytarr(nT,nT)
outflx=dblarr(nT,nL) & outeflx=outflx
dTmatrix=dblarr(nT,nT)

;	shtep thru the temperatures and compute the fluxes
for i=0L,nT-1L do begin
  fx=fltarr(nL)
  for k=0L,nL-1L do fx[k]=lineflx(emis[i,k],i,wave[k],ZZ[k],DEM=EMdef,abund=abund, _extra=e)
  efx=abs(fsig) & for k=0L,nL-1L do if fsig[k] gt 0 then efx[k]=fx[k]*fsig[k]
  outflx[i,*]=fx & outeflx[i,*]=efx
endfor

;	construct the discriminability matrix
for i=0L,nT-1L do for j=0L,nT-1L do $
  dTmatrix[i,j]=total((outflx[i,*]-outflx[j,*])^2/(outeflx[i,*]^2+outeflx[j,*]^2),/nan)

;	renormalize
dTmatrix=dTmatrix/marknorm^2
if keyword_set(normark) then begin
  for i=0L,nT-1L do begin
    arr=0
    if i gt 0L then begin
      if keyword_set(arr) then arr=[arr,dTmatrix[i,i-1],dTmatrix[i-1,i]] else arr=[dTmatrix[i,i-1],dTmatrix[i-1,i]]
    endif
    if i lt nT-1L then begin
      if keyword_set(arr) then arr=[arr,dTmatrix[i,i+1],dTmatrix[i+1,i]] else arr=[dTmatrix[i,i+1],dTmatrix[i+1,i]]
    endif
    dTmatrix[i,i]=mean(arr,/nan)
  endfor
endif

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,dTmatrix
end
