function evtflux,swvl,seffar,nbkg,bkgscal,nsim=nsim,$
	flxerr=flxerr,verbose=verbose,_extra=e
;+
;function	evtflux
;	given a set of events with known wavelengths and
;	effective areas, computes and returns the flux and
;	uncertainty on the flux
;
;syntax
;	flx=evtflux(swvl,seffar,nbkg,bkgscal,nsim=nsim,$
;	flxerr=flxerr,verbose=verbose)
;
;parameters
;	swvl	[INPUT; required] wavelength of photon in [Ang]
;	seffar	[INPUT; required] effective area at each photon
;		* size of SWVL and SEFFAR must match
;	nbkg	[INPUT; required] number of background events
;		* if not scalar, then number of elements in the
;		  supplied vector is taken to be NBKG
;	bkgscal	[INPUT; required] ratio of background to source area
;
;keywords
;	nsim	[INPUT; default=100] number of draws to make to
;		account for background in estimate of uncertainty
;	flxerr	[OUTPUT] error on flux, computed by combining Poisson
;		variance with variance due to drawing out background
;		events
;
;history
;	Vinay Kashyap (2014oct)
;-

;	usage
ok='ok' & np=n_params() & nph=n_elements(swvl) & nea=n_elements(seffar)
mbkg=n_elements(nbkg) & nscal=n_elements(bkgscal)
if np lt 4 then ok='Insufficient parameters' else $
 if nph eq 0 then ok='SWVL is not defined' else $
  if nea eq 0 then ok='SEFFAR is not defined' else $
   if nph ne nea then ok='SWVL and SEFFAR are incompatible' else $
    if mbkg eq 0 then ok='NBKG is not defined' else $
     if nscal eq 0 then ok='BKGSCAL is not defined'
if ok ne 'ok' then begin
  print,'Usage: flx=evtflux(swvl,seffar,nbkg,bkgscal,nsim=nsim,$
  print,'       flxerr=flxerr,verbose=verbose)
  print,'  compute and return flux and error on flux for given list of events'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
numbkg=mbkg & if mbkg eq 1 then numbkg=nbkg[0]
backscal=bkgscal[0]
;
numsim=100L & if keyword_set(nsim) then numsim=long(nsim[0])>1

;	convert ph to erg/cm^2
sevt = (1.9864776d-8/swvl)/seffar

;	background draws
simflx=dblarr(nsim)
ibkg=randomu(seed,numsim,poisson=numbkg*backscal)
for i=0L,numsim-1L do begin
  isrc=(nph-ibkg[i]) > 0
  if isrc eq 1 then begin
    ii=long(randomu(seed)*nph)
    simflx[i]=evt[ii]
  endif
  if isrc gt 1 then begin
    ;	this part does sampling with replacement
    ;ii=randomu(seed,isrc)*nph
    ;simflx[i]=mean(evt[ii],/double,/nan)
    ;	this part does sampling without replacement
    rr=randomu(seed,nph) & ii=(sort(rr))[0L:isrc-1L]
    simflx[i]=mean(evt[ii],/double,/nan)
  endif
endfor

;	this is the background subtracted flux
flx=mean(simflx,/double,/nan)

;	the uncertainty on the flux
;	first part is Poisson, second is due to background
flxerr=sqrt((flx/sqrt(nph))^2+stddev(simflx)^2)

return,flx
end
