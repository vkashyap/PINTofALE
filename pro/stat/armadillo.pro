function armadillo,src,bkg,quantil=quantil,quaup=quaup,quadn=quadn,$
	quasig=quasig,clev=clev,asrc=asrc,abkg=abkg,nsim=nsim,verbose=verbose,$
	simsrc=simsrc,simcdf=simcdf,simxx=simxx, _extra=e
;+
;function	armadillo
;	returns the quantiles of a background contaminated histogram
;
;	they say that there is nothin in the middle of the road except
;	yellow stripes and dead armadillos.  here however are medians
;	with error bars.
;
;syntax
;	x=armadillo(src,bkg,quantil=quantil,quaup=quaup,quadn=quadn,$
;	quasig=quasig,clev=clev,asrc=asrc,abkg=abkg,nsim=nsim,$
;	verbose=verbose,simsrc=simsrc,simcdf=simcdf,simxx=simxx,$
;	funit=funit,agamma=agamma,bgamma=bgamma,priorbg=priorbg,$
;	nsgrid=nsgrid,srcmin=srcmin,srcmax=srcmax)
;
;parameters
;	src	[INPUT; required] source spectrum
;	bkg	[INPUT] background spectrum
;		* size must match SRC.  if it does not,
;		-- if 0-element, assumed to be 0
;		-- if 1-element, assumed to be constant
;		-- if >1 and .NE. N(SRC), returns with error
;
;keywords
;	quantil	[INPUT] the quantile(s) to compute
;		* default is the median, quantile=0.5
;		* may be an array
;		* must be in the range 0..1. if not,
;		-- if >1 and <100, taken to be percentage number
;		-- if >100, taken to be (1 - reciprocal fraction)
;		-- if <0, abs value is taken to be reciprocal fraction
;	quaup	[OUTPUT] the upper bound on the quantile(s) at
;		confidence level CLEV
;	quadn	[OUTPUT] the lower bound on the quantile(s) at
;		confidence level CLEV
;	quasig	[OUTPUT] the std.dev. on the quantile(s)
;	clev	[INPUT] confidence level at which to get QUAUP and QUADN
;	asrc	[INPUT] area of the source region
;		* default=1
;	abkg	[INPUT] area of the background region
;		* default=1
;	nsim	[INPUT] number of simulations to run
;	verbose	[INPUT] controls chatter
;	simsrc	[OUTPUT] stores the counts from all the simulations,
;		returned in an array of size (NSIM,N(SRC))
;	simcdf	[OUTPUT] stores the cdfs from all the simulations,
;		returned in an array of size (NSIM,N(SRC))
;	simxx	[OUTPUT] stores the quantiles from all the simulations,
;		returned in an array of size (NSIM,N(QUANTIL))
;	_extra	[INPUT ONLY] pass defined keywords to
;		PPD_SRC: FUNIT,AGAMMA,BGAMMA,PRIORBG,NSGRID,SRCMIN,SRCMAX
;
;description
;	given counts histograms from a source region and an estimate
;	of the background contamination, computes the posterior probability
;	distribution at each bin, generates sample counts from the ppd,
;	and uses these monte carlo samples to estimate an average value
;	of the quantiles and errors on them.
;
;subroutines
;	PPD_SRC
;	PROB_GAMMADIST
;
;history
;	vinay kashyap (Mar2005)
;-

;	usage
ok='ok' & np=n_params() & ns=n_elements(src) & nb=n_elements(bkg)
if np eq 0 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='SRC undefined' else $
  if ns lt 3 then ok='SRC too small an array'
if ok ne 'ok' then begin
  print,'Usage: x=armadillo(src,bkg,quantil=quantil,quaup=quaup,quadn=quadn,$
  print,'       quasig=quasig,clev=clev,asrc=asrc,abkg=abkg,nsim=nsim,$
  print,'       verbose=verbose,simsrc=simsrc,simcdf=simcdf,simxx=simxx,$
  print,'       funit=funit,agamma=agamma,bgamma=bgamma,priorbg=priorbg,$
  print,'       nsgrid=nsgrid,srcmin=srcmin,srcmax=srcmax)'
  print,'  returns the quantiles of a background contaminated histogram'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
bg=0.*src
if nb gt 0 then begin
  if nb eq 1 then bg[*]=abs(bkg[0]) else begin
    if ns ne nb then begin
      message,'SRC and BKG are incompatible; returning',/informational
      return,-1L
    endif else bg=bkg
  endelse
endif
;
qq=0.5 & nq=n_elements(quantil)
if nq gt 0 then begin
  qq=fltarr(nq)+0.5
  for i=0L,nq-1L do begin
    tmp=quantil[i]
    if tmp lt 0 then qq[i]=abs(1./tmp) else $
     if tmp gt 100 then qq[i]=1.-(1./tmp) else $
      if tmp gt 1 and tmp lt 100 then qq[i]=tmp/100. else $
       qq[i]=tmp
    if vv ge 5 then message,'computing quantile '+strtrim(qq[i],2),$
	/informational
  endfor
endif
nq=n_elements(qq)
;
cl=0.68
if keyword_set(clev) then cl=clev[0]
if cl lt 0 then cl=abs(cl)
if cl gt 1 and cl lt 100 then cl=cl/100.
if cl gt 100 then begin
  message,'CLEV='+strtrim(clev[0],2)+' not understood; using 0.68',$
	/informational
  cl=0.68
endif
;
srcar=1. & if keyword_set(asrc) then srcar=asrc[0]
if srcar lt 0 then begin
  message,'ASRC cannot be negative; resetting to 1',/informational
  srcar=1.
endif
;
bkgar=1. & if keyword_set(abkg) then bkgar=abkg[0]
if bkgar lt 0 then begin
  message,'ABKG cannot be negative; resetting to 1',/informational
  bkgar=1.
endif
;
numsim=100L & if keyword_set(nsim) then numsim=long(nsim[0])>1L

;	outputs
xx=fltarr(nq) & quaup=xx & quadn=xx & quasig=xx
simcdf=fltarr(numsim,ns)
simsrc=fltarr(numsim,ns)
simxx=fltarr(numsim,nq)

;	now, for each bin, compute the posterior probability
;	distribution for the source intensity, and draw NSIM
;	samples from it
for i=0L,ns-1L do begin		;{step through SRC array

  ppd=ppd_src(src[i],bg[i],sgrid,asrc=srcar,abkg=bkgar,verbose=vv, _extra=e)
  ng=n_elements(sgrid)
  cdf=fltarr(ng)+ppd[0] & for j=1L,ng-1L do cdf[j]=cdf[j-1L]+ppd[j]
  cdf=cdf/cdf[ng-1L]

  r=randomu(seed,numsim)
  rg=interpol(sgrid,cdf,r)
  simsrc[*,i]=rg

endfor				;I=0,NS-1}

;	construct cdfs for each simulations and figure out the quantiles
simxx=fltarr(numsim,nq)
for i=0L,numsim-1L do begin
  tmp=simsrc[i,*]
  cdf=fltarr(ns)+tmp[0] & for j=1L,ns-1L do cdf[j]=cdf[j-1L]+tmp[j]
  cdf=cdf/max(cdf)
  simcdf[i,*]=cdf
  simxx[i,*]=interpol(lindgen(ns),cdf,qq)
endfor

;	get mean and errors
for i=0L,nq-1L do begin
  tmp=simxx[*,i]
  xx[i]=mean(tmp)
  quasig[i]=stddev(tmp)
  tmp=tmp[sort(tmp)] & cdf=findgen(numsim)/(numsim-1.)
  tmp2=interpol(tmp,cdf,0.5+[-cl/2.,cl/2.])
  quadn[i]=tmp2[0] & quaup[i]=tmp2[1]
endfor

return,xx
end
