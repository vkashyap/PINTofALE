function puplim,bkg,nsig,asrc=asrc,abkg=abkg,nsim=nsim,ulsim=ulsim,$
	gaussy=gaussy,verbose=verbose, _extra=e
;+
;function	puplim
;	compute and return the counts upper limit for
;	detection of a source at a given significance,
;	given the background counts.
;
;syntax
;	ul=puplim(bkg,nsig,asrc=asrc,abkg=abkg,nsim=nsim,ulsim=ulsim,$
;	/gaussy,verbose=verbose)
;
;parameters
;	bkg	[INPUT; required] counts in the background
;	nsig	[INPUT] the Gaussian-equivalent sigma at which
;		to compute upper limit to detection
;		* if not given, assumed to be 3
;
;keywords
;	asrc	[INPUT; default=1] area in which source counts are collected
;	abkg	[INPUT; default=1] area in which background counts are collected
;	nsim	[INPUT; default=0] number of Monte Carlo simulations to run
;		to account for error in background
;		* setting this results in computing the upper limit for
;		  a number of realizations of BKG; the resulting 1-sigma
;		  range in the value of computed upper limits is reported
;		  in ULSIM
;	ulsim	[OUTPUT] the 1-sigma error on the upper limit, estimated
;		by bootstrapping BKG
;	gaussy	[INPUT] if set, computes the limit corresponding to the
;		significance matching the location of the NSIG-sigma
;		*intercept* of a Gaussian, rather than matching the total
;		area under the curve.
;	verbose	[INPUT] controls chatter
;
;description
;	compute the cumulative significance of obtaining a specified
;	number of counts given the background, and assume that a source
;	would be considered detected if the counts were to exceed the
;	NSIG threshold.
;
;subroutines
;	LNPOISSON()
;	KILROY
;
;history
;	vinay kashyap (Apr2004)
;-

message,'OBSOLETE!  Use DETECT_LIMIT() instead',/informational

;	usage
ok='ok' & np=n_params() & nb=n_elements(bkg) & ns=n_elements(nsig)
if np eq 0 then ok='Insufficient parameters' else $
 if nb eq 0 then ok='BKG is undefined'
if ok ne 'ok' then begin
  print,'Usage: ul=puplim(bkg,nsig,asrc=asrc,abkg=abkg,nsim=nsim,ulsim=ulsim,$'
  print,'       /gaussy,verbose=verbose)'
  print,'  return counts upper limit for detection'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
bb=[bkg[*]]
;
msig=fltarr(nb)+3.
if ns gt 0 then begin
  msig[*]=nsig[ns-1L]
  if ns lt nb then msig[0L:ns-1L]=nsig[*] else msig[*]=nsig[0L:nb-1L]
endif
;
areas=fltarr(nb)+1. & nas=n_elements(asrc)
if nas gt 0 then begin
  areas[*]=asrc[nas-1L]
  if nas lt nb then areas[0L:nas-1L]=asrc[*] else areas[*]=asrc[0L:nb-1L]
endif
;
areab=fltarr(nb)+1. & nab=n_elements(abkg)
if nab gt 0 then begin
  areab[*]=abkg[nab-1L]
  if nab lt nb then areab[0L:nab-1L]=abkg[*] else areab[*]=abkg[0L:nb-1L]
endif
;
msim=fltarr(nb)+3. & ksim=n_elements(nsim)
if ksim gt 0 then begin
  msim[*]=nsim[ksim-1L]
  if ksim lt nb then msim[0L:ksim-1L]=nsim[*] else msim[*]=nsim[0L:nb-1L]
endif

;	output
ul=fltarr(nb) & ulsim=ul

;	get upper limits
for i=0L,nb-1L do begin			;{for each given BKG
  numsim=msim[i] > 0
  if numsim gt 0 then tmpsim=fltarr(numsim)
  for j=0L,numsim do begin
    if vv gt 5 then kilroy
    if areab[i] gt 0 then area_ratio=areas[i]/areab[i] else area_ratio=1.
    bg=bb[i]*area_ratio
    if j gt 0 then bg=randomu(seed,poisson=bb[i])*area_ratio
    nx=long(5*bg)>20L & x=findgen(nx+1L)
    dpr=exp(lnpoisson(x,bg))
    cpr=0.D*dpr+dpr[0] & for k=1L,nx do cpr[k]=cpr[k-1L]+dpr[k]
    ;cpr=0.D*dpr & for k=1L,nx do cpr[k]=cpr[k-1L]+dpr[k-1L]
    cpr=cpr/max(cpr)
    cgauss=errorf(msig[i]/sqrt(2.D))
    if keyword_set(gaussy) then begin
      cgauss=cgauss/2.D + 0.5D
      if vv gt 5 then print,$
	'Computing the limit corresponding to a significance of ',cgauss
    endif
    tmp=interpol(x,cpr,cgauss)	;Gaussian NSIG equivalent
    if j eq 0 then ul[i]=tmp[0] else tmpsim[j-1L]=tmp[0]
  endfor
  if numsim gt 0 then ulsim[i]=stddev(tmpsim)
  if vv gt 10 then plot,x,dpr,xtitle='ct',ytitle='p(ct|bkg)'
  if vv gt 10 then oplot,tmp*[1,1],[1e-30,1]
  if vv gt 100 then stop,'HALTing; type .CON to continue'
endfor					;I=0,NB-1}

return,ul
end
