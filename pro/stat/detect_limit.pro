function detect_limit,bkg,nsig,asrc=asrc,abkg=abkg,bgalt=bgalt,abgalt=abgalt,$
	nsim=nsim,ulsig=ulsig,ulsim=ulsim,gaussy=gaussy,nxbin=nxbin,$
	verbose=verbose, _extra=e
;+
;function	detect_limit
;	compute and return the counts upper limit for
;	detection of a source at a given significance,
;	given the background counts.
;
;description
;	compute the cumulative significance of obtaining a specified
;	number of counts given the background, and assume that a source
;	would be considered detected if the counts were to exceed the
;	NSIG threshold.
;
;	see Pease, Drake, & Kashyap (2006, ApJ 636, 426) for full description.
;	briefly, computes the probability that as many as D counts can be
;	observed for a given background b, p(=<D|b), and thence the probability
;	that a given number of counts can be obtained in an observation simply
;	due to the background.  The number of counts required for a detection
;	at a specified probability is the upper limit.  Note also that this
;	goes only halfway towards a full upper limit (see Kashyap et al., 2008,
;	AAS-HEAD, 2008.03.03) in that this produces an upper limit in counts
;	space and not in intrinsic flux space.
;
;syntax
;	ul=detect_limit(bkg,nsig,asrc=asrc,abkg=abkg,bgalt=bgalt,abgalt=abgalt,$
;	nsim=nsim,ulsig=ulsig,ulsim=ulsim,/gaussy,nxbin=nxbin,verbose=verbose)
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
;	bgalt	[INPUT; default=0] alternative set of background
;		contamination, say from a different area of the instrument,
;		or from a model, or from an extended source, etc.
;		* may be an array
;		* NOTE: the _same_ BGALT is appended to _all_ elements of BKG
;	abgalt	[INPUT; default=asrc] area in which BGALT is "collected"
;		* size must match BGALT
;		* if size < size(BGALT), first element is assumed to be
;		  the default
;	nsim	[INPUT; default=0] number of Monte Carlo simulations to run
;		to account for error in background
;		* setting this results in computing the upper limit for
;		  a number of realizations of BKG; the resulting 1-sigma
;		  range in the value of computed upper limits is reported
;		  in ULSIG, and a conservative upper limit based on
;		  combining all of the simulations is returned in ULSIM[*,0]
;	ulsig	[OUTPUT] the 1-sigma error on the upper limit, estimated
;		by bootstrapping BKG
;	ulsim	[OUTPUT] a 2D array of size (NBKG,NSIM+1) which contains
;		all the simulated limits
;		* ULSIM[*,0] is identical to the primary output if NSIM=0
;		* for NSIM>0, ULSIM[*,0] is the conservative limit that
;		  is derived from the coadded probability distributions
;		  that take into account the variations in the background
;	gaussy	[INPUT] if set, computes the limit corresponding to the
;		significance matching the location of the NSIG-sigma
;		*intercept* of a Gaussian, rather than matching the total
;		area under the curve.
;	nxbin	[INPUT] number of bins to use in the integration
;		* the integration is carrid out over a range of
;		  0..5*E(bg) or 20, whichever is greater.
;		  by default, the number of bins is set by the step size,
;		  which is set to 1 count,
;		  -- unless E(bg) < 1, in which case a bin width of 0.05
;		  is used by default
;		* changing NXBIN does not change the range, but only
;		  changes the bin size.
;		* a hard lower limit of 20 is set -- cannot use a bin
;		  width larger than 1 count
;	verbose	[INPUT] controls chatter
;
;subroutines
;	LNPOISSON()
;	KILROY
;
;history
;	vinay kashyap (Apr2004)
;	added keyword ULSIG; changed name from PUPLIM to DETECT_LIMIT
;	  (VK; May2004)
;	added keyword NXBIN (VK; Sep2004)
;	added keywords BGALT and ABGALT (VK; Dec2004)
;	modified output behavior of ULSIM[*,0]; now BGALT can be 0
;	  (VK; Mar2005)
;	fixed bug for when NXBIN is set; made the numerical precision
;	  problem at high NSIG explicit by making the code return -1 as
;	  the UL (VK; Mar2009)
;-

message,'Obsolete.  Look at ULIM_SIMPLEPOI() instead.',/informational

;	usage
ok='ok' & np=n_params() & nb=n_elements(bkg) & ns=n_elements(nsig)
if np eq 0 then ok='Insufficient parameters' else $
 if nb eq 0 then ok='BKG is undefined'
if ok ne 'ok' then begin
  print,'Usage: ul=detect_limit(bkg,nsig,asrc=asrc,abkg=abkg,nsim=nsim,$'
  print,'       bgalt=bgalt,abgalt=abgalt,ulsig=ulsig,ulsim=ulsim,$'
  print,'       /gaussy,nxbin=nxbin,verbose=verbose)'
  print,'  return counts limit for detection in the presence of background'
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
nalt1=n_elements(bgalt) & nalt2=n_elements(abgalt)
if nalt1 gt 0 then begin
  altbg=[bgalt[*]]
  altareab=fltarr(nalt1)-1
  if nalt2 gt 0 then begin
    altareab[*]=abgalt[0]
    altareab[0L:(nalt1<nalt2)-1L]=abgalt[*]
  endif
endif
;
msim=intarr(nb) & ksim=n_elements(nsim)
if ksim gt 0 then begin
  msim[*]=nsim[ksim-1L]
  if ksim lt nb then msim[0L:ksim-1L]=nsim[*] else msim[*]=nsim[0L:nb-1L]
endif

;	output
ul=fltarr(nb) & ulsig=ul & ulsim=fltarr(nb,max(msim)+1L)

;	get upper limits
for i=0L,nb-1L do begin			;{for each given BKG
  numsim=msim[i] > 0
  if numsim gt 0 then tmpsim=fltarr(numsim)
  cgauss=errorf(msig[i]/sqrt(2.D))
  if keyword_set(gaussy) then begin
    cgauss=cgauss/2.D + 0.5D
    if vv gt 5 then print,$
	'Computing the limit corresponding to a significance of ',cgauss
  endif
  for j=0L,numsim do begin		;{for each realization
    if vv gt 5 then kilroy
    if areab[i] gt 0 then area_ratio=areas[i]/areab[i] else area_ratio=1.
    bg=bb[i]*area_ratio
    if nalt1 gt 0 then begin
      ob=where(altareab gt 0,mob) & alt_area_ratio=fltarr(nalt1)+1.
      if mob gt 0 then alt_area_ratio[ob]=areas[i]/altareab[ob]
      alt_bg=altbg*alt_area_ratio
      bg=bg+total(alt_bg)
    endif
    if j eq 0 then begin
      nx=0
      if keyword_set(nxbin) then begin	;(user defined number of bins
	if long(nxbin[0]) ge 20L then nx=long(nxbin[0])
      endif				;NXBIN)
      if nx eq 0 then begin 		;(default grid
	dx=1. & if bg lt 1 then dx=0.05
	nmax=long(5*bg)>20L & nx=long(nmax/dx)+1
	;if bg lt 1 then nx=400L else nx=long(5*bg)>20L
	x=findgen(nx+1L)*dx
      endif else x=findgen(nx+1L)	;default grid)
    endif
    if j gt 0 then begin
      bg=randomu(seed,poisson=bb[i])*area_ratio
      if nalt1 gt 0 then begin
	alt_bg_sim=fltarr(nalt1)
	for k=0L,nalt1-1L do $
	 if altbg[k] gt 0 then $
	  alt_bg_sim[k]=randomu(seed,poisson=altbg[k])*alt_area_ratio[k]
	bg=bg+total(alt_bg_sim)
      endif
    endif
    dpr=exp(lnpoisson(x,bg))
    if j eq 0 then mcdpr=dpr else mcdpr=mcdpr+dpr
    cpr=0.D*dpr+dpr[0] & for k=1L,nx do cpr[k]=cpr[k-1L]+dpr[k]
    ;cpr=0.D*dpr & for k=1L,nx do cpr[k]=cpr[k-1L]+dpr[k-1L]
    cpr=cpr/max(cpr)
    tmp=interpol(x,cpr,cgauss)	;Gaussian NSIG equivalent
    if finite(tmp) eq 0 then tmp=-1	;a temporary hack
    if j eq 0 then ul[i]=tmp[0] else tmpsim[j-1L]=tmp[0]
  endfor				;J=0,NUMSIM}
  mccpr=0.*mcdpr+mcdpr[0] & for k=1L,nx do mccpr[k]=mccpr[k-1L]+mcdpr[k]
  mccpr=mccpr/max(mccpr) & tmp=interpol(x,mccpr,cgauss)
  if finite(tmp) eq 0 then tmp=0.
  ulsim[i,0]=tmp[0]
  if numsim gt 0 then ulsim[i,1L:numsim]=tmpsim[*]
  if numsim gt 1 then ulsig[i]=stddev(tmpsim)
  if vv gt 10 then plot,x,dpr,xtitle='ct',ytitle='p(ct|bkg)'
  if vv gt 10 then oplot,tmp*[1,1],[1e-30,1]
  if vv gt 100 then stop,'HALTing; type .CON to continue'
endfor					;I=0,NB-1}

return,ul
end
