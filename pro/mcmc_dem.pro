function mcmc_dem,wvl,flx,emis,Z=Z,logt=logt,diffem=diffem,abund=abund,$
	onlyrat=onlyrat,sysdev=sysdev,sysgrp=sysgrp,$
	fsigma=fsigma,ulim=ulim,demrng=demrng,abrng=abrng,$
	xdem=xdem,xab=xab,xnsig=xnsig,$
	nsim=nsim,nbatch=nbatch,nburn=nburn,nosrch=nosrch,savfil=savfil,$
	bound=bound,hwhm=hwhm,demerr=demerr,aberr=aberr,smoot=smoot,$
	smooscl=smooscl,$
	loopy=loopy,spliny=spliny,sampenv=sampenv,$
	storpar=storpar,storidx=storidx,stordem=stordem,bestpar=bestpar,$
	simprb=simprb,simdem=simdem,simabn=simabn,simflx=simflx,simprd=simprd,$
	mixsstr=mixsstr,mixlstr=mixlstr,verbose=verbose, _extra=e
;+
;function	mcmc_dem
;	performs Markov-Chain Monte-Carlo using a Metropolis algorithm
;	on a set of supplied line fluxes and returns an estimate of the
;	DEM that generates observed fluxes.
;
;syntax
;	dem=mcmc_dem(wavelengths,fluxes,emissivities,Z=Z,logt=logt,$
;	diffem=initial_DEM,abund=abundances,fsigma=flux_errors,$
;	onlyrat=onlyrat,sysdev=sysdev,sysgrp=sysgrp,$
;	sampenv=sampenv,storpar=storpar,storidx=storidx,$
;	stordem=stordem,bestpar=bestpar,$
;	simprb=simprb,simdem=simdem,simabn=simabn,simflx=simflx,simprd=simprd,$
;	ulim=upper_limit_flag,softlim=softlim,demrng=DEMrange,abrng=abrng,$
;	smoot=smooth,loopy=loopy,/spliny,verbose=verbose, xdem=DEM_bins,xab=ABUND_bins,$
;	nsim=nsim,nbatch=nbatch,/steps,weight=weight,/gauss,type=type,/boxcar,$
;	nburn=nburn,/nosrch,savfil=savfil,bound=bound,demerr=demerr,$
;	aberr=aberr,/temp,/noph,effar=effar,wvlar=wvlar,/ikev,/regrid,$
;	mproxy=mproxy,mplet=mplet, epsroot=epsroot,thrfrac=thrfrac,thrnum=thrnum,$
;	/dvdt,cbox=cbox,cwhisk=cwhisk,cdots=cdots,cbest=cbest,cmode=cmode,cgrey=cgrey,cnum=cnum,$
;	charthick=charthick,charsize=charsize,xrange=xrange,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title)
;
;parameters
;	wvl	[INPUT; required] wavelengths [Ang] at which fluxes are
;		observed
;		* NOTE: these are used only if the keyword NOPH is _not_ set,
;		  i.e., when the input is to be converted from [erg/...] to
;		  [ph/...].  In cases when the EMIS is precomputed with the
;		  effective area and abundances (e.g., for filter data from
;		  solar observatories), this is a dummy variable.
;		. in other words, if the input is in ergs/..., and the
;		  emissivities are also in ergs/..., then set /NOPH
;		. if input is in ph/..., and emissivities are in ergs/...,
;		  do not set NOPH
;		. when NOPH is set, then the conversion is not done within
;		  LINEFLX, so WVL does not get used
;	flx	[INPUT; required] flux [ph/s] at each WVL ([ph/s/cm^2]
;		if EFFAR and WVLAR are not passed to LINEFLX)
;	emis	[INPUT; required] line emissivities [1e-23 ergs cm^3/s]
;		* EMIS==EMIS(LOGT,WVL)
;		* DEM is computed at all temperatures that EMIS are defined.
;		* NOTE: the units are assumed to be [1e-23 ergs cm^3/s],
;		  but of course the code doesn't care what they are exactly.
;		  So it is perfectly possible to have them in units of, e.g.,
;		  [DN cm^5 s^-1 pix^-1] if one so chooses (see SOLAR_TRESP.PRO)
;
;keywords
;	Z	[INPUT] atomic numbers of the elements generating the
;		lines at WVL.
;		* if not given, assumed to be Fe
;		* if EMIS are precalculated including abundances, then
;		  set this to *1* (to mimic H)
;	logt	[INPUT] log_10(Temperature [K]) at which EMIS are given
;		* if not supplied, assumed to go from 4 to 8
;		* if there is a size mismatch, logT will get stretched/shrunk
;		  to cover EMIS
;	diffem	[INPUT] initial guess for DEM(LOGT); default is 1e14 [cm^-5]
;	abund	[I/O] abundances.  default is to use Anders & Grevesse (see
;		GETABUND)
;		* if ABRNG is non-trivial, will be allowed to vary and on
;		  exit will contain the MAP estimates of the abundances.
;		* minimum is hardcoded at 1e-20
;	fsigma	[INPUT] error on FLX; if not given, taken to be
;		1+sqrt(abs(FLX)+0.75)
;	ulim	[INPUT] long-integer array specifying which of the fluxes
;		are upper limits (1: UL, 0: not)
;	onlyrat	[INPUT] string array describing which of the input fluxes
;		must be considered only as ratios
;		* each element of ONLYRAT corresponds to the fate of the
;		  corresponding element of FLX
;		* basic format is: "sP#[,sP#[,...]]" where
;		  -- "#" is an integer flag describing the ratio being
;		     constructed
;		  -- "P" is a positional descriptor and can take on
;		     values N (for numerator) or D (for denominator)
;		  -- "s" stands for the sign with which the flux
;		     goes into the ratio "+" or "-"
;		* e.g., to construct a simple ratio FLUXES[2]/FLUXES[1],
;		  ONLYRAT=['','+D1','+N1']
;		* e.g., to construct two hardness ratios
;			FLUXES[3]/FLUXES[1] and
;			(FLUXES[3]-FLUXES[1])/(FLUXES[3]+FLUXES[1]),
;		  ONLYRAT=['','+D1,-N2,+D2','','+N1,+N2,+D2']
;	sysdev	[INPUT] systematic deviations in the underlying atomic
;		physics data, multiplies the emissivity for a given line
;		by a random number in each iteration to mimic systematic
;		uncertainty in EMIS
;		* it is assumed that these modifications are small -- if
;		  you do use this option, please check the output traces
;		  to make sure that the chain is stable
;		* should be an array of the same size as FLX, and each
;		  number defines the range over which EMIS[*,WVL] varies,
;		  10.^(randomn(seed)*SYSDEV[i])*EMIS[*,i]
;		* if size is smaller than FLX, the missing elements are
;		  set to 0, i.e., no modifications are made
;	sysgrp	[INPUT] group different rows in EMIS to vary the same way
;		in each iteration, is simply an array of 1-based grouping
;		indices
;		* e.g., if N(WVL)=6, and SYSGR=[0,1,2,1,0,2], then rows
;		  1 and 3 are modified in conjunction (i.e., if #1 is
;		  changed by a factor 10^(r*SYSDEV[1]), then #3 is changed
;		  by a factor 10^(r*SYSDEV[3]), where r is the same for both),
;		  2 and 5 are modified in conjunction, but 0 and 4 are
;		  _not_ tied together and are independently varied
;		* if SYSDEV is not defined, SYSGRP is ignored
;		* size must match that of FLX, but if smaller, the missing
;		  elements are set to -1, i.e., are all tied together
;		  - so if SYSGRP is not defined at all, all the rows are
;		    tied together, and only SYSDEV[0] is used
;	demrng	[INPUT] array containing the allowed range of variation
;		in DEM: DEMRNG(i,0) < DEM(i) < DEMRNG(i,1)
; 		* if not supplied, DEMRNG(*,0)=DEM/1e5 & DEMRNG(*,1)=DEM*1e5
;		* if DEMRNG(k,0)=DEMRNG(k,1), this component is frozen
;	abrng	[INPUT] as DEMRNG, for ABUND.
;		* if not supplied, ABRNG(*,0)=ABUND(*)=ABRNG(*,1)
;	xdem	[INPUT] bin boundaries for accumulating the realized DEMs
;		* if scalar or 1-element vector, taken to be bin-width in log10
;		* if not specified, runs from min(DEMRNG) to max(DEMRNG)
;		  in steps of 0.1 in log10
;	xab	[INPUT] bin boundaries for accumulating the realized ABUNDs
;		* if scalar or 1-element vector, taken to be bin-width
;		* if not specified, runs from min(ABRNG) to max(ABRNG) in
;		  steps of 0.1 in log10
;	xnsig	[INPUT; default=3] force a minimum search range of
;			XNSIG*sigma_determined_internally
;	nsim	[INPUT; default: 100; min=10] # of batches
;	nbatch	[INPUT; default: 10; min=5] # simulations per batch
;	nburn	[INPUT; default=NSIM/10] max # of simulations to run
;		before assuming that solution has stabilized
;		* note that this is a maximum -- if burn-in is achieved
;		  earlier, the algorithm will move on to standard MCMC
;	nosrch	[INPUT] if set, does *not* carry out a (crude) grid search
;		for a best guess initial starting point.
;		* if set to a IDL savefile name, restores the best guess
;		  from said file
;	savfil	[INPUT] if set, saves all variables to an IDL save file
;		prior to exit
;		* if not a string, or is an array, saves to 'mcmc.sav'
;	bound	[INPUT; default=0.9] confidence interval of interest
;		* minimum=0.1
;	hwhm	[INPUT] if set, returns the upper and lower Half-Widths
;		at Half-Maxmimum in DEMERR and ABERR.
;	demerr	[OUTPUT] confidence bounds on MAP estimates of DEM
;		* DEMERR(*,0) are lower bounds, DEMERR(*,1) are upper bounds
;	aberr	[OUTPUT] confidence bounds on MAP estimates of ABUND
;		* ABERR(*,0) are lower bounds, ABERR(*,1) are upper bounds
;	smoot	[INPUT] if set, multiplies the local length scale determined
;		using FINDSCALE by SMOOT.
;	smooscl	[INPUT] if set, bypasses FINDSCALE and uses the scales input
;		via this keyword as the local length scales
;		* _must_ be in units of [bins], not delta_T, dlogT, or whatever
;		* if vector and size matches LOGT, then uses the input values
;		  as the [bins] over which to smooth
;		* otherwise assumes constant smoothing scale corresponding to
;		  first element
;	loopy	[INPUT] if set, "smooths" by converting to T^3/2, loop-like DEM
;		* use keyword SLOOP to change the index
;		* set to a number >1 to limit the maximum number of
;		  temperature components
;	spliny	[INPUT] if set, "smooths" by making a spline interpolation
;	sampenv	[INPUT] envelope from which to sample the temperature bins.
;		this is a way to avoid T bins which may have no data and thus
;		work more efficiently to explore just that part of the parameter
;		space that is interesting.
;		* if not given, or is not understandable, assumed to be the
;		  simple envelope of all the emissivities in EMIS
;		* if set to:
;		  'ONLYRAT' - then the envelope is computed taking
;		  	into account the ratios
;		  'INCABUND' - then use emissivities weighted by the abundance
;		  'FLAT' - then make no distinction between different T bins
;		  'INCALLT' - then force a look at all T bins at a 1% level
;		  'STEPPED' - then zero everywhere EMIS does not cover the
;		  	T range, and at least 25% of maximum where it does
;		  'COVER' - then zero everywhere EMIS does not cover the
;		  	T range, and constant where it does
;		* NOTE: multiple SAMPENV settings may be specified, and
;		  they will be applied in the above sequence -- whether the
;		  sequence makes sense or not is up to the user to determine
;		* if set to a numerical array the same size as LOGT, then
;		  simply use that array as the envelope
;		* NOTE: recommend using 'COVER' or 'STEPPED' for EMIS
;		  derived from broad-band filters such as for solar data
;	storpar	[OUTPUT] parameters from the MCMC chain
;	storidx	[OUTPUT] index of the parameters stored in STORPAR
;	stordem	[OUTPUT] DEMs calculated from the parameters for given STORIDX
;	bestpar	[OUTPUT] best values of the parameters as found in the chain
;	simprb	[OUTPUT] the metric at the end of each batch, an array of size
;		(NSIM+1), with last element containing the value for the best-fit
;	simdem	[OUTPUT] the DEM at the end of each batch, an array of size
;		(NT,NSIM+1), with last column containing the best-fit
;	simabn	[OUTPUT] same as SIMDEM, but for abundances
;	simflx	[OUTPUT] same as SIMDEM, but for the fluxes with each [DEM,ABUND]
;	simprd	[OUTPUT] same as SIMFLX, but for the ratios, only if ONLYRAT is used
;       mixsstr [I/O]  Set to turn on line deblending via MIXIE().  
;                Line (and/or conitnuum) structure of
;                RD_LIST() format containing info (wvl,emis,z) on 
;                lines (and/or continuum) for which possible contaminants will be identified 
;                and for which correction factors  will be
;                calculated. Assumes that ion equilibria are
;                included. (see RD_LIST() keyword INCIEQ). 
;       mixlstr [I/O] Line (and/or conitnuum) structure of
;                RD_LIST() format containing info (wvl,emis,z) on 
;                lines (and/or continuum) for which possible contaminants will be identified 
;                and for which correction factors  will be
;                calculated. Assumes that ion equilibria are
;                included. (see RD_LIST() keyword INCIEQ)
;                set to turn on line deblending via MIXIE().  
;	verbose	[INPUT] explicitly set to 0 to avoid making plots
;		* THINK TWICE before doing this!!
;	_extra	[INPUT] pass defined keywords to subroutines
;               ties    [ADJUSTIE] tie abundance parameters 
;		temp	[LINEFLX] assume input LOGT are actually T [K]
;		noph	[LINEFLX] compute fluxes in [ergs/s/...] not [ph/s/...]
;		effar	[LINEFLX] effective area in [cm^2]
;		wvlar	[LINEFLX] wavelengths at which EFFAR are defined [Ang]
;		nhne	[LINEFLX] ratio of N(H)/n_e, assumed to be 0.83
;		softlim	[LIKELI] allow for upper limits to be softly defined
;		steps	[VARSMOOTH] assume stepped averaging from smoothing scale
;		weight	[VARSMOOTH] allow weighting of adjacent points in boxcar smooth
;		type	[VARSMOOTH] functional form to smooth with, see X3MODEL/LIBMODEL
;		boxcar	[VARSMOOTH] controls smoothing behavior at endpoints
;		betap	[MK_LORENTZ/MK_SLANT] index of beta-profile
;		angle	[MK_SLANT] angle of tilt
;		fwhm	[MK_GAUSS/MK_LORENTZ/MK_ROGAUSS] return FWHM rather than
;			sigma or core width
;		sloop	[LOOPEM] slope of loop DEM
;		dfx_mul	[GENERATIO] multiplicative delta_Fx to compute partials
;		dfx_add	[GENERATIO] additive delta_Fx to compute partials
;		choice	[FINDSCALE] flag to choose algo to find the scales
;		pick	[FINDSCALE] flag on hoe to collapse scales from 2D to 1D
;		crunch	[FINDSCALE] collapse 2D to 1D prior to finding scales
;		half	[FINDSCALE] return half-scale
;		eps	[FINDSCALE] a small number
;               mproxy  [MIXIE]  see MIXIE() 
;               mplet   [MIXIE]  see MIXIE()
;               
;restrictions
;	requires subroutines
;	-- KILROY
;	-- GETABUND
;	-- MCMC_DEM_ONLY [LINEFLX, LIKELI, MK_DEM, VARSMOOTH, LOOPEM]
;	-- LINEFLX [GETABUND, WHEE]
;	-- LIKELI [SOFTLIM]
;	-- FINDSCALE [WVLT_SCALE [ROOFN]]
;	-- MK_DEM
;	-- LOOPEM
;	-- VARSMOOTH
;	-- GENERATIO
;       -- ADJUSTIE
;	-- PEASECOLR
;	-- MCMC_DEM_WHISKERPLOT [EPSROOT,THRFRAC,THRNUM,DVDT,
;	   CBOX,CWHISK,CDOTS,CBEST,CMODE,CGREY,CNUM, plot keywords]
;
;history
;	vinay kashyap (Apr 97)
;	added keyword SMOOT (VK; May97)
;	added keyword LOOPY and replaced call to VARSMOOTH by call to
;	  MK_DEM (VK; Nov01)
;	added keyword ONLYRAT and call to GENERATIO (VK; Nov'02)
;	numerous bug fixes (VK/LL; Dec'02)
;	ratio stuff was not being passed on to MCMC_DEM_ONLY; now stores
;	  simulated DEMs, ABUNDs, fluxes, and ratios (VK/LL; Jan'03)
;	added keywords SAMPENV,STORPAR,STORIDX,SIMPRB,SIMDEM,SIMABN,
;	  SIMFLX,SIMPRD,NBURN (VK; Feb'03)
;	added keyword SMOOSCL (VK; May'03)
;	added keyword SPLINY (VK; Jun'03)
;	added _extra to call to PRED_FLX (VK; Jul'03)
;       bugfix sampenv ONLYRAT keyword handler included (VK/LL:Dec'03)
;       bugfix sampenv was ignored if input as array (LL: Dec'03)
;       bugfix nemt is undefined in logt keyword handler(LL:Dec'03)
;       added BRNSTR in call to MCMC_DEM_ONLY() to enable deblending
;         via MIXIE() (LL:Jan'04)
;       bugfix abndupdt, a MIXIE keyword  now set by default and 
;         turned off if abrng keyword is set (LL:Mar'04) 
;	added keyword XNSIG (VK; Jun'04)
;	bug fix: compute DEMs for /LOOPY and /SPLINY on correct T grid
;	  (VK; Jan05)
;	allow limiting number of T components for LOOPY (VK; Mar05)
;       added call to adjustie and MIXLSTR keyword (LL; Jun05) 
;	added LCOMP to call to MK_DEM, check to make sure FSIGMA
;	  are all +ve definite, and a warning if ABUND[0] is not 1,
;	  force single T deviate per iteration for LOOPY (VK; Jun05)
;	bug correction: if no update found in a batch, was still
;	  storing DEM from that batch; now corrected (VK; Jul05)
;	added options COVER and STEPPED for SAMPENV (VK; Jan06)
;	include extra kick at NPAR+1 to unfreeze parameters in case
;	  of LOOPY getting stuck; also bug correction for NPAR+1
;	  selections in regular case if some T bins are frozen
;	  (VK; Apr06)
;	added keyword VERBOSE (VK; Jun07)
;	added keyword BESTPAR (VK; Nov08)
;	added keywords SYSDEV and SYSGRP (VK; Apr09)
;	changed how stddev of proposal distr. is calculated, from using
;	  histograms to using unbinned distr.; added keyword STORDEM and
;	  now STORPAR outputs unsmoothed parameters (VK; Apr12)
;	added call to PEASECOLR and MCMC_DEM_WHISKERPLOT (VK; Sep15)
;-

;	usage
nw=n_elements(wvl) & nf=n_elements(flx) & sem=size(emis)
if nw eq 0 or nf eq 0 or sem(1) eq 0 then begin
  print,'Usage: dem=mcmc_dem(wavelengths,fluxes,emissivities,Z=Z,logt=logt,$'
  print,'  diffem=initial_DEM,abund=abundances,fsigma=flux_errors,$'
  print,'  onlyrat=use_fluxes_in_ratios,sampenv=sampenv,$'
  print,'  sysdev=sysdev,sysgrp=sysgrp,$'
  print,'  storpar=storpar,storidx=storidx,stordem=stordem,simprb=simprb,simdem=simdem,$'
  print,'  simabn=simabn,simflx=simflx,simprd=simprd,ulim=upper_limit_flag,$'
  print,'  demrng=DEMrange,abrng=abrng,smoot=smooth,smooscl=smooscl,$
  print,'  loopy=loopy,/spliny,verbose=verbose, xdem=DEM_bins,xab=ABUND_bins,xnsig=xnsig,$'
  print,'  nsim=nsim,nbatch=nbatch,nburn=nburn,/steps,weight=weight,$'
  print,'  /gauss,type=type,/boxcar,softlim=softlim,$'
  print,'  /nosrch,savfil=savfil,bound=bound,demerr=demerr,aberr=aberr,$
  print,'  /temp,/noph,effar=effar,wvlar=wvlar,/ikev,/regrid,$
  print,'  mproxy=mproxy,mplet=mplet, epsroot=epsroot,thrfrac=thrfrac,thrnum=thrnum,$
  print,'  /dvdt,cbox=cbox,cwhisk=cwhisk,cdots=cdots,cbest=cbest,cmode=cmode,cgrey=cgrey,cnum=cnum,$
  print,'  charthick=charthick,charsize=charsize,xrange=xrange,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title)'
  print,'returns best estimate of DEM(logT) and ABUND given flux@wavelengths'
  return,-1L
endif

;	cast variables
;  verbosity
vv=1L & if n_elements(verbose) ne 0 then begin
  if verbose[0] eq 0 then vv=0 else vv=long(verbose[0])>1
endif
;  wavelength, flux, line emissivity
ww=[wvl] & fx=[flx] & line=[emis]
;  atomic number
zz=intarr(nw)+26 & if keyword_set(Z) then begin
  if n_elements(z) le nw then zz(0)=Z else zz(*)=Z(0:nw-1)
endif
;  error on flux
sig=1.+sqrt(abs(flx)+0.75) & if keyword_set(fsigma) then begin
  if n_elements(fsigma) le nf then sig(0)=fsigma else sig(*)=fsigma(0:nf-1)
endif & sig=double(sig)
	;error in error?
o0=where(sig le 0,mo0)
if mo0 gt 0 then begin
  message,'some FSIGMA are below zero; replacing with Poisson errors',/informational
  sig(o0)=1.+sqrt(abs(flx(o0))+0.75)
endif
;  upper limits
if n_elements(ulim) ne nw then begin
  if n_elements(ulim) gt 0 then message,'unsetting upper limits',/informational
  ulim=lonarr(nw)
endif & oul=where(ulim eq 0,moul)
if moul gt 0 then tfx=total(fx(oul)) else tfx=0.
if moul gt 0 then tsig=sqrt(total(sig(oul)^2)) else tsig=1.
;  ratios
nrat=n_elements(onlyrat) & nrats=0 & idxflx=lindgen(nf) & obsdat=fx & obssig=sig
if nrat eq nw then begin		;(parse ONLYRAT
  generatio,fx,onlyrat,rx,fxerr=sig,rxerr=rsigma,verbose=vv, _extra=e
  nrats=n_elements(rx)
  if nrats gt 0 then begin
    idxflx=where(onlyrat eq '',midxflx)
    if midxflx gt 0 then begin
      obsdat=[fx(idxflx),rx] & obssig=[sig(idxflx),rsigma]
    endif else begin
      obsdat=[rx] & obssig=[rsigma]
    endelse
  endif
endif else begin			;ONLYRAT)(not kosher
  if keyword_set(onlyrat) then message,'ONLYRAT must be a string array of same size as FLUXES; ignoring.',/info
endelse					;ONLYRAT not kosher)
;  logT
nemt=sem(1)
if not keyword_set(logt) then begin
  if nemt gt 1 then logt=findgen(nemt)*(8.-4.)/(nemt-1)+4. else logt=[6.]
endif & logt=[logt] & nt=n_elements(logt)
;  DEM
dem=0*logt+1e14 & if keyword_set(diffem) then begin
  if n_elements(diffem) le nt then dem(0)=diffem else dem(*)=diffem(0:nt-1)
endif
;  abundances
nab=n_elements(abund) & abnd=getabund('anders & grevesse') & abmn=1e-20
if nab ne 0 then begin
  if nab le n_elements(abnd) then abnd(0:nab-1)=abund
  if nab gt n_elements(abnd) then abnd=abund(0:n_elements(abnd)-1)
endif
oo=where(abnd le 0) & if oo(0) ne -1 then abnd(oo)=abmn	;absolute min
if abnd(0) ne 1 then begin
  print,abnd
  message,'WARNING: H abundance is '+strtrim(abnd(0),2),/informational
  stop,'HALTing; type .CON to continue'
endif
;  number of simulations in each batch
if keyword_set(nsim) then nsim=(long(nsim(0))>10) else nsim=100L
;  number of batches
if keyword_set(nbatch) then nbatch=(long(nbatch(0))>5) else nbatch=10L
;  number of burn-ins
if keyword_set(nburn) then nburn=(long(nburn(0))) else nburn=long(nsim/10.+1)
;  systematic deviations
devsys=fltarr(nf) & nsd=n_elements(sysdev)
if nsd gt 0 then devsys[0L:(nsd<nf)-1L]=sysdev[0L:(nsd<nf)-1L]
grpsys=lonarr(nf)-1L & nsg=n_elements(sysgrp)
if nsg gt 0 then grpsys[0L:(nsg<nf)-1L]=sysgrp[0L:(nsg<nf)-1L]
ogp=where(grpsys eq 0,mogp) & if mogp gt 0 then grpsys[ogp]=max(grpsys)+lindgen(mogp)+1L
ugrp=grpsys[uniq(grpsys,sort(grpsys))] & nugrp=n_elements(ugrp)	;there are NUGRP groups

;	catch size mismatches
;nemt=sem(1) & ; nemt is needed above in logt handler
nemw=sem(2) & ok='ok'
if nw eq 1 and nf eq 1 and sem(0) eq 0 then begin	;single wvl, logT
  nemt=1 & nemw=1
endif
if nw gt 1 and sem(0) eq 1 and nw eq nemt then begin	;many wvls, single logT
  nemt=1 & nemw=nw
endif
if nemt eq 1 then begin					;reset logT and DEM
  if nt ne 1 then logt=[6.]
  nt=1L & dem=[dem(0)]
endif
if nw ne nf then ok='*** flux v/s wavelength mismatch ***'
if nw ne nemw then ok='*** line emissivity v/s wavelength mismatch ***'
if sem(0) gt 2 then ok='*** line emissivity array in strange format ***'
if strmid(ok,0,1) eq '*' then begin
  message,ok,/info & return,0.*dem+1.
endif

;	check parameter ranges
rngD=double([dem/1e5,dem*1e5])
if not keyword_set(demrng) then tmp=[rngD(*)] else tmp=[demrng(*)]
if n_elements(tmp) eq nt then begin
  oo=where(tmp gt 0,moo)
  if moo gt 0 then begin
    rngD(oo)=dem(oo)/tmp(oo) & rngD(oo+nt)=dem(oo)*tmp(oo)
  endif
  tmp=[rngD]
endif
if n_elements(tmp) le 2*nt then rngD(0)=tmp else rngD(*)=tmp(0:2*nt-1)
rngD=reform(rngD,nt,2)
;
nab=n_elements(abnd) & rngA=[abnd,abnd] & abndupdt = 1
if not keyword_set(abrng) then tmp=[rngA(*)] else begin 
 tmp=[abrng(*)] & abndupdt = 0 ;abndupdt toggles mixie() abundnace updating
endelse  
if n_elements(tmp) le 2*nt then rngA(0)=tmp else rngA(*)=tmp(0:2*nab-1)
rngA=reform(rngA,nab,2)

;	catch trivial errors
oo=where(dem le 0,moo) & if moo gt 0 then dem(oo)=1.		;(we'll
oo=where(rngD(*,0) le 0,moo) & if moo gt 0 then rngD(oo,0)=1.	;operate
oo=where(rngD(*,1) le 0,moo) & if moo gt 0 then rngD(oo,1)=1.	;in log)
oo=where(rngD(*,0) gt rngD(*,1),moo)
if moo gt 0 then begin						;LowLim > UL?
  tmp=rngD(oo,1) & rngD(oo,1)=rngD(oo,0) & rngD(oo,0)=tmp
endif
oo=where(abnd le 0,moo) & if moo gt 0 then abnd(oo)=abmn	;(we'll
oo=where(rngA(*,0) le 0,moo) & if moo gt 0 then rngA(oo,0)=abmn	;operate
oo=where(rngA(*,1) le 0,moo) & if moo gt 0 then rngA(oo,1)=abmn	;in log)
oo=where(rngA(*,0) gt rngA(*,1),moo)				;LL > UpLim?
if moo gt 0 then begin
  tmp=rngA(oo,1) & rngA(oo,1)=rngA(oo,0) & rngA(oo,0)=tmp
endif

;	figure out which parameters are to be thawed
allpar=[alog10(dem),alog10(abnd)] & ipar=lindgen(n_elements(allpar))
aparmx=[alog10(rngD(*,1)),alog10(rngA(*,1))]
aparmn=[alog10(rngD(*,0)),alog10(rngA(*,0))]
oo=where(aparmx eq aparmn,moo) & if moo gt 0 then ipar(oo)=-1
opar=where(ipar ge 0,npar)
if npar eq 0 then begin			;give the parameters their freedom
  message,'kuma kuma pende zende',/info
  return,dem
endif
opad=where(opar lt nt,npad) & if npad gt 0 then opad=opar(opad)

;	histogram bins
wdem=0.1 & wab=0.1
mindem=min(rngD,max=maxdem,/nan) & minab=min(rngA,max=maxab,/nan)
mindem=alog10(mindem) & maxdem=alog10(maxdem)
minab=alog10(minab) & maxab=alog10(maxab)
n_dem=n_elements(xdem) & n_ab=n_elements(xab)
if n_dem le 1 then begin
  if n_dem eq 0 then n_dem=long((maxdem-mindem)/wdem+0.5)+1 else $
	n_dem=long(xdem(0))>2
  xdem=findgen(n_dem)*wdem+mindem
endif
if n_ab le 1 then begin
  if n_ab eq 0 then n_ab=long((maxab-minab)/wab+0.5)+1 else $
	n_ab=long(xab(0))>2
  xab=findgen(n_ab)*wab+minab
endif
if max(xdem) lt maxdem then begin
  xdem=[xdem,max(xdem)+wdem] & n_dem=n_dem+1
endif
if max(xab) lt maxab then begin
  xab=[xab,max(xab)+wab] & n_ab=n_ab+1
endif

;	get the emissivity envelope
tmp=line & emenv=dblarr(nt) & for i=0,nt-1 do emenv(i)=total(tmp(i,*))
if keyword_set(sampenv) then begin
  nsampenv=n_elements(sampenv)
  szse=size(sampenv) & nsze=n_elements(szse)
  if szse(nsze-2L) eq 7 then begin
    cc=strupcase(sampenv(0))
    if strpos(cc,'ONLYRAT',0) ge 0 then begin 
     if szse(1) eq nt and nsze gt 1 and nrats gt 0 then begin
       emenv=dblarr(nt)
       for it=0L,nt-1L do begin
  	        generatio,reform(tmp(i,*)),onlyrat,rx,verbose=vv, _extra=e
        	if midxflx eq 0 then emenv(i)=total(rx) else $
		emenv(i)=total(rx)+total(tmp(i,idxflx))
       endfor
     endif   
    endif
    if strpos(cc,'INCABUN',0) ge 0 then begin
     emenv=dblarr(nt) & for i=0,nt-1 do emenv(i)=total(tmp(i,*)*abnd(zz-1))
    endif
    if strpos(cc,'FLAT',0) ge 0 then emenv(*)=1.
    if strpos(cc,'INCALLT',0) ge 0 then emenv=emenv+max(emenv)/100.
    if strpos(cc,'STEPPED',0) ge 0 then begin
      o1=where(emenv lt 0.25*max(emenv),mo1)
      o0=where(emenv le max(emenv)/1d15,mo0)
      if mo1 gt 0 then emenv(o1)=0.25*max(emenv)
      if mo0 gt 0 then emenv(o0)=0.
    endif
    if strpos(cc,'COVER',0) ge 0 then begin
      o0=where(emenv le max(emenv)/1d10,mo0) & emenv(*)=1.
      if mo0 gt 0 then emenv(o0)=0.
    endif
  endif else begin
    if szse(1) eq nt then emenv=sampenv*1.D
  endelse
endif
;;for i=0,nw-1 do tmp(*,i)=tmp(*,i)*abnd(zz(i)-1)		;incl. abund
;for i=0,nw-1 do tmp(*,i)=tmp(*,i)/max(tmp(*,i))		;w/o incl. abund
;emenv=dblarr(nt) & for i=0,nt-1 do emenv(i)=total(tmp(i,*))	;envelope
;;emenv=emenv+max(emenv)/100.		;force looks at "bad" points too
;
oo=where(aparmx(0:nt-1) eq aparmn(0:nt-1),moo)	;these T points frozen
if moo gt 0 then emenv(oo)=0.	;don't EVER consider frozen T points
oo=where(aparmx(0:nt-1) ne aparmn(0:nt-1),moo)	;these T points thawed
cenv=dblarr(moo)+emenv(0)
for i=1,moo-1 do cenv(i)=cenv(i-1)+emenv(oo[i])		;cumulative
cenv=(cenv/cenv(moo-1))>0

;	local smoothing scales
nsmtscl=n_elements(smooscl)
if nsmtscl eq 0 then lscal=findscale(line,_extra=e)>1 else begin
  if nsmtscl eq nt then lscal=smooscl else $
	lscal=intarr(nt)+smooscl[0]
endelse

;	reset the scales if necessary
if n_elements(smoot) ne 0 then lscal=long(lscal*smoot(0)+0.5)

;	initialize
isim=0L				;number of simulations
ksim=long(nsim)*long(npar)	;number of simulations saved
istor=0L			;index of saved simulation
iburn=0L			;batches used to burn-in
prb=0. 				;-log(probability)
ebound=0.9			;fraction of p(M|D) to include in error bar
if not keyword_set(bound) then ebound=0.9 else ebound=(double(bound(0))>0.1)<1
sigpar=(aparmx-aparmn)/5.	;std.dev. of parameters
if not keyword_set(xnsig) then xnsig=3.
;
parcnt=lonarr(npar)		;keep count of # times parameter is simulated
storpar=dblarr(ksim)		;store generated parameters
stordem=dblarr(ksim)		;store generated DEM values
storidx=lonarr(ksim)-1		;parameter index in above
;hdem=lonarr(nt,n_dem-1)		;histogram of all DEM values
;hab=lonarr(nab,n_ab-1)		;histogram of all ABUND values
bestpar=allpar & bestsig=sigpar			;begin with this
;
simprb=fltarr(nsim+1L)		;store NSIM "quality of fit" probabilities
				;and the best-fit PRB
simdem=dblarr(nt,nsim+1L)	;store NSIM simulations of DEMs and the best-fit
simabn=dblarr(nab,nsim+1L)	;store NSIM simulations of abundances and the best-fit
simflx=fltarr(nf,nsim+1L)	;store NSIM simulations of predicted fluxes and the
 				;flux produced by the best-fit 
if keyword_set(nrats) then simprd=fltarr(n_elements(obsdat),nsim+1L)
				;store NSIM simulations of predicted fluxes and ratios
				;and those produced by the best-fit

;       Check on mixsstr
call_mixie = 0 
if keyword_set(mixsstr) then begin 
  call_mixie=1 
  ;     construct lstr for mixie from wvl,emis,z 
if not keyword_set(mixlstr) then  $
  mixlstr=create_struct('LINE_INT',emis,'LOGT',logt,'WVL',wvl,'z',z,'ION',0*z, $
  'DESIG', strarr(2,nw),'CONFIG',strarr(2,nw),'SRC',0*z, 'JON',0*z)
  xfrac = mixie(mixlstr, mixsstr, dem=dem, dlogt = logt, obsflx=flx, ofsig=fsigma,brnstr=brnstr, abund=abnd,abndupdt=abndupdt,_extra=e)
endif else xfrac=0*wvl*1.0

;	figure out best regions to search in
if not keyword_set(nosrch) then begin		;{do search
  for i=0,npar-1 do begin			;(shtep thru parameters
    j=opar(i) & kilroy,dot='('+strtrim(j,2)+')'
    if j lt nt then xx=xdem else xx=xab
    nx=n_elements(xx) & prbpar=fltarr(nx)
    for k=0,nx-1 do begin			;{step thru range
      tmpar=allpar & tmpar(j)=xx(k)
      mcmc_dem_only,tmpar,fx,sig,ulim,line,logT,ww,zz,lscal,ff,prb,dem,abnd,$
	/chi2,/rchi,loopy=loopy,spliny=spliny,nrats=nrats,onlyrat=onlyrat,$
	obsdat=obsdat,obssig=obssig,idxflx=idxflx,prddat=prddat,$
	mixlstr=mixlstr,mixsstr=mixsstr,brnstr=brnstr,abndupdt=abndupdt,$
	xfrac=xfrac,verbose=vv,_extra=e
	;compute predicted flux and likelihood
      tf=total(ff)
      prbpar(k)=prb
    endfor					;k=0,nx-1}
    ;if i ge nt then stop
    prbpar=-prbpar & prbpar=((prbpar-max(prbpar))<69)>(-69)
    prbpar=exp(prbpar) & prbpar=prbpar/total(prbpar)
    oo=where(prbpar eq max(prbpar),moo)
    if moo gt 1 then begin
      bestpar(j)=total(prbpar(oo)*xx(oo))/total(prbpar(oo))
    endif else begin
      if moo eq 1 then bestpar(j)=xx(oo)
    endelse
    minsig=xnsig*abs(min(xx(1:*)-xx))
    bestsig(j)=(total(prbpar*xx*xx)-(total(prbpar*xx))^2) > minsig^2
    bestsig(j)=sqrt(bestsig(j))
    if finite(bestsig(j)) eq 0 then bestsig(j)=(max(xx)-min(xx))/2.
    kilroy,dot=strtrim(bestpar(j),2)+'+-'+strtrim(bestsig(j),2)+' '
    allpar(j)=bestpar(j)
    ;if i ge nt then stop
  endfor					;i=0,npar-1)
endif else begin				;}{don't waste time...
  szsrch=size(nosrch) & nszsrch=n_elements(szsrch)
  if szsrch(nszsrch-2) eq 7 then restore,nosrch(0)
endelse						;NOSRCH}
allpar=bestpar & sigpar=bestsig
if vv gt 0 then begin
  window,0 & !p.multi=[0,2,3] & window,1 & wset,0
endif

bestpar=allpar & bestsig=sigpar		;begin with this
;	ok, what's the bad news?  how far off are we?
dem=[allpar(0:nt-1)] & abnd=[allpar(nt:*)]
dem=10.D^(dem) & abnd=10.^(abnd)

if call_mixie then xfrac=mixie(mixlstr,mixsstr,dem=dem, dlogt = logt, obsflx = flx,brnstr=brnstr, abund=abund,abndupdt=abndupdt,_extra=e)

if nt eq 1 then begin		;because LINEFLX can't handle LINE(WVL)
  ff=fltarr(nw)
  for iw=0,nw-1 do ff(iw)=lineflx(line(iw),logT,ww(iw),zz(iw),DEM=DEM,$
	abund=abnd, _extra=e)/xfrac(IW)
endif else ff=lineflx(line,logT,ww,zz,DEM=DEM,abund=abnd, _extra=e)/xfrac
tf=total(ff)
oldprb=likeli(fx,ff,dsigma=sig,/chi2,ulim=ulim,_extra=e)
bestprb=oldprb & kilroy,dot=' ['+strtrim(bestprb,2)+']'

;	begin simulation
while isim lt nsim do begin			;{start

  kilroy,dot='('+strtrim(isim,2)+')'
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if vv gt 0 then plot,allpar(opar),yr=[min(aparmn(opar)),max(aparmx(opar))],col=150,$
		charsize=1.7
	if vv gt 0 then for i=0,npar-1 do oplot,[i,i],allpar(opar(i))+$
		sigpar(opar(i))*[1,-1],col=100
  for ibtc=0,nbatch-1 do begin			;(one batch
    rn=randomn(seed,npar)
    ru=randomu(seed,2*npar+1)
    hitpar=lonarr(npar)
    nnpar=npar+1L ;& if keyword_set(loopy) then nnpar=npar
    ;	NOTE: the way this works normally is, go through NPAR
    ;	steps, picking a random T param at each time.  at the
    ;	end of NPAR picks, usually there will be some that have
    ;	not been picked.  pick all those at one shot at the
    ;	NPAR+1'th step.  this bit was originally removed for
    ;	loopy because what's the point, ey?  but it turns out
    ;	that loopy gets stuck in local minima without recourse
    ;	to getting out.  so, in the NPAR+1'th time, simply try
    ;	the case with all the already picked parameters
    loopar=bestpar
    for jpar=0,nnpar-1 do begin			;{NPAR+1 steps/batch (unless LOOPY)
      tmpar=allpar				;change from
      ;	this is a hack to force loopy out of the comfortable
      ;	local minima it keeps finding itself in
      if keyword_set(loopy) then tmpar(0:nt-1)=$
	shift(varsmooth(allpar(0:nt-1),lscal),fix(randomn(seed)))
      if jpar lt npar then begin
        if opar(jpar) lt nt then begin		;pick parameter to change
	  ;kpar=jpar
	  oo=where(cenv ge ru(jpar)) & kpar=oo(0)
	  hitpar(kpar)=1
        endif else kpar=jpar
      endif else kpar=jpar

      ;	make the change
      if kpar lt npar then begin		;(NPAR steps
	tmpar(opar(kpar))=tmpar(opar(kpar))+rn(jpar)*sigpar(opar(kpar))
	lcomp=opar(kpar)
	loopar(opar(kpar))=tmpar(opar(kpar))
      endif else begin				;)(get those that were missed
	ono=where(hitpar(opar) eq 0,mno)
	if mno gt 0 then tmpar(opar(ono))=tmpar(opar(ono))+$
		randomn(seed,mno)*sigpar(opar(ono))
	oya=where(hitpar ne 0,mya)
	if keyword_set(loopy) then begin
	  lcomp=-1L ;& if mno gt 0 then lcomp=opar(ono)
	  if mya gt 0 then tmpar(oya)=loopar(oya)
	endif
      endelse					;KPAR=NPAR)
      ;	and force compliance with hard limits
      oo=where((tmpar gt aparmx or tmpar lt aparmn) and ipar ge 0,moo)
      while moo gt 0 do begin
        tmpar(oo)=allpar(oo)+randomn(seed,moo)*sigpar(oo)
        oo=where((tmpar gt aparmx or tmpar lt aparmn) and ipar ge 0,moo)
      endwhile
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if kpar lt npar and vv gt 0 then oplot,[kpar],[tmpar(opar(kpar))],psym=3,col=150

      ; add jitter to the emissivities
      jitline=line
      if nsd gt 0 then begin		;(jitter EMIS
	rsys=10.^(randomn(seed,nf)*devsys)
	for igp=0L,nugrp-1L do begin
	  ogp=where(grpsys eq ugrp[igp],mogp) & if mogp eq 0 then message,'BUG!'
	  for kgp=0L,mogp-1L do jitline[*,ogp[kgp]]=$
		jitline[*,ogp[kgp]]*10.^(rsys[ogp[0]]*devsys[ogp[kgp]])
	endfor
      endif				;NSD>0)

      ;	compute predicted flux and likelihood  
      mcmc_dem_only,tmpar,fx,sig,ulim,jitline,logT,ww,zz,lscal,tmpff,prb,tmpdem,tmpabnd,$
	/chi2,loopy=loopy,lcomp=lcomp,spliny=spliny,nrats=nrats,onlyrat=onlyrat,$
	obsdat=obsdat,obssig=obssig,idxflx=idxflx,prddat=tmpprd, mixlstr=mixlstr,mixsstr=mixsstr,$
        xfrac=xfrac,brnstr=brnstr,abndupdt=abndupdt,verbose=vv, _extra=e ;compute predicted flux and likelihood

      ;	NOTE: when priors are implemented, multiply exp(-prb) by the product of
      ;	the priors of the parameters here.
      ;	*bad* e.g. -- match the total counts in observed v/s model spectrum:
      ;if moul gt 0 then tf=total(tmpff(oul)) else tf=0.
      ;prb=prb+(tfx-tf)^2/tsig^2/2.

      ;	Metropolis check -- keep new value for parameter or not?
      ;	(remember that PRB is actually CHI^2)
      ;	{accept_fn=(prb/oldprb)<1 & metro_check=randomu(seed)}
      	;NOTE:	if new pars make for a better fit, PRB will be <OLDPRB
	;	(because PRB is -log(likelihood)) so ACCEPT_FN will be
	;	capped at 0.  Then, any value of metro_check, which runs
	;	from -\infty to <0 will fit the condition, so it will
	;	always be accepted.
	;	if the new pars make for a very poor fit, then ACCEPT_FN
	;	will be a large negative number, and the probability of
	;	accepting is appropriately reduced since there is a smaller
	;	chance that METRO_CHECK will now meet the condition.
      accept_fn=(oldprb-prb)<0 & metro_check=alog(ru(jpar+npar))

      if metro_check lt accept_fn then begin	;accept new value

	kilroy,dot='*'+strtrim(prb,2)
        allpar=tmpar & dem=tmpdem & oldprb=prb & abnd=tmpabnd
	ff=tmpff & prddat=tmpprd
	
	if vv gt 0 then begin
	  if kpar lt npar then oplot,[kpar],[allpar(opar(kpar))],psym=3 else $
		oplot,allpar(opar),col=150
	  if kpar lt npar then oplot,alog10(dem),col=1,psym=10 else $
		oplot,alog10(dem),col=3,psym=10
	endif
      endif

      ;	have we found the best parameters yet?
      if prb lt bestprb then begin	;remember that PRB is really CHI^2
	bestprb=prb & bestpar=allpar & bestdem=dem & bestabnd=abnd
	bestflx=ff
	if keyword_set(nrats) then bestprd=prddat
	kilroy,dot=' ['+strtrim(bestprb,2)+']'
	simprb(nsim)=bestprb
	simdem(*,nsim)=bestdem
	simabn(*,nsim)=bestabnd
	simflx(*,nsim)=bestflx
	if keyword_set(nrats) then simprd(*,nsim)=bestprd
	;stop,prb
      endif

      ;	store every NBATCH'th simulation after stabilization
      if kpar lt npar then begin
	storpar(istor)=tmpar(opar(kpar)) & storidx(istor)=opar(kpar)
        ;	**************
        ;	temporarily...
        ;	**************
        if opar(kpar) lt nt then begin
          ko=opar(kpar)
          stordem(istor)=alog10(abs(dem(ko)))<aparmx(ko)>(aparmn(ko))
        endif
        ;	**************
        ;	temporarily...
        ;	**************
        if iburn eq nburn then begin
          parcnt(kpar)=parcnt(kpar)+1
          if parcnt(kpar) eq long(parcnt(kpar)/nbatch)*nbatch then $
		istor=istor+1L
        endif else istor=istor+1L
      endif

    endfor					;JPAR=0,NNPAR-1}
  endfor					;ibtc=0,nbatch-1)
  simprb(isim)=oldprb
  simdem(*,isim)=dem
  simabn(*,isim)=abnd
  simflx(*,isim)=ff
  if keyword_set(nrats) then simprd(*,isim)=prddat
  dem=[allpar(0:nt-1)] & abnd=[allpar(nt:*)]
  ;dem=varsmooth(dem,lscal,_extra=e)				;smooth
  if keyword_set(loopy) then $
    dem=simdem(*,isim) else $
    dem=alog10(simdem(*,isim))
    ;dem=mk_dem('loop',indem=10.D^(dem),pardem=logT,logT=logT,loopy=loopy,$
    ;	lcomp=lcomp, _extra=e) else $
    ;if keyword_set(spliny) then $
    ;  dem=mk_dem('spline',logT=logT,indem=dem,pardem=logT,logT=logT, _extra=e) else $
    ;    dem=mk_dem('varsmooth',logT=logT,indem=dem,pardem=lscal, _extra=e)
  if not keyword_set(loopy) then dem=10.D^(dem) & abnd=10.^(abnd)
  adjustie, abnd, _extra = e 
  if vv gt 0 then oplot,alog10(dem),psym=10,col=4
  if vv gt 0 then wset,1
  tmp=pred_flx(line,logT,ww,dem,Z=zz,fobs=fx,fsigma=sig,abund=abnd,ulim=ulim,verbose=vv,_extra=e)
  if vv gt 0 then wset,0
  kilroy,dot=string(likeli(fx,tmp,dsigma=sig,/chi2,ulim=ulim, _extra=e))

  if iburn lt nburn then begin		;{burn-in until "stabilization"
    jstor=istor-nbatch*npar & norm=0L & tmpsig=sigpar
    if jstor lt 1 then dsig=0.D else begin	;(how different is new set?
      opars=storpar(0:jstor-1) & pars=storpar(0:istor-1)
      opari=storidx(0:jstor-1) & pari=storidx(0:istor-1)
      for i=0,npar-1 do begin		;(shtep through parameters
	j=opar(i)
	oo0=where(opari eq j,moo0) & oo1=where(pari eq j,moo1)
	if moo1 lt moo0 then message,'	insect infestation!'
	if moo0 gt 0 then begin		;(any reason to bother?
	  if j lt nt then xx=xdem else xx=xab
	  minx=min(xx,max=maxx) & xx=0.5*(xx(1:*)+xx)
	  f0=opars(oo0) & f1=pars(oo1) & binx=abs(xx(1)-xx(0))
	  hf0=float(histogram(f0,min=minx,max=maxx,binsize=binx))
	  hf1=float(histogram(f1,min=minx,max=maxx,binsize=binx))
	  hf0=hf0/total(hf0) & hf1=hf1/total(hf1)
	  dsig=dsig+total(hf0*hf1)^2/total(hf0^2)/total(hf1^2) & norm=norm+1
	  ;
	  os=sort(f1) & f1=f1(os) & xf1=findgen(moo1)/float(moo1)
	  if moo1 gt 3 then begin
	    tmp=interpol(f1,xf1,[0.16,0.84])
	    if tmp(0) lt tmp(1) then tmpsig(j)=(3.*(tmp(1)-tmp(0))) > (0.3*binx) else $
	  	tmpsig(j)=3*binx ;in case the distribution is too narrow,
				 ;or there are too few values
	  endif else tmpsig(j)=(max(f1)-min(f1))>(3.0*binx)
	  ;{VK Apr2012:	the code below calculated the stddev via binned histograms. so mediveal.
	  ;mpar=total(hf1*xx) & vpar=total(hf1*xx*xx)-mpar^2
	  ;tmpsig(j)=sqrt(vpar)>(xnsig*binx)		;new sigma
	  ;if finite(tmpsig(j)) eq 0 then tmpsig(j)=(maxx-minx)/2.
	  ;
	  ;cf1=0*xx+hf1(0)
	  ;for ii=1,n_elements(xx)-1 do cf1(ii)=cf1(ii-1)+hf1(ii)
	  ;cf1=cf1/max(cf1) & oh=where(cf1 ge 0.16 and cf1 le 0.84,moh)
	  ;if moh gt 0 then begin			;new sigmas/pars
	  ;  tmpsig(j)=0.5*(xx(oh(moh-1))-xx(oh(0))) > (3*binx)
	  ;endif else begin
	  ;  tmpsig(j)=3*binx
	  ;endelse
	  ;VK Apr2012}
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	  ;plot,xx,hf0,col=150,title=strtrim(j,2)+' '+strtrim(mpar,2)+' '+$
	  ;  strtrim(tmpsig(j),2),/xs,charsize=2
	  ;oplot,xx,hf1,line=1 & wait,0.5
	endif				;moo0>0)
      endfor				;I=0,NPAR-1)
    endelse					;jstor>0)
    dsig=dsig/(norm>1)					;comparison metric
    kilroy,dot='<'+strtrim(dsig,2)+'>'
    if dsig gt 0.98 then iburn=nburn else iburn=iburn-1	;stopping rule #1
    if abs(iburn) eq nburn then iburn=nburn		;stopping rule #2
    if istor ge ksim then iburn=nburn			;stopping rule #3
    if iburn eq nburn then begin
      istor=0			;discard all saved parameters
      storpar(*)=0.D & storidx(*)=-1L
      sigpar=tmpsig		;reset sigmas
    endif
  endif else isim=isim+1		;then increment}
endwhile						;ISIM<NSIM}

;	outputs
dem=reform(simdem(*,nsim)) & demerr=dblarr(nt,2)
for i=0L,nt-1L do begin
  demerr(i,*)=dem(i)
  tmp=reform(simdem(i,0L:nsim)) & os=sort(tmp) & tmp=tmp(os)
  if total(tmp) gt 0 then begin
    ou=where(tmp ge dem(i),mou) & tmpu=tmp[ou]
    ol=where(tmp le dem(i),mol) & tmpl=tmp[ol]
    if mol gt 1 then demerr(i,0)=interpol(tmpl,reverse(findgen(mol)),ebound*mol)
    if mou gt 1 then demerr(i,1)=interpol(tmpu,findgen(mou),ebound*mou)
  endif
endfor
;
abund=reform(simabn(*,nsim)) & aberr=fltarr(nab,2)
for i=0L,nab-1L do begin
  aberr(i,*)=abund(i)
  tmp=reform(simabn(i,0L:nsim)) & os=sort(tmp) & tmp=tmp(os)
  if total(tmp) gt 0 then begin
    ou=where(tmp ge abund(i),mou) & tmpu=tmp[ou]
    ol=where(tmp le abund(i),mol) & tmpl=tmp[ol]
    if mol gt 1 then aberr(i,0)=interpol(tmpl,reverse(findgen(mol)),ebound*mol)
    if mou gt 1 then aberr(i,1)=interpol(tmpu,findgen(mou),ebound*mou)
  endif
endfor
;
sflx=reform(simflx(*,nsim)) & sflxerr=dblarr(nf,2)
for i=0L,nf-1L do begin
  sflxerr(i,*)=sflx(i)
  tmp=reform(simflx(i,0L:nsim)) & os=sort(tmp) & tmp=tmp(os)
  if total(tmp) gt 0 then begin
    ou=where(tmp ge sflx(i),mou) & tmpu=tmp[ou]
    ol=where(tmp le sflx(i),mol) & tmpl=tmp[ol]
    if mol gt 1 then sflxerr(i,0)=interpol(tmpl,reverse(findgen(mol)),ebound*mol)
    if mou gt 1 then sflxerr(i,1)=interpol(tmpu,findgen(mou),ebound*mou)
  endif
endfor
;
if keyword_set(nrats) then begin
  sprd=reform(simprd(*,nsim)) & sprderr=dblarr(nf,2)
  nfp=n_elements(obsdat)
  for i=0L,nfp-1L do begin
    sprderr(i,*)=sprd(i)
    tmp=reform(simprd(i,0L:nsim)) & os=sort(tmp) & tmp=tmp(os)
    if total(tmp) gt 0 then begin
      ou=where(tmp ge sprd(i),mou) & tmpu=tmp[ou]
      ol=where(tmp le sprd(i),mol) & tmpl=tmp[ol]
      if mol gt 1 then sprderr(i,0)=interpol(tmpl,reverse(findgen(mol)),ebound*mol)
      if mou gt 1 then sprderr(i,1)=interpol(tmpu,findgen(mou),ebound*mou)
    endif
  endfor
endif

;	show a plot of the DEM with "errors"
if vv gt 0 then wset,0 & !p.multi=0
tvlct,rr,gg,bb,/get & peasecolr	;set color table
mcmc_dem_whiskerplot,'none',logT=logT,simdem=simdem,storidx=storidx, _extra=e
tvlct,rr,gg,bb			;reset color table
;
;{	this was the original plotting code.  whisker plot above is wayyy better
;icols=255-fix(255*(simprb-min(simprb))/(max(simprb)-min(simprb)))
;os=sort(icols)
;plot,logt,simdem(*,nsim),psym=10,/ylog,/ynoz,xtitle='logT',ytitle='DEM'
;for i=0L,nsim-1L do oplot,logt,simdem(*,os(i)),psym=10,col=icols(os(i))
;oplot,logt,simdem(*,nsim),psym=10,thick=5
;}

;;	outputs
;;mapar=double(varsmooth(allpar,lscal,_extra=e))	;MAP estimates of parameters
;if not keyword_set(loopy) then $
;  mapar=double(mk_dem('varsmooth',logT=logT,indem=allpar(0:nt-1),pardem=lscal,_extra=e)) else $
;  mapar=double(mk_dem('loop',indem=allpar(0:nt-1),pardem=logT,_extra=e))	;MAP estimates of parameters
;;mapar=double(allpar)	;MAP estimates of parameters
;upar=mapar				;upper confidence level
;upar=allpar				;upper confidence level
;lpar=upar				;lower confidence level
;fupar=upar				;upper HalfWidth@HalfMax
;flpar=lpar				;lower HalfWidth@HalfMax
;for i=0,npar-1 do begin			;{shetp thru parameters
;  j=opar(i)
;  if j lt nt then x=xdem else x=xab
;  if j lt nt then xx=0.5*(xdem(1:*)+xdem) else xx=0.5*(xab(1:*)+xab)
;  ;
;  oh=where(storidx eq j,moh)
;  if moh gt long(0.25*nsim) then begin	;(any reason to bother?
;    g=storpar(oh) & nx=n_elements(xx) & gf=fltarr(nx)
;    for k=0,nx-1 do begin		;probability distribution for parameter
;      oo=where(g ge x(k) and g lt x(k+1),moo) & gf(k)=moo
;    endfor
;    f=gf
;    ;smth=long(0.1*nx)>3 & f=smooth(gf,smth,/edge)
;    f=float(f)/total(f)
;    cf=f & for k=1,n_elements(f)-1 do cf(k)=cf(k-1)+f(k)
;    cf=cf/max(cf)				;cumulative distribution
;    ;
;    tmp=max(f,ip) & mapar(j)=xx(ip)	;MAP estimate
;    ;
;    ;	error bounds
;    ebu=(cf(ip)+ebound/2.)<1 & ebl=(cf(ip)-ebound/2.)>0
;    ipu=min(where(cf ge ebu)) & ipl=max(where(cf le ebl))>0
;    if ipu eq -1 then ipu=nx-1
;    upar(j)=xx(ipu) & lpar(j)=xx(ipl)
;    ;
;    ;	HWHMs
;    oo=where(f ge f(ip)/2.)
;    ifl=min(oo,max=ifu)
;    fupar(j)=xx(ifu) & flpar(j)=xx(ifl)
;  endif					;MOO>0.25*NSIM)
;endfor					;I=0,NPAR-1}
;
;;	recast outputs
;;dem=mapar(0:nt-1) & abund=mapar(nt:*)
;dem=allpar(0:nt-1) & abund=allpar(nt:*)
;demerr=dblarr(nt,2) & demerr(*,0)=lpar(0:nt-1) & demerr(*,1)=upar(0:nt-1)
;aberr=fltarr(nab,2) & aberr(*,0)=lpar(nt:*) & aberr(*,1)=upar(nt:*)
;if keyword_set(hwhm) then begin
;  demerr(*,0)=flpar(0:nt-1) & demerr(*,1)=fupar(0:nt-1)
;  aberr(*,0)=flpar(nt:*) & aberr(*,1)=fupar(nt:*)
;endif
;dem=10.D^(dem) & abund=10.^(abund) & demerr=10.D^(demerr) & aberr=10.^(aberr)

if keyword_set(savfil) then begin	;save variables
  szsv=size(savfil) & nszsv=n_elements(szsv)
  if szsv(nszsv-2) ne 7 or szsv(0) ge 1 then svfil='mcmc.save' else svfil=savfil
  message,'saving variables to IDL save file '+svfil,/info
  save,file=savfil
endif

return,dem
end
