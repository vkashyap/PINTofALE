function mcmc_dem_post,savfil,clev=clev,verbose=verbose, _extra=e
;+
;function	mcmc_dem_post
;	returns the best-fit values and the error ranges for the
;	various parameters in a structure of the form:
;	{LOGT, NUEFF, MINCHISQ, BOUND,$
;	DEM_BEST,DEM_MODE,DEM_HPD,DEM_EQT,DEM_MEDIAN,DEM_RNG,DEM_ENV,DEM_FRAC,$
;	AB_BEST,AB_MODE,AB_HPD,AB_EQT,AB_MEDIAN,AB_RNG,AB_ENV}
;
;syntax
;	mcmcstr=mcmc_dem_post(savfil,clev=clev,verbose=verbose)
;
;parameters
;	savfil	[INPUT; required] full pathnmae to the save file
;		that is written out by MCMC_DEM()
;
;keywords
;	clev	[INPUT] set this to override the confidence level at
;		which to determine the bounds on the parameters
;		* default is to use EBOUND in SAVFIL
;		* if < 0, abs value is used
;		* if > 1 and < 100, then assumed to be given as a percentage
;		* if > 100, then 1-1/CLEV is used as the true value
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;output description
;	LOGT:	the temperature grid
;	NUEFF:	the effective degrees of freedom for adopted smoothing
;	MINCHISQ:	the minimum chisq achieved
;	BOUND:	the confidence level at which to compute errors
;	*_BEST: best-fit parameter value
;	*_MODE: mode of parameter
;	*_MEDIAN:	median of parameter
;	*_HPD: highest probability-density range at confidence BOUND
;	*_EQT: equal-tail range at confidence BOUND
;	*_RNG: half-tail range at confidence BOUND measured from best-fit
;	*_ENV:	envelope range of best BOUND fraction (by chisq)
;	DEM_FRAC:	number of hits at each temperature bin
;	HELP:	a description of the output
;
;subroutines
;	MODALPOINT
;	HIPD_INTERVAL
;	EQT_INTERVAL
;
;history
;	vinay kashyap (Apr2006)
;	bug fix: simprb was being compared to chisqcvf rather than half that
;	  (LL; May2006)
;-

;	usage
ok='ok' & np=n_params() & nf=n_elements(savfil)
sz=size(savfil) & nsz=n_elements(sz)
if np eq 0 then ok='Insufficient parameters' else $
 if nf eq 0 then ok='SAVFIL is undefined' else $
  if nf gt 1 then ok='SAVFIL cannot be an array' else $
   if sz[nsz-2] ne 7 then ok='SAVFIL is not a character string'
if ok ne 'ok' then begin
  print,'Usage: mcmcstr=mcmc_dem_post(savfil,clev=clev,verbose=verbose)'
  print,'  computes and returns best-fit values and error ranges'
  print,'  based on MCMC_DEM simulations'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	verbosity
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	restore the savfil
fil=findfile(savfil,count=nfil)
if nfil eq 0 then begin
  message,SAVFIL[0]+': file not found',/informational
  return,-1L
endif
restore,savfil[0],verbose=(vv gt 500)
;	restores LOGT, SIMPRB, EBOUND, SIMDEM, SIMABN, STORIDX
;	and also NT, NAB, NSIM, LSCAL, LOOPY, SPLINY, RNGA, RNGD

;	figure out preliminary output
;NUEFF
nueff=nT
if keyword_set(loopy) then begin
    if loopy gt 1 then nueff=fix(loopy)
  endif else begin
    if not keyword_set(spliny) then $
	tmp=varsmooth(dem,lscal,nueff=nueff,nutilde=nutilde, _extra=e)
endelse
;MINCHISQ
minchisq=min(simprb)*2
;EBOUND
credlev=ebound & if keyword_set(clev) then credlev=0.0+clev[0]
if credlev lt 0 then credlev=abs(credlev)
if credlev gt 1 and credlev lt 100 then credlev=credlev/100.
if credlev gt 100 then credlev=1.0D - 1.0D/credlev
;HELP
help=[	'These are the summary statistics from the MCMC_DEM run',$
	'stored in the IDL save file '+savfil,$
	'',$
	'The variables are:',$
	'LOGT:   	the temperature grid',$
	'NUEFF:   	effective degrees of freedom if smoothing',$
	'MINCHISQ:	minimum chisq achieved',$
	'BOUND:   	confidence level at which to compute errors',$
	'*_BEST: 	best-fit parameter value',$
	'*_MODE: 	mode of parameter',$
	'*_MEDIAN:	median of parameter',$
	'*_HPD:  	highest probability-density range at confidence BOUND',$
	'*_EQT:  	equal-tail range at BOUND',$
	'*_RNG:  	half-tail range at BOUND measured from best-fit',$
	'*_ENV:  	envelope range of best BOUND fraction (by chisq)',$
	'DEM_FRAC:	number of hits at each temperature bin',$
	'HELP:   	this description']

;	the "best" chisqs
chisqrcvf=chisqr_cvf(1.-credlev,nueff)
ochi=where(simprb[0L:nsim-1L] le chisqrcvf/2.,mochi)
if mochi lt credlev*nsim*0.5 then begin
  message,'WARNING ***	this looks like a bad fit!',/informational
  ochi=lindgen(long(credlev*nsim+0.5))
endif
osim=sort(simprb)

;	figure out the best-fit abundances
ab_best=reform(SIMABN[*,NSIM]) & ab_mode=ab_best & ab_median=ab_best
ab_hpd=fltarr(NAB,2) & ab_hpd[*,0]=ab_best & ab_hpd[*,1]=ab_best
ab_eqt=ab_hpd & ab_rng=ab_hpd & ab_env=ab_hpd
for i=0L,NAB-1L do begin
  if RNGA[i,0] ne RNGA[i,1] then begin
    ab_mode[i]=modalpoint(simabn[i,0L:nsim-1L])
    ab_median[i]=median(simabn[i,0L:nsim-1L])
    ab_hpd[i,*]=hipd_interval(simabn[i,0L:nsim-1L],/fsample,clev=credlev)
    ab_eqt[i,*]=eqt_interval(simabn[i,0L:nsim-1L],/fsample,clev=credlev)
    ab_rng[i,*]=eqt_interval(simabn[i,0L:nsim-1L],/fsample,clev=credlev,$
	xaround=ab_best[i])
    ab_env[i,*]=minmax(simabn[i,osim[ochi]])
    if vv gt 100 then print,i,[ab_best[i],ab_mode[i],ab_median[i]],$
	reform(ab_hpd[i,*]-ab_mode[i]),$
	reform(ab_eqt[i,*]-ab_median[i]),$
	reform(ab_rng[i,*]-ab_best[i]),$
	reform(ab_env[i,*])
  endif
endfor

;	figure out the best-fit DEMs
DEM_best=reform(SIMDEM[*,NSIM]) & DEM_mode=DEM_best & DEM_median=DEM_best
DEM_frac=lonarr(NT)
DEM_hpd=fltarr(NT,2) & DEM_hpd[*,0]=DEM_best & DEM_hpd[*,1]=DEM_best
DEM_eqt=DEM_hpd & DEM_rng=DEM_hpd & DEM_env=DEM_hpd
for i=0L,NT-1L do begin
  if RNGD[i,0] ne RNGD[i,1] then begin
    DEM_mode[i]=modalpoint(simDEM[i,0L:nsim-1L])
    DEM_median[i]=median(simDEM[i,0L:nsim-1L])
    DEM_hpd[i,*]=hipd_interval(simDEM[i,0L:nsim-1L],/fsample,clev=credlev)
    DEM_eqt[i,*]=eqt_interval(simDEM[i,0L:nsim-1L],/fsample,clev=credlev)
    DEM_rng[i,*]=eqt_interval(simDEM[i,0L:nsim-1L],/fsample,clev=credlev,$
	xaround=DEM_best[i])
    DEM_env[i,*]=minmax(simDEM[i,osim[ochi]])
    ok=where(storidx eq i,mok) & DEM_frac[i]=mok
    if vv gt 200 then print,i,DEM_frac[i],$
	[DEM_best[i],DEM_mode[i],DEM_median[i]],$
	reform(DEM_hpd[i,*]-DEM_mode[i]),$
	reform(DEM_eqt[i,*]-DEM_median[i]),$
	reform(DEM_rng[i,*]-DEM_best[i]),$
	reform(DEM_env[i,*])
  endif
endfor

;	define the output structure
mcmcstr=create_struct('LOGT',LOGT, 'NUEFF',NUEFF, 'MINCHISQ', MINCHISQ,$
	'BOUND',CREDLEV, 'DEM_BEST',DEM_BEST, 'DEM_MODE',DEM_MODE,$
	'DEM_HPD',DEM_HPD, 'DEM_EQT',DEM_EQT, 'DEM_MEDIAN',DEM_MEDIAN,$
	'DEM_RNG',DEM_RNG, 'DEM_ENV',DEM_ENV, 'DEM_FRAC',DEM_FRAC,$
	'AB_BEST',AB_BEST, 'AB_MODE',AB_MODE, 'AB_HPD',AB_HPD,$
	'AB_EQT',AB_EQT, 'AB_MEDIAN',AB_MEDIAN, 'AB_RNG',AB_RNG,$
	'AB_ENV',AB_ENV, 'HELP',HELP)

return,mcmcstr
end
