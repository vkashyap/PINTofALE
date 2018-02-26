pro mcmc_dem_only,pars,fobs,fsigma,ulim,line,logT,wvl,Z,scale,$
    ff,prb,dem,abnd,loopy=loopy,lcomp=lcomp,spliny=spliny,nrats=nrats,$
    onlyrat=onlyrat,obsdat=obsdat,obssig=obssig,idxflx=idxflx,$
    prddat=prddat,mixlstr=mixlstr,mixsstr=mixsstr, brnstr=brnstr,$
    xfrac=xfrac,verbose=verbose, _extra=e
;+
;procedure	mcmc_dem_only
;	dedicated subroutine to MCMC_DEM, written only to make MCMC_DEM
;	more readable.  this part decodes the parameters, calls LINEFLX
;	to get predicted fluxes and LIKELI to get the likelihood.
;
;syntax
;	mcmc_dem_only,pars,fobs,fsigma,ulim,line,logT,wvl,Z,scale,ff,prb,$
;	dem,abnd,loopy=loopy,/spliny,nrats=nrats,onlyrat=onlyrat,obsdat=obsdat,$
;	obssig=obssig,idxflx=idxflx,prddat=prddat, _extra=e
;
;parameters
;	pars	[I/O] parameters {DEM(T),ABUND}, ABUND can be updated
;	              via MIXIE() 
;	fobs	[INPUT] observed fluxes
;	fsigma	[INPUT] errors on FOBS
;	ulim	[INPUT] integer array specifying upper limits (1: UL, 0: not)
;	line	[INPUT] emissivities of lines identified with FOBS
;	logT	[INPUT] temperatures at which LINE is defined
;	wvl	[INPUT] wavelengths at which fluxes are observed
;	Z	[INPUT] atomic numbers of each id'd line
;	scale	[INPUT] smoothing scale
;	ff	[OUTPUT] predicted fluxes
;	prb	[OUTPUT] likeihood (could be chi^2)
;	dem	[OUTPUT] DEM
;	abnd	[OUTPUT] abundances
;
;keywords
;	loopy	[INPUT] if set, "smoothing" by converting to T^3/2 DEM
;	spliny	[INPUT] if set, "smooth" by using spline interpolation
;	onlyrat	[INPUT]
;	obsdat	[INPUT]
;	obssig	[INPUT]
;	idxflx	[INPUT]
;	nrats	[INPUT]
;	verbose	[INPUT]
;	_extra	[INPUT] pass defined keywords to LINEFLX (TEMP, EFFAR,
;		WVLAR, KEV, REGRID), LIKELI (SOFTLIM, CHI2, RCHI, BINOM),
;		VARSMOOTH (STEPS, WEIGHT, TYPE), LOOPEM (SLOOP), or
;               ADJUSTIE (TIES) 
;
;restrictions
;	no error checking is done because it is assumed that this is
;	called from and only from MCMC_DEM
;
;history
;	vinay kashyap (May97)
;	changed SIGMA to FSIGMA; added DEM and ABND to parameter list
;	  (VK; AugMM)
;	added keyword LOOPY, and switched out VARSMOOTH for MK_DEM
;	  (VK; Nov01)
;	added keywords ONLYRAT,OBSDAT,OBSSIG,IDXFLX,NRATS and call to
;	  GENERATIO (VK; Dec'02)
;	numerous bug fixes (liwei lin/VK; Dec'02)
;	added keyword PRDDAT (VK; Jan'03)
;	added keyword SPLINY (VK; Jun'03)
;       added keywords BRNSTR,MIXLSTR,MIXSSTR,BRNSTR,XFRAC, as well as 
;         abundance updating via parameter  PARS. PARS now
;         contains updated abundances from MIXIE() on
;         output if MIXIE() keyword ABNDUPDT is set (LL; Jan'04) 
;	bug fix: compute DEMs for /LOOPY and /SPLINY on correct T grid (VK; Jan05)
;	allow limiting number of T components for LOOPY (VK; Mar05)
;       added call to ADJUSTIE to tie ABUND values (LL; Jun05)  
;	added keyword VERBOSE (VK; Jun07)
;-

np=n_params()
if np lt 11 then begin
  print,'Usage: mcmc_dem_only,pars,fobs,fsigma,ulim,line,logT,wvl,Z,scale,$'
  print,'       ff,prb,dem,abnd'
  print,'  computes flux and likelihood for given parameters'
  return
endif

;	deduce
nw=n_elements(wvl) & nt=n_elements(logT)

;	decode parameter array
dem=[pars(0:nt-1)] & abnd=[pars(nt:*)]

;	smooth DEM locally
if keyword_set(loopy) then $
  dem=mk_dem('loop',indem=10.D^(dem),pardem=logT,logT=logT,loopy=loopy,$
  	lcomp=lcomp,verbose=verbose, _extra=e) else $
  if keyword_set(spliny) then $
    dem=mk_dem('spline',logT=logT,indem=dem,pardem=logT,logT=logT,verbose=verbose, _extra=e) else $
      dem=mk_dem('varsmooth',logT=logT,indem=dem,pardem=lscal,verbose=verbose, _extra=e)

;	unlog
if not keyword_set(loopy) then dem=10.D^(dem) & abnd=10.^(abnd)
adjustie, abnd, _extra = e 

if n_tags(mixsstr) gt 0 then begin 
  xfrac=mixie(mixlstr,mixsstr,dem=dem, dlogt = logt, abund=mabnd,obsflx = fobs, $
              ofsig=fsigma,brnstr=brnstr, _extra = e) 
pars[nt:*]=alog10(mabnd) 
;wait, 1 
endif else xfrac=0*wvl+1

;	compute fluxes
if nt eq 1 then begin		;because LINEFLX can't handle LINE(WVL)
  ff=fltarr(nw)
  for iw=0,nw-1 do ff(iw)=lineflx(line(iw),logT,wvl(iw),Z(iw),DEM=DEM,$
	abund=abnd, _extra=e)/xfrac(IW)
endif else ff=lineflx(line,logT,wvl,Z,DEM=DEM,abund=abnd, _extra=e)/xfrac

if keyword_set(nrats) then begin
  generatio,ff,onlyrat,rr,verbose=verbose, _extra=e
  if keyword_set(idxflx) then begin
    if idxflx(0) ne -1 then begin
      prddat=[ff(idxflx),rr]
      ulimnew=[ulim,lonarr(nrats)]
    endif else begin
      prddat=[rr]
      ulimnew=lonarr(nrats)
    endelse
  endif
endif else begin
  obsdat=fobs & prddat=ff & obssig=fsigma & ulimnew=ulim
endelse

;	compute probabilities
;prb=likeli(fobs,ff,dsigma=fsigma,ulim=ulim, _extra=e)
prb=likeli(obsdat,prddat,dsigma=obssig,ulim=ulim, _extra=e)

return
end
