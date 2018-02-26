;+
;function	timevarks
;	computes the probability that the given set of photon arrival
;	times are obtainable as a random instance from a constant source.
;
;	smaller values imply higher significances of variability.
;
;syntax
;	vsigni=timevarks(times,dobs,pobs,tstart=tstart,tstop=tstop,$
;	nsim=nsim,verbose=verbose)
;
;parameters
;	times	[INPUT; required] photon arrival times
;		* WARNING: this routine does not handle dead-time
;		  intervals, primbsching, or even GTIs
;	dobs	[OUTPUT] the deviation found in the cdfs between
;		data and the flat lightcurve model
;	pobs	[OUTPUT] significance of DOBS, as calculated for
;		the one-sample Kolmogorov-Smirnoff test
;
;keywords
;	tstart	[INPUT] beginning times of each segment of the
;		good time intervals
;	tstop	[INPUT] ending times of each segment of the GTIs
;		* array sizes of TSTART and TSTOP must match, or
;		  else they will be ignored
;		* if in any case TSTOP < TSTART, they will be flipped
;	nsim	[INPUT; default=1000] number of simulations
;		to run to determine the significance
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	FLAT_LC_CDF (included)
;	KSONE (from IDLAstro)
;	PROB_KS (from IDLAstro)
;	TI_CLEAN (in $PINTofALE/pro/timing/ti_clean.pro)
;
;restrictions
;	keyword values are passed on to FLAT_LC_CDF via KSONE, so
;	make sure that KSONE allows that (post Apr 2005)
;
;description
;	first does an one-sample K-S test on the data, compared to
;	a flat lightcurve model, to determine the deviation and its
;	significance.  the observed deviation is then calibrated
;	with one-sample K-S tests using randomly generated events
;	from a flat lightcurve (also compared to a flat lightcurve
;	model) and comparing the distribution of the deviations thus
;	found with the observed deviation.
;
;history
;	vinay kashyap (MMVI.VIII)
;	added keywords TSTART and TSTOP (VK; MMVI.IX)
;	bug: GTI filtering broken, tzero offsets incorrect, and
;	  FLAT_LC_CDF was changing TMIN,TMAX (VK; MMVII.IV)
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function flat_lc_cdf,x,xmin=xmin,xmax=xmax,verbose=verbose, _extra=e
;{
;return the cdf assuming that the inputs are drawn
;from a flat distribution
;}

if n_elements(x) lt 2 then return,-1L

if n_elements(xmin) eq 0 then xmin=min(x)
if n_elements(xmax) eq 0 then xmax=max(x)
if xmin ge xmax then return,-1L

return,(x-xmin)/(xmax-xmin)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function timevarks,times,dobs,pobs,tstart=tstart,tstop=tstop,$
	nsim=nsim,verbose=verbose, _extra=e

;	usage
ok='ok' & np=n_params() & nt=n_elements(times)
if np eq 0 then ok='Insufficient parameters' else $
 if nt eq 0 then ok='TIMES are undefined' else $
  if nt lt 3 then ok='TIMES does not have sufficient elements'
if ok ne 'ok' then begin
  print,'Usage: vsigni=timevarks(times,dobs,pobs,tstart=tstart,tstop=tstop,$
  print,'       nsim=nsim,verbose=verbose)'
  print,'  compute the significance that a given set of events'
  print,'  come from a varying source'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
numsim=1000L & if keyword_set(nsim) then numsim=long(nsim[0])>1

;	if there are GTIs, then collapse the data
ok='ok' & nstart=n_elements(tstart) & nstop=n_elements(tstop)
if nstart eq 0 then ok='TSTART is undefined' else $
 if nstop eq 0 then ok='TSTOP is undefined' else $
  if nstart ne nstop then ok='TSTART does not match TSTOP'
if ok eq 'ok' then begin
  starttim=tstart & stoptim=tstop
  ti_clean,starttim,stoptim
  tzero=min(starttim) & texp=total(stoptim-starttim) & tt=times
  ;	remove the gaps from the GTIs
  mgti=n_elements(starttim)
  for it=mgti-2L,1,-1L do begin
    ot=where(times ge stoptim[it],mot)
    if mot gt 0 then tt[ot]=tt[ot]-(starttim[it+1L]-stoptim[it])
  endfor
  tt=tt-tzero
  tmin=min(starttim)-tzero & tmax=max(stoptim)-tzero
  ot=where(tt ge 0,mot)
  if mot eq 0 then begin
    message,'No events within the GTIs',/informational
    return,-1L
  endif
  tt=tt[ot]
endif else begin
  tzero=min(times) & tt=times-tzero
  tmin=min(tt,max=tmax)
endelse

;	the input
ot=sort(tt) & tt=tt[ot]
cft=[0.,1.] & xcft=[tmin,tmax]

;	first, run KS test for constant model
ksone,tt,'flat_lc_cdf',dobs,pobs,xmin=tmin,xmax=tmax

;	now run KS test for realizations of constant lightcurve
;	vs constant model and check what fraction of the d lies
;	above the observed value
dsim1=fltarr(numsim) ;& dsim2=dsim1
psim1=fltarr(numsim) ;& psim2=psim1
for isim=0L,numsim-1L do begin
  if vv gt 0 and isim eq 100*long(isim/100) then $
	print,strtrim(isim,2)+'..',form='(a,$)'
  r=randomu(seed,nt) & ttrand=interpol(xcft,cft,r)
  ksone,ttrand,'flat_lc_cdf',d1,p1,xmin=tmin,xmax=tmax
  dsim1[isim]=d1 & psim1[isim]=p1
  ;kstwo,tt,ttrand,d2,p2 & dsim2[isim]=d2 & psim2[isim]=p2
endfor

ok=where(dsim1 ge dobs,mok)
vsigni=float(mok)/float(numsim)
if vv gt 2 then begin
  print,''
  print,'probability of obtaining observed data as a random fluctuation'
  print,'from a constant lightcurve = ',strtrim(vsigni,2)
  if vv gt 10 then begin
    plot,dsim1,psym=1,xtitle='SIM #',ytitle='K-S deviation'
    oplot,dobs+0*dsim1
    if vv gt 100 then stop,'HALTING; type .CON to continue'
  endif
endif

return,vsigni
end
