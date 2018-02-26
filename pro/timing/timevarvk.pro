;+
;function	timevarvk
;	computes the probability that the given set of photon arrival
;	times are obtainable as a random instance from a constant source.
;
;	smaller values imply higher significances of variability.
;
;	uses a form of Kolmogorov-Smirnoff test, as in TIMEVARKS(), but
;	differs in an important feature: the deviations between the data
;	cdf and the model cdf are _averaged_ over a predefined number of
;	consecutive photons.  the reasoning behind this is that the K-S
;	statistic is a point estimate, whereas a real variability feature
;	always shows some width, i.e., the intensities in the adjacent bins
;	in a lightcurve are correlated.  the averaging of the distance
;	differences will suppress the random Poisson deviations but not
;	the actual signal.
;
;syntax
;	vsigni=timevarvk(times,dobs,pobs,tstart=tstart,tstop=tstop,$
;	nsim=nsim,evtavg=evtavg,verbose=verbose)
;
;parameters
;	times	[INPUT; required] photon arrival times
;		* WARNING: this routine does not handle dead-time
;		  intervals, primbsching
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
;	evtavg	[INPUT] the number of events on either side of the
;		current event over which the distance deviations
;		must be averaged over; i.e., a boxcar smoothing scale
;		of 2*EVTAVG+1 is used
;		* if not given, assumed to be 1
;		* set explicitly to 0 to avoid doing the averaging;
;		  this should produce results identical to that of
;		  TIMEVARKS()
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	FLAT_LC_CDF (included)
;	PROB_KS (from IDLAstro)
;	TI_CLEAN (in $PINTofALE/pro/timing/ti_clean.pro)
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
;	vinay kashyap (MMVI.X; modified from TIMEVARKS())
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
function timevarvk,times,dobs,pobs,tstart=tstart,tstop=tstop,$
	nsim=nsim,verbose=verbose, _extra=e

;	usage
ok='ok' & np=n_params() & nt=n_elements(times)
if np eq 0 then ok='Insufficient parameters' else $
 if nt eq 0 then ok='TIMES are undefined' else $
  if nt lt 3 then ok='TIMES does not have sufficient elements'
if ok ne 'ok' then begin
  print,'Usage: vsigni=timevarvk(times,dobs,pobs,tstart=tstart,tstop=tstop,$
  print,'       nsim=nsim,evtavg=evtavg,verbose=verbose)'
  print,'  compute the significance that a given set of events'
  print,'  come from a varying source'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
numsim=1000L & if keyword_set(nsim) then numsim=long(nsim[0])>1
phavg=1 & if n_elements(evtavg) ne 0 then phavg=long(evtavg[0])
sclavg=2*phavg+1

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
;ksone,tt,'flat_lc_cdf',dobs,pobs,xmin=tmin,xmax=tmax
f0=findgen(nt)/nt & fn=(findgen(nt)+1.)/nt
ff=call_function('flat_lc_cdf',tt,xmin=tmin,xmax=tmax)
d0=(f0-ff) & dn=(fn-ff) & dobs=max(abs(d0)) > (max(abs(dn)))
prob_ks,dobs,nt,pobs
if sclavg gt 1 then begin
  sd0=smooth(d0,sclavg,/edge_truncate)
  sdn=smooth(dn,sclavg,/edge_truncate)
endif else begin
  sd0=d0 & sdn=dn
endelse
sdobs=max(abs(sd0)) > (max(abs(sdn)))

;	now run KS test for realizations of constant lightcurve
;	vs constant model and check what fraction of the d lies
;	above the observed value
sdsim1=fltarr(numsim)
for isim=0L,numsim-1L do begin
  if vv gt 0 and isim eq 100*long(isim/100) then $
	print,strtrim(isim,2)+'..',form='(a,$)'
  r=randomu(seed,nt) & os=sort(r) & r=r[os] & ttrand=interpol(xcft,cft,r)
  ;ksone,ttrand,'flat_lc_cdf',d1,p1,xmin=tmin,xmax=tmax
  f00=findgen(nt)/nt & fn0=(findgen(nt)+1.)/nt
  ff0=call_function('flat_lc_cdf',ttrand,xmin=tmin,xmax=tmax)
  d00=f00-ff0 & dn0=fn0-ff0
  if sclavg gt 1 then begin
    sd00=smooth(d00,sclavg,/edge_truncate)
    sdn0=smooth(dn0,sclavg,/edge_truncate)
  endif else begin
    sd00=d00 & sdn0=dn0
  endelse
  sdsim1[isim]=max(abs(sd00)) > (max(abs(sdn0)))
endfor

ok=where(sdsim1 ge sdobs,mok)
vsigni=float(mok)/float(numsim)
if vv gt 2 then begin
  print,''
  print,'probability of obtaining observed data as a random fluctuation'
  print,'from a constant lightcurve = ',strtrim(vsigni,2)
  if vv gt 10 then begin
    plot,sdsim1,psym=1,xtitle='SIM #',ytitle='smoothed K-S deviation'
    oplot,sdobs+0*sdsim1
    if vv gt 100 then stop,'HALTING; type .CON to continue'
  endif
endif

return,vsigni
end
