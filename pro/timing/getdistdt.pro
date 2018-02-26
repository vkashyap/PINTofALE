function getdistdt,lc,ldt,ldtmin=ldtmin,ldtmax=ldtmax,ldtbin=ldtbin,$
	ihdt=ihdt,frate=frate,nolog=nolog,verbose=verbose, _extra=e
;+
;function	getdistdt
;	compute a theoretical distribution of arrival time differences
;	for a supplied light curve.
;
;description
;	the Poisson probability of finding k events in the interval dt,
;	for rate R,
;		p(k|R,dt) = (R*dt)^k * exp(-R*dt) / Gamma(k+1)
;
;	when we consider a list of photon arrival times, the time differences
;	between each photon constitute the case where exactly one event is
;	observed in that duration.  these time differences must therefore be
;	distributed as
;		p(1|R,dt) = (R*dt) * exp(-R*dt)
;
;	If R is varying, R=R(t_i), the counts at any given time [t_i,t_i+tau] is
;	iid Poisson, so the overall distribution is the sum of the distributions
;	at each time, weighted by the expected number of events
;		Sum_i R_i*tau * p(1|R_i,dt)
;
;	this routine takes count rates assumed to be binned at constant bin
;	sizes, computes a distribution for each of the rates, adds them up
;	with suitable weighting, and returns the result.
;
;warning
;	Remember to multiply the output by the time bin width at which
;	the input light curve is binned, if you want to normalize it to
;	the total number of photons.  If the input light curve is not
;	computed with constant bins, then divide by the total of the
;	rates and multiply by the number of expected events.
;
;syntax
;	hdt=getdistdt(lc,ldt,ldtmin=ldtmin,ldtmax=ldtmax,ldtbin=ldtbin,$
;	ihdt=ihdt,frate=frate,/nolog,verbose=verbose)
;
;parameters
;	lc	[INPUT; required] the light curve for which the differences
;		in photon arrival times needs to be computed
;		* must be given in count rates [counts/sec]
;		* negative values will be silently ignored
;	ldt	[OUTPUT] the log10(deltaT) grid over which the output
;		is calculated
;		* output will be a grid of plain deltaT if NOLOG is set
;		* note that this will have one more element than the output
;		  histogram, will have all bin beginnings and endings
;
;keywords
;	ldtmin	[INPUT] the minimum time difference to consider
;		* default is -5
;		* if NOLOG is set, assumed to be not in log10,
;		  and default is changed to 0.
;	ldtmax	[INPUT] the maximum time difference to consider
;		* default is 2
;		* if NOLOG is set, assumed to be not in log10,
;		  and default is changed to 1e3.
;	ldtbin	[INPUT] the bin size for the output LDT grid
;		* default is 0.01
;		* if NOLOG is set, assumed to be not in log10,
;		  and default is changed to 1.
;	ihdt	[OUTPUT] the distributions computed for individual FRATEs
;	frate	[OUTPUT] the count rates for which the IHDT are computed
;	nolog	[INPUT] if set, computes the arrival time differences
;		in regular sec, not in log10.
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example usage
;	.run getdistdt
;
;history
;	vinay kashyap (Nov2009)
;	clarified documentation (Jan2016)
;-

;	usage
ok='ok' & np=n_params() & nlc=n_elements(lc)
if np eq 0 then ok='Insufficient parameters' else $
 if nlc eq 0 then ok='LC is undefined'
if ok ne 'ok' then begin
  print,'Usage: hdt=getdistdt(lc,ldt,ldtmin=ldtmin,ldtmax=ldtmax,ldtbin=ldtbin,$'
  print,'       ihdt=ihdt,frate=frate,/nolog,verbose=verbose)'
  print,'  return the distribution of arrival time differences for given count'
  print,'  rate light curve'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
minldt=-5. & if keyword_set(nolog) then minldt=0.
if keyword_set(ldtmin) then minldt=ldtmin[0]
;
maxldt=2. & if keyword_set(nolog) then maxldt=1e3
if keyword_set(ldtmax) then maxldt=ldtmax[0]
;
binldt=0.01 & if keyword_set(nolog) then binldt=1.
if keyword_set(ldtbin) then binldt=abs(ldtbin[0])
;
if minldt gt maxldt then begin
  if vv gt 0 then message,'reversing LDTMIN and LDTMAX',/informational
  tmp=maxldt & maxldt=minldt & minldt=tmp
endif
if binldt eq 0 then begin
  if vv gt 0 then message,'cannot have zero bin size; returning',$
  	/informational
  return,-1L
endif
if keyword_set(nolog) then begin
  if minldt le 0 then begin
    if vv gt 0 then message,'cannot deal with -ve time intervals',/informational
    minldt=0.
  endif
  if maxldt le 0 then begin
    if vv gt 0 then message,'cannot deal with -ve time intervals',/informational
    maxldt=binldt
  endif
endif

;	remove 0 or -ve rates from light curve
o0=where(lc gt 0,mo0) & rate=lc
if mo0 ne 0 then rate=lc[o0]
;urate=rate[uniq(rate,sort(rate))]
urate=rate & frate=urate;/total(urate)
nr=n_elements(urate)

;	compute grid
nbin=long((maxldt-minldt)/binldt+0.5)
ldt=findgen(nbin)*binldt+minldt
midt=0.5*(ldt[1:*]+ldt)	;mid bin values in log10
if keyword_set(nolog) then midt=alog10(midt)	;make sure it is in log10

;	compute distribution
hdt=dblarr(nbin-1L)
ihdt=fltarr(nr,nbin-1L)
ulc=lc[uniq(lc,sort(lc))] & nulc=n_elements(ulc)
;for i=0L,nr-1L do begin
;  tmp = alog(urate[i])+midt/alog(10) - urate[i]*10.^(midt)
;  oo = where(tmp gt -69,moo)
;  if moo gt 0 then hdt[oo]=hdt[oo]+exp(tmp[oo])
;endfor
for i=0L,nr-1L do begin
  ihdt[i,*]=exp(alog(urate[i])+midt*alog(10) - urate[i]*10.^(midt))
  hdt = hdt + frate[i]*ihdt[i,*]
endfor
;	NOTE: each of the IHDT are normalized to unity, so to match
;	an observed HDT, they should really be weighted by the counts
;	that correspond to each bin in the LC.  But there is no real
;	way to get the bin width into this program, so just remember
;	to multiply the output HDT by the time bin width with which the
;	input light curve was generated, to get the right weighting.
;	Or divide by total(FRATE) and multiply by expected counts.

if vv gt 1000 then stop,'HALTing; type .CON to continue'

;	normalize to unity
dldt=(ldt[1:*]-ldt)
if not keyword_set(nolog) then dldt=dldt/alog10(exp(1))
hdt=hdt*dldt
for i=0,nr-1 do ihdt[i,*]=ihdt[i,*]*dldt

return,hdt
end
;...............................................................................

;	example usage
hdt=getdistdt()

;	example call
lc=fltarr(2000) & lc[500]=exp(-findgen(1500)/500.)
if not keyword_set(verbose) then verbose=10
if not keyword_set(ldtmin) then ldtmin=-5.0
if not keyword_set(ldtmax) then ldtmax=2.0
if not keyword_set(ldtbin) then ldtbin=0.1

message,'LDTMIN='+strtrim(ldtmin,2),/informational
message,'LDTMAX='+strtrim(ldtmax,2),/informational
message,'LDTBIN='+strtrim(ldtbin,2),/informational

hdt=getdistdt(lc,ldt,ldtmin=ldtmin,ldtmax=ldtmax,ldtbin=ldtbin,$
	ihdt=ihdt,frate=frate,nolog=nolog,verbose=verbose)

plot,ldt,hdt,psym=10,xtitle='log!d10!n(!4d!Xt)',/ylog
peasecolr & loadct,3 & peasecolr
sz=size(ihdt)
for i=0,sz[1]-1 do oplot,ldt,frate[i]*ihdt[i,*],col=(i mod 155)+100
oplot,ldt,hdt,psym=10,thick=2
oplot,ldt,frate#ihdt,thick=2,col=3
print,total(frate),total(hdt)

end
