pro fluxcurve,time,chan,charf,chnrg,fluxlc,tgrid,flxerr,$
	ctlc=ctlc,adapbin=adapbin,tmin=tmin,tmax=tmax,tbin=tbin,$
	tstart=tstart,tstop=tstop,shift1=shift1,verbose=verbose, _extra=e
;+
;procedure	fluxcurve
;	compute a fluxed light curve from low-spectral-resolution
;	events data
;
;syntax
;	fluxcurve,time,chan,charf,chnrg,fluxlc,tgrid,flxerr,$
;	ctlc=ctlc,/adapbin,tmin=tmin,tmax=tmax,tbin=tbin,$
;	tstart=tstart,tstop=tstop,/shift1,verbose=verbose,$
;	/slowok,snrthr=snrthr
;
;parameters
;	time	[INPUT; required] photon arrival times
;	chan	[INPUT; required] the PI or PHA values of each photon
;	charf	[INPUT; required] the average effective area in each channel
;	chnrg	[INPUT; required] the average photon energy in each channel
;	fluxlc	[OUTPUT; required] the fluxed light curve, in [ergs/s/cm^2]
;	tgrid	[OUTPUT; required] the bin boundaries over which the light
;		curve is computed
;	flxerr	[OUTPUT; optional] the standard errors on the fluxed light
;		curve, computed as stddev(fluxes)/nphotons in each bin
;
;keywords
;	ctlc	[OUTPUT] the counts light curve, in [ct/s]
;	adapbin	[INPUT] if set, rebins adaptively to ensure a
;		minimum number of counts in each bin, by pushing
;		the counts histogram through SMOOTHIE
;		* assumes errors to be sqrt(CTLC*TBIN)
;		* for each group, the fluxes from the individual
;		  bins are averaged
;	tmin	[INPUT] the minimum time value to consider
;		* default: min(TIME) < min(TSTART)
;	tmax	[INPUT] the maximum time value to consider
;		* default: max(TIME) > max(TSTOP)
;	tbin	[INPUT] the bin size
;		* default: TMAX-TMIN
;	tstart	[INPUT] array of start times of GTI's
;	tstop	[INPUT] array of stop times of GTI's
;	shift1	[INPUT] if set, assumes that the channel numbers
;		start from 1, not from 0
;		* see the FIRSTCHAN field in the response matrix
;		  structure read in via RD_OGIP_RMF()
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		SMOOTHIE : SNRTHR
;		REBINW : SLOWOK
;
;subroutines
;	TI_COVER
;	REBINW
;	FINDEX
;	SMOOTHIE
;
;history
;	vinay kashyap (MMVI.IX)
;-

;	usage
ok='ok' & np=n_params()
nt=n_elements(time) & nc=n_elements(chan)
nea=n_elements(charf) & nen=n_elements(chnrg)
if np lt 4 then ok='Insufficient inputs' else $
 if np lt 6 then ok='No place to return the output' else $
  if nt eq 0 then ok='TIME is undefined' else $
   if nc eq 0 then ok='CHAN is undefined' else $
    if nt ne nc then ok='TIME and CHAN are incompatible' else $
     if nea eq 0 then ok='CHARF is undefined' else $
      if nen eq 0 then ok='CHNRG is undefined' else $
       if nea ne nen then ok='CHARF and CHNRG are incompatible'
if ok ne 'ok' then begin
  print,'Usage: fluxcurve,time,chan,charf,chnrg,fluxlc,tgrid,flxerr,$'
  print,'       ctlc=ctlc,/adapbin,tmin=tmin,tmax=tmax,tbin=tbin,$'
  print,'       tstart=tstart,tstop=tstop,/shift1,verbose=verbose,$'
  print,'       /slowok,snrthr=snrthr'
  print,'  compute a fluxed light curve in [ergs/s/cm^2] from low spectral'
  print,'  resolution events data'
  if np ne 0 then message,ok,/informational
  return
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
mint=min(time,max=maxt)
okgti='ok' & ngti=n_elements(tstart) & ngtii=n_elements(tstop)
if ngti eq 0 then ok='TSTART not defined' else $
 if ngtii eq 0 then ok='TSTOP not defined' else $
  if ngti ne ngtii then ok='TSTART and TSTOP are incompatible'
if okgti ne 'ok' then begin
  tstart=[mint] & tstop=[maxt]
endif
;
ntmin=n_elements(tmin) & if ntmin ne 0 then mint=tmin[0] < min(tstart)
ntmax=n_elements(tmax) & if ntmax ne 0 then maxt=tmax[0] > max(tstop)
if mint gt maxt then begin
  jnk=mint & mint=maxt & maxt=jnk
  if vv gt 0 then message,$
  	'FYI: TMIN and TMAX were inverted',/informational
endif
bint=maxt-mint+1.
if keyword_set(tbin) then bint=abs(tbin[0])
nbin=long((maxt-mint)/bint+0.5)
;
if keyword_set(shift1) then firstchan=1 else firstchan=0

;	compute ergs/cm^2 for each photon
h=6.6261760e-27	;[ergs s]
c=2.9979246e+10	;[cm/s]
keV=interpol(chnrg,findgen(nen)+firstchan,chan)	;[keV/ph]
ang=12.3985/keV	;[Ang/ph]
nrg=h*c*1e8/ang	;[ergs/ph]
cm2=interpol(charf,findgen(nen)+firstchan,chan)	;[cm^2 ct/ph]
cecf=nrg/cm2	;[ergs/cm^2/ct]

;	make the counts light curve
ctlc=float(histogram(time,min=mint,max=maxt,binsize=bint,locations=tmid,reverse_indices=ri))
nh=n_elements(ctlc)
tgrid=[tmid,tmid[nh-1L]+bint]	;time grid bin boundaries

if vv gt 4000 then stop

;	make the flux light curve
fluxlc=fltarr(nh) & flxerr=fltarr(nh)
rvec=ri[0L:nh] & drvec=rvec[1:*]-rvec & isee=where(drvec ne 0,misee)
for i=0L,misee-1L do begin
  j=isee[i]
  fluxlc[j]=total(cecf[ri[ri[j]:ri[j+1L]-1L]])	;[ergs/cm^2/bin]
  if ctlc[j] gt 1 then $
    flxerr[j]=stddev(cecf[ri[ri[j]:ri[j+1L]-1L]])/sqrt(ctlc[j])
endfor
oo=where(flxerr eq 0 and fluxlc gt 0,moo)
if moo ne 0 then flxerr[oo]=stddev(cecf[oo])/sqrt(moo)

if vv gt 3000 then stop

;	adaptively bin to improve S/N
if nh ge 2 and keyword_set(adapbin) then begin
  smoothie,ctlc,sctlc,sectlc,igrp,yerr=sqrt(ctlc),verbose=(vv<1), _extra=e
  grp=igrp[uniq(igrp,sort(igrp))] & ngrp=n_elements(grp)
  sfluxlc=fluxlc
  sefluxlc=flxerr
  for i=0L,ngrp-1L do begin
    ok=where(igrp eq grp[i],mok)
    if mok gt 0 then sfluxlc[ok]=mean(fluxlc[ok])
    if mok gt 0 then sefluxlc[ok]=sqrt(total(flxerr[ok]^2)/mok)
  endfor
  ctlc=sctlc
  fluxlc=sfluxlc
  flxerr=sefluxlc
endif

;	correct for GTIs
fcover=fltarr(nh)+1.
if okgti ne 'ok' then begin
  if vv gt 5 then message,'TSTART and TSTOP are not used',/informational
endif else fcover=ti_cover(tstart,tstop,tgrid, _extra=e)

if vv gt 2000 then stop

o0=where(fcover eq 0,mo0) & o1=where(fcover gt 0,mo1)
if mo0 ne 0 then ctlc[o0]=0.
if mo0 ne 0 then fluxlc[o0]=0.
if mo1 ne 0 then ctlc[o1]=ctlc[o1]/float(bint)/fcover[o1]	;[ct/s]
if mo1 ne 0 then fluxlc[o1]=fluxlc[o1]/float(bint)/fcover[o1]	;[ergs/s/cm^2]
if mo1 ne 0 then flxerr[o1]=flxerr[o1]/float(bint)/fcover[o1]	;[ergs/s/cm^2]

if vv gt 1000 then stop

return
end
