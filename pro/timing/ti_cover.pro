function ti_cover,tstart,tstop,tgrid,tmin=tmin,tmax=tmax,tbin=tbin,$
	verbose=verbose, _extra=e
;+
;function	ti_cover
;	compute and return the good time coverage fraction for each
;	bin of a time grid, such as those used to generate lightcurves
;
;syntax
;	tcover=ti_cover(tstart,tstop,tgrid,tmin=tmin,tmax=tmax,tbin=tbin,$
;	verbose=verbose,/slowok)
;
;parameters
;	tstart	[INPUT; required] array of start times of GTI's
;	tstop	[INPUT; required] array of stop times of GTI's
;	tgrid	[I/O] the bin boundaries of the grid for which the
;		coverage must be calculated
;		* if undefined, or has less then 2 elements, will
;		  be calculated internally using the keywords
;		  TMIN, TMAX, and TBIN, and will be overwritten on output
;		* if legal, TMIN, TMAX, and TBIN will be ignored and
;		  the values will not have been changed on output
;
;keywords
;	tmin	[INPUT; default=min(TSTART)] time binning grid minimum
;	tmax	[INPUT; default=max(TSTOP)] time binning grid maximum
;	tbin	[INPUT; default=TMAX-TMIN] the binsize of a regular grid
;		for binning
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		REBINW: SLOWOK (skip call to FINDEX and use INTERPOL instead)
;
;related
;	TI_CLEAN	- cleans up intervals
;	TI_OR		- merges intervals
;	TI_AND		- intersection of intervals
;	TI_WRITE	- writes intervals to file
;	TI_FILTER()	- filter an array by GTI
;
;subroutines
;	REBINW
;	FINDEX
;
;history
;	vinay kashyap (MMVI.IX)
;-

ok='ok' & np=n_params()
ngti=n_elements(tstart) & ngtii=n_elements(tstop) & ntg=n_elements(tgrid)
if np lt 3 then ok='Insufficient parameters' else $
 if ngti eq 0 then ok='TSTART is undefined' else $
  if ngtii eq 0 then ok='TSTOP is undefined' else $
   if ngti ne ngtii then ok='TSTART and TSTOP are incompatible'
if ok ne 'ok' then begin
  print,'Usage: tcover=ti_cover(tstart,tstop,tgrid,tmin=tmin,tmax=tmax,tbin=tbin,$'
  print,'       verbose=verbose)'
  print,'  compute and return the good time coverage fraction for each'
  print,'  bin of a time grid, such as those used to generate lightcurves'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
mint=min(tstart,max=maxt)
ntmin=n_elements(tmin) & if ntmin ne 0 then mint=tmin[0]
ntmax=n_elements(tmax) & if ntmax ne 0 then maxt=tmax[0]
if mint gt maxt then begin
  jnk=mint & mint=maxt & maxt=jnk
  if vv gt 0 then message,$
	'FYI: TMIN and TMAX were inverted',/informational
endif
bint=maxt-mint
if keyword_set(tbin) then bint=abs(tbin[0])
nbin=long((maxt-mint)/bint+0.5)

;	parameters
if ntg lt 2 then tgrid=findgen(nbin+1L)*bint+mint

;	now construct a nominal coverage array
tt=[tstart,tstop]
os=sort(tt) & tt=tt[os] & ntt=n_elements(tt)
gtcov=fltarr(2*ngti-1L)+1
for i=0,ngti-2 do gtcov[2*i+1]=0
if min(tgrid) lt min(tstart) then begin
  tt=[min(tgrid),tt] & gtcov=[0.,gtcov]
endif
if max(tgrid) gt max(tstop) then begin
  tt=[tt,max(tgrid)] & gtcov=[gtcov,0.]
endif
if n_elements(gtcov) eq 1 then begin
  if vv gt 0 then message,'All bins are OK',/informational
  return,0*tgrid[1:*]+1
endif

;	rebin that to the required grid
if vv gt 1000 then stop
gtcovr=rebinw(gtcov,tt,tgrid,slowok=slowok)

;	et voila
return,gtcovr
end
