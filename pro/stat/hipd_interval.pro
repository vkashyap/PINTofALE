function hipd_interval,f,x,fsample=fsample,clev=clev,pdfnorm=pdfnorm,$
	fmode=fmode,verbose=verbose, _extra=e
;+
;function	hipd_interval
;	computes and returns the interval [lower_bound,upper_bound] at a
;	specified confidence level that includes the highest probability
;	densities.  by definition, this is the shortest possible interval.
;
;syntax
;	hpd=hipd_interval(f,x,/fsample,clev=clev,pdfnorm=pdfnorm,$
;	fmode=fmode,verbose=verbose)
;
;parameters
;	f	[INPUT; required] the array for which the confidence interval
;		must be computed
;		* assumed to be a density function unless FSAMPLE is set
;	x	[INPUT; optional] abscissae values
;		* if not given, and F is a density function, then taken
;		  to be the array indices
;		* ignored if FSAMPLE is set
;		* if N(X).GT.1 but N(X).NE.N(F), X is ignored and
;		  FSAMPLE is set automatically
;
;keywords
;	fsample	[INPUT] if set, assumes that F is a set of samples from a
;		density function, as opposed to being the density function
;		itself
;	clev	[INPUT] confidence level at which to compute the intervals
;		* default is 0.68
;		* if < 0, abs(CLEV) is used
;		* if > 1 and < 100, then assumed to be given as a percentage
;		* if > 100, then 1-1/CLEV is used
;	pdfnorm	[INPUT] if set, forces F to integrate to abs(PDFNORM)
;		* if explicitly set to 0, does not normalize F at all
;		* if not set, normalizes to 1
;		* ignored if FSAMPLE is set
;		* WARNING: do not use this keyword unless you know
;		  what you are doing
;	fmode	[OUTPUT] the mode of the distribution
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		MODALPOINT: EPS
;
;subroutines
;	MODALPOINT
;
;description
;	* if density function, find the cumulative integral around the
;	  mode and interpolate at CLEV to derive the HPD interval
;	* if array of sample values, find all intervals corresponding to
;	  CLEV that contain the mode and pick the smallest of the lot.
;	  this is a method devised by Vinay K and Taeyoung Park during
;	  the course of developing BEHR.
;	* works well only for unimodal distributions (and for multimodal
;	  distributions where all the modes are within the interval),
;	  but that's better than nothing.
;	* note that this interval is not invariant under transformations.
;	  for that, one must use equal-tail intervals, see eqt_interval()
;
;example
;	for i=1,20 do print,hipd_interval(randomn(seed,10000L)*i,/fsample)
;
;history
;	vinay kashyap (Mar2006)
;	bug correction with FSAMPLE cdf (VK; Apr2006)
;	added keyword FMODE (VK; Nov2006)
;	bug correction with F(X) case (VK; Apr2007)
;	now handles NaNs in input (VK; Apr2014)
;-

;	usage
ok='ok' & np=n_params() & nf=n_elements(f) & nx=n_elements(x)
if np eq 0 then ok='Insufficient parameters' else $
 if nf eq 0 then ok='F is not defined' else $
  if nf lt 2 then ok='array size too small'
if ok ne 'ok' then begin
  print,'Usage: hpd=hipd_interval(f,x,/fsample,clev=clev,pdfnorm=pdfnorm,$'
  print,'       fmode=fmode,verbose=verbose)'
  print,'  compute highest probability density interval'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	figure out density function or array
dens=1
if keyword_set(fsample) then dens=0
if nx gt 1 and nx ne nf then dens=0
;
if keyword_set(dens) then begin
  xx=findgen(nf) & if nx eq nf then xx=x
endif
;
ff=f

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0]) > 1
;
crlev=0.68 & if keyword_set(clev) then crlev=0.0+clev[0]
if crlev lt 0 then crlev=abs(crlev)
if crlev ge 1 and crlev lt 100 then crlev=crlev/100.
if crlev ge 100 then crlev = 1.0D - 1.0D/crlev

;	find the mode
fmax=max(ff,imx,min=fmin,/nan)
if keyword_set(dens) then fmode=xx[imx] else $
	fmode=modalpoint(ff,verbose=vv, _extra=e)

;	compute interval
if keyword_set(dens) then begin		;(if prob density
	;sort
  os=sort(xx) & xx=xx[os] & ff=ff[os]
	;get cdf
  cff=total(ff,/cumul) & ii=lindgen(nf)
	;normalize to 1, or not, or to PDFNORM
  if not arg_present(pdfnorm) then cff=cff/max(cff,/nan) else begin
    if keyword_set(pdfnorm) then cff=abs(pdfnorm[0])*cff/max(cff,/nan)
  endelse
	;this is where the mode lies on the cdf
  cfmode=interpol(cff,xx,fmode)
	;given the mode's level, how far back can we go?
  cfmin=cfmode-crlev > 0
	;and how far up can we go?
  cfmax=cfmode+crlev < 1
	;the cdf levels translate to these bin indices
  ixmode=long(interpol(ii,cff,cfmode)) > 0
  ixmin=long(interpol(ii,cff,cfmin)) > 0
	;set up the loop to find the smallest range
  go_on=1 & k=ixmin & cfmax=cff[k]+crlev & drng0=max(xx,/nan)-min(xx,/nan)
  while go_on do begin		;{check every interval and pick the smallest
	;the highest index given the current lower index
    ixmax=long(interpol(ii,cff,cfmin+crlev)) < (nf-1L)
	;and the range that corresponds to these indices
    drng=abs(xx[k]-xx[ixmax])
	;check if current interval is smaller
    if drng lt drng0 then begin
	;update smallest interval
      hpdm=xx[k] & hpdp=xx[ixmax] & drng0=drng
      if vv gt 1000 then print,k,hpdm,hpdp,drng0,cff[k],cff[ixmax],cfmode
    endif				;DRNG<DRNG0)
	;next step
    k=k+1L
    cfmax=cff[k]+crlev
    if cff[k] gt cfmode then go_on=0	;quit if hit the mode
    if cfmax ge 1 then go_on=0		;quit if bumped up to the end
    if vv gt 500 and go_on eq 0 then stop,'HALTing; type .CON to continue'
  endwhile			;GO_ON}
	;and done
  return,[hpdm,hpdp]

  ;	;need the reverse() because we want to start from the peak
  ;os=reverse(sort(ff))
  ;	;subtract fmin to account for cases where ff drops below zero
  ;if fmin lt 0 then cff=total(ff[os]-fmin,/cumul) else $
  ;	cff=total(ff[os],/cumul)
  ;	;normalize to 1, or not, or to PDFNORM
  ;if not arg_present(pdfnorm) then cff=cff/max(cff) else begin
  ;  if keyword_set(pdfnorm) then cff=abs(pdfnorm[0])*cff/max(cff)
  ;endelse
  ;	;sort the indices
  ;xx=xx[os]
  ;	;pick out those indices that fall on either side of mode
  ;om=where(xx le fmode,mom) & op=where(xx ge fmode,mop)
  ;	;in case mode is at extreme, peg the range to mode
  ;hpdm=fmode & hpdp=fmode
  ;	;interpolate on integrated function
  ;if mom gt 1 then hpdm=interpol(xx[om],cff[om],crlev)
  ;if mop gt 1 then hpdp=interpol(xx[op],cff[op],crlev)
  ;	;and done
  ;if vv gt 100 then print,hpdm,hpdp,fmode
  ;if vv gt 500 then stop,'HALTing; type .CON to continue'
  ;return,[hpdm,hpdp]
endif else begin			;DENS)(if array of values
	;first make sure everything is sorted
  os=sort(ff) & ii=lindgen(nf)
	;get cdf
	;no PDFNORM nonsense with samples
  cff=dindgen(nf)/(nf-1.)
	;the mode is at this cumulative level
  cfmode=interpol(cff,ff[os],fmode)
	;given the mode's level, how far back can we go?
  cfmin=cfmode-crlev > 0
	;and how far up can we go?
  cfmax=cfmode+crlev < 1
	;the levels translate to these indices
  ixmode=long(interpol(ii,cff,cfmode)) > 0
  ixmin=long(interpol(ii,cff,cfmin)) > 0
	;set up the loop to find the smallest range
  go_on=1 & k=ixmin & cfmax=cff[k]+crlev & drng0=fmax-fmin
  while go_on do begin		;{check every interval and pick the smallest
	;the highest index given the current lower index
    ixmax=k+long(crlev*nf+0.5) < (nf-1L)
	;and the range that corresponds to these indices
    drng=abs(ff[os[k]]-ff[os[ixmax]])
	;check if current interval is smaller
    if drng le drng0 then begin		;(found smaller interval
	;update smallest interval
      hpdm=ff[os[k]] & hpdp=ff[os[ixmax]] & drng0=drng
      if vv gt 100 then print,k,hpdm,hpdp,drng0,cff[k],cfmode
    endif				;DRNG<DRNG0)
	;next step
    k=k+1L
    cfmax=cff[k]+crlev
    if cff[k] gt cfmode then go_on=0	;quit if hit the mode
    if cfmax ge 1 then go_on=0		;quit if bumped up to the end
    if vv gt 500 and go_on eq 0 then stop,'HALTing; type .CON to continue'
  endwhile			;GO_ON}
	;and done
  return,[hpdm,hpdp]
endelse					;not DENS)

return,!values.F_NAN	;should never get here
end
