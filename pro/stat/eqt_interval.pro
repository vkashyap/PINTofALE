function eqt_interval,f,x,fsample=fsample,clev=clev,pdfnorm=pdfnorm,$
	xaround=xaround,verbose=verbose, _extra=e
;+
;function	eqt_interval
;	computes and returns the double-sided equal-tail interval
;	[lower_bound,upper_bound] at a specified confidence level.
;	interval includes the central mass of the probability
;	distribution unless explicitly required otherwise.
;
;syntax
;	hpd=eqt_interval(f,x,/fsample,clev=clev,pdfnorm=pdfnorm,$
;	xaround=xaround,verbose=verbose)
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
;	xaround	[INPUT] by default, the central part of the distribution
;		is used, i.e., an equal area is left out on both ends.  if
;		this keyword is defined, a fraction CLEV/2 of the area is
;		used on either side of XAROUND.
;		* if F is a sample, XAROUND must be in the same units as F
;		* you can return single-sided intervals by setting XAROUND
;		  to min(X) or max(X) (min/max of F if sample)
;		* if not set, XAROUND is set to the median
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;
;example
;	for i=1,20 do print,eqt_interval(randomn(seed,10000L)*i,/fsample)
;
;history
;	vinay kashyap (Apr2006)
;-

;	usage
ok='ok' & np=n_params() & nf=n_elements(f) & nx=n_elements(x)
if np eq 0 then ok='Insufficient parameters' else $
 if nf eq 0 then ok='F is not defined' else $
  if nf lt 2 then ok='array size too small'
if ok ne 'ok' then begin
  print,'Usage: hpd=eqt_interval(f,x,/fsample,clev=clev,pdfnorm=pdfnorm,$'
  print,'       xaround=xaround,verbose=verbose)'
  print,'  compute double-sided equal-tail confidence interval'
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
fmin=min(ff,max=fmax)

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0]) > 1
;
crlev=0.68 & if keyword_set(clev) then crlev=0.0+clev[0]
if crlev lt 0 then crlev=abs(crlev)
if crlev ge 1 and crlev lt 100 then crlev=crlev/100.
if crlev ge 100 then crlev = 1.0D - 1.0D/crlev

;	get the cdf and the median level corresponding to XAROUND
if keyword_set(dens) then begin		;(if prob density
	;sort the absissae
  os=sort(xx)
	;compute cdf, getting rid of the -ves if any
  if fmin lt 0 then cff=total(ff[os]-fmin,/cumul) else $
	cff=total(ff[os],/cumul)
	;reset the integral, if required
  if not arg_present(pdfnorm) then cff=cff/max(cff) else begin
    if keyword_set(pdfnorm) then cff=abs(pdfnorm[0])*cff/max(cff)
  endelse
	;find the median
  cfaround=0.5 & fmed=interpol(xx[os],cff,cfaround)
	;change the central point, if required
  if arg_present(xaround) then begin
    if n_elements(xaround) gt 0 then cfaround=interpol(cff,xx[os],xaround[0])
  endif
endif else begin			;DENS)(samples from pdf
	;sort the samples
  os=sort(ff)
	;get the cdf
  cff=dindgen(nf)/(nf-1.)
	;find the median
  fmed=median(ff) & cfaround=0.5
	;change the central point, if required
  if arg_present(xaround) then begin
    if n_elements(xaround) gt 0 then cfaround=interpol(cff,ff[os],xaround[0])
  endif
endelse					;not DENS)

;	compute interval
	;confidence levels at which bounds are to be estimated
cfp=cfaround+crlev*(1.-cfaround)
cfm=cfaround-crlev*(cfaround)
	;note that CFP-CFM = CRLEV
if cfp gt 1 then message,'BUG: upper confidence level > 1'
if cfm lt 0 then message,'BUG: lower confidence level < 0'
	;pick out the points which lie on either side
om=where(cff le cfaround,mom) & op=where(cff ge cfaround,mop)
if keyword_set(dens) then begin		;(if prob density
	;default is to peg to the ends
  eqtm=min(xx,max=eqtp)
	;interpolate in cdf to get interval
  if mom gt 1 then eqtm=interpol((xx[os])[om],cff[om],cfm)
  if mop gt 1 then eqtp=interpol((xx[os])[op],cff[op],cfp)
  	;and done
  return,[eqtm,eqtp]
endif else begin			;DENS)(samples from pdf
	;default is to peg to the ends
  eqtm=min(ff,max=eqtp)
	;interpolate in cdf to get interval
  if mom gt 1 then eqtm=interpol((ff[os])[om],cff[om],cfm)
  if mop gt 1 then eqtp=interpol((ff[os])[op],cff[op],cfp)
  if vv gt 100 then stop,'HALTing.  type .CON to continue'
	;and done
  return,[eqtm,eqtp]
endelse					;not DENS)

return,!values.F_NAN	;should never get here
end
