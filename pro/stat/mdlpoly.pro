function mdlpoly,ydat,xdat,ysig=ysig,yfit=yfit,kord=kord,$
	pfit=pfit,psig=psig,kmin=kmin,kmax=kmax,verbose=verbose,$
	_extra=e
;+
;function	mdlpoly
;	compute and return the minimum description length statistic obtained
;	by fitting polynomials of different orders
;
;syntax
;	mdl=mdlpoly(ydat,xdat,ysig=ysig,yfit=yfit,kord=kord,pfit=pfit,psig=psig,$
;	kmin=kmin,kmax=kmax,verbose=verbose)
;
;parameters
;	ydat	[INPUT; required] data to be fit
;		* must have at least 4 points
;	xdat	[INPUT] abscissae for YDAT
;		* if not given, assumeed to be lindgen(n_elements(YDAT))
;		* if size does not match YDAT, then assumed to be not given
;
;keywords
;	ysig	[INPUT] 1-sigma uncertainties on YDAT
;		* if not given, set equal to stddev(YDAT) for all points
;		* if size does not match YDAT, set equal to mean of input for all points
;		* if scalar, set equal to that value for all points
;	yfit	[OUTPUT] best-fit function corresponding to KORD and PFIT
;	kord	[OUTPUT] polynomial order at which minimum MDL is achieved
;	pfit	[OUTPUT] parameters for the best-fit function
;	psig	[OUTPUT] error bars on PFIT
;	kmin	[INPUT] minimum polynomial order to consider
;		* default is 0, corresponding to a constant
;	kmax	[INPUT] maximum polynomial order to consider
;		* default is the smaller of 10 or half the number of data points
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;requires
;	POLY_FIT
;
;about
;	The MDL (minimum description length) statistic is a way to penalize
;	complex models so that the correct model may be chosen.  The example
;	case (the one you get when you .run mdlpoly) starts with a simulated
;	curve generated from a 5th order polynomial, then uses MDL to pick
;	the best polynomial fit.  You can see how it behaves by changing the
;	coefficients of the generating model (a#) and the assigned error (sigma).
;	For more details on MDL, see
;	http://www.scholarpedia.org/article/Minimum_description_length
;	or the MDL tutorial by Peter Grünwald at arxiv.math:0406077
;
;history
;	Vinay Kashyap (2014sep)
;-

;	usage
ok='ok' & np=n_params() & ny=n_elements(ydat) & nx=n_elements(xdat)
if np eq 0 then ok='Insufficient parameters' else $
 if ny eq 0 then ok='YDAT is undefined' else $
  if ny lt 4 then ok='YDAT too small in size to bother'
if ok ne 'ok' then begin
  print,'Usage: mdl=mdlpoly(ydat,xdat,ysig=ysig,yfit=yfit,kord=kord,pfit=pfit,psig=psig,$
  print,'       kmin=kmin,kmax=kmax,verbose=verbose)'
  print,'  compute and return the minimum description length statistic obtained'
  print,'  by fitting polynomials of different orders'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
yy=ydat & xx=findgen(ny)
if nx eq ny then xx=xdat

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
if not keyword_set(kmin) then kmin=0
if not keyword_set(kmax) then kmax=10<(ny/2)
;
nsig=n_elements(ysig) & sig=fltarr(ny)+stddev(yy)
if nsig eq 1 then sig[*]=ysig[0]
if nsig gt 1 then sig[*]=mean(ysig,/nan)
if nsig eq ny then sig=ysig

;	now do polynomial fits and compute MDL for each case and find the minimum
if vv gt 2 and !d.name eq 'X' then begin
  plot,xx,yy,psym=-4
  if vv gt 3 then for i=0L,ny-1L do oplot,xx[i]*[1,1],yy[i]+sig[i]*[-1,1]
endif
for k=kmin,kmax,1 do begin	;{K=KMIN,KMAX,1
  tmp=poly_fit(xx,yy,k,chisq=chisq,measure_errors=sig,yfit=yfitk,sigma=sigma)
  mdl=chisq + (k/2.)*alog(ny) + lngamma(ny+1)-lngamma(k+1)-lngamma(ny-k+1)
  if vv gt 0 then print,'order,MDL=chisq+penalty',k,mdl,chisq,mdl-chisq
  if vv gt 3 and !d.name eq 'X' then begin
    if k gt 0 then oplot,xx,yfitk,color=(150+5*k mod 255) else $
    	oplot,xx,0*xx+tmp[0],line=1
  endif
  if k eq 0 then begin
    mdlmin=mdl
    yfit=yfitk & kord=k & pfit=tmp & psig=sigma
  endif else begin
    if mdl lt mdlmin then begin
      mdlmin=mdl
      yfit=yfitk & kord=k & pfit=tmp & psig=sigma
    endif
  endelse
endfor				;K=KMIN,KMAX,1}

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,mdlmin
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; example

;	let's start with a 5th order polynomial
if not keyword_set(a0) then a0=1.
if not keyword_set(a1) then a1=-a0/10.
if not keyword_set(a2) then a2=-a1/10.
if not keyword_set(a3) then a3=-a2/10.
if not keyword_set(a4) then a4=-a3/10.
if not keyword_set(a5) then a5=-a4/10.
npt=30. & xx=findgen(npt)-npt/2. & yy=a0+a1*xx+a2*xx^2+a3*xx^3+a4*xx^4+a5*xx^5
if not keyword_set(sigma) then sigma=0.5
ydat=randomn(seed,npt)*sigma+yy & xdat=xx
if not keyword_set(verbose) then verbose=10
ysig=0.*yy+sigma

plot,xx,ydat,psym=-4

;	run
mdl=mdlpoly(ydat,xdat,ysig=ysig,yfit=yfit,kord=kord,pfit=pfit,psig=psig,kmin=kmin,kmax=kmax,verbose=verbose)
print,'mdl=mdlpoly(ydat,xdat,ysig=ysig,yfit=yfit,kord=kord,pfit=pfit,psig=psig,kmin=kmin,kmax=kmax,verbose=verbose)'
oplot,xx,yfit,line=2
help,kord,pfit,psig
print,'This is the best-fit polynomial coefficients and errors:'
for i=0,kord do print,strtrim(pfit[0,i],2)+'	+- '+strtrim(psig[i],2)

;	calling sequence
print,'This is the calling sequence:'
jnk=mdlpoly()

end
