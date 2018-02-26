function logit,x,sigx,sigl=sigl, _extra=e
;+
;function	logit
;	compute and return the logit function.
;
;	LOGIT() maps a variable in the [0,1] range to the
;	real number line [-\infty,+\infty] by applying the
;	transformation ln(x/(1-x)).
;	conversely, UNLOGIT() maps the real number line to
;	[0,1] by applying the transformation e^y/(1+e^y).
;
;syntax
;	lx=logit(x,sigx,sigl=sigl)
;
;parameters
;	x	[INPUT; required] a variable in the range [0,1]
;		* can be an array
;		* value is set to NaN in output if input falls
;		  outside valid range
;	sigx	[INPUT] standard error on X
;		* if <0 and >-1, then ABS(SIGX) is taken to be the
;		  fractional error
;		* if <-1 and >-100, then ABS(SIGX) is taken to be
;		  the percentage error
;		* if <-100, then 1/ABS(SIGX) is taken to be the
;		  fractional error
;
;keywords
;	sigl	[OUTPUT] standard error on LOGIT(X)
;		* note that SIGL = SIGX/(X*(1-X))
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (3may04; based on conversation with Xiao-li Meng)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & nsx=n_elements(sigx)
if np eq 0 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is undefined'
if ok ne 'ok' then begin
  print,'Usage: lx=logit(x,sigx,sigl=sigl)'
  print,'  Apply logit function, lx=ln(x/(1-x))'
  if np gt 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
sx=x*0.
if nsx gt 0 then begin
  sx[*]=sigx[nsx-1L]
  if nsx lt nx then sx[0L:nsx-1L]=sigx[*] else sx[*]=sigx[0L:nx-1L]
endif
osx=where(sx lt 0,mosx)
for i=0L,mosx-1L do begin
  j=osx[i] & xsig=sx[j]
  if xsig gt -1 and xsig lt 0 then sx[j]=x[j]*abs(xsig)
  if xsig gt -100 and xsig le -1 then sx[j]=x[j]*abs(xsig)/100.
  if xsig le -100 then sx[j]=x[j]/abs(xsig)
endfor

o0=where(x[*] eq 0,mo0)
o1=where(x[*] eq 1,mo1)
oo=where(x[*] lt 0 or x[*] gt 1,moo)
ok=where(x[*] gt 0 and x[*] lt 1,mok)

;	outputs
lp=x*0. & sigl=lp

;	compute logit function
if mo0 gt 0 then lp[o0]=-!values.F_INFINITY
if mo1 gt 0 then lp[o1]=!values.F_INFINITY
if moo gt 0 then lp[oo]=!values.F_NAN
if mok gt 0 then lp[ok]=alog(x[ok]/(1.-x[ok]))
;	and errors..
if mok gt 0 then sigl[ok]=sx[ok]/(x[ok]*(1.-x[ok]))

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,lp
end
