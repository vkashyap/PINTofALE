function unlogit,lx,sigl,sigx=sigx, _extra=e
;+
;function	unlogit
;	compute and return the inverse of the logit function.
;
;	LOGIT() maps a variable in the [0,1] range to the
;	real number line [-\infty,+\infty] by applying the
;	transformation ln(x/(1-x)).
;	conversely, UNLOGIT() maps the real number line to
;	[0,1] by applying the transformation e^y/(1+e^y).
;
;syntax
;	x=unlogit(lx,sigl,sigx=sigx)
;
;parameters
;	lx	[INPUT; required] a variable in the range [-\infty,+\infty]
;		* can be an array
;		* value is set to NaN in output if input falls
;		  outside valid range
;	sigl	[INPUT] standard error on LX
;		* if <0 and >-1, then ABS(SIGL) is taken to be the
;		  fractional error
;		* if <-1 and >-100, then ABS(SIGL) is taken to be
;		  the percentage error
;		* if <-100, then 1/ABS(SIGL) is taken to be the
;		  fractional error
;
;keywords
;	sigx	[OUTPUT] standard error on UNLOGIT(X)
;		* note that
;		  SIGX = SIGL*exp(LX)/(1+exp(LX))^2 == SIGL*X/(1+exp(LX))
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (6jun05)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(lx) & nsx=n_elements(sigl)
if np eq 0 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='LX is undefined'
if ok ne 'ok' then begin
  print,'Usage: x=unlogit(lx,sigl,sigx=sigx)'
  print,'  Apply inverse logit function, x=exp(lx)/(1+exp(lx))'
  if np gt 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
sl=lx*0.
if nsx gt 0 then begin
  sl[*]=sigl[nsx-1L]
  if nsx lt nx then sl[0L:nsx-1L]=sigl[*] else sl[*]=sigl[0L:nx-1L]
endif
osx=where(sl lt 0,mosx)
for i=0L,mosx-1L do begin
  j=osx[i] & lxsig=sl[j]
  if lxsig gt -1 and lxsig lt 0 then sl[j]=lx[j]*abs(lxsig)
  if lxsig gt -100 and lxsig le -1 then sl[j]=lx[j]*abs(lxsig)/100.
  if lxsig le -100 then sl[j]=lx[j]/abs(lxsig)
endfor

o0=where(finite(lx,/infinity) eq 1 and lx lt 0,mo0)
o1=where(finite(lx,/infinity) eq 1 and lx gt 0,mo1)
oo=where(finite(lx,/nan) eq 1,moo)
ok=where(finite(lx) ne 0,mok)
oko=where(abs(lx[ok]) lt 69,moko)
ok0=where(lx[ok] lt -69,mok0)
ok1=where(lx[ok] gt 69,mok1)

;	outputs
x=0.*lx & sigx=x

;	compute inverse logit function
if mo0 gt 0 then x[o0]=0
if mo1 gt 0 then x[o1]=1
if moo gt 0 then x[oo]=-1
if mok0 gt 0 then x[ok[ok0]]=0
if mok1 gt 0 then x[ok[ok1]]=1
if moko gt 0 then x[ok[oko]]=exp(lx[ok[oko]])/(1.+exp(lx[ok[oko]]))
;	and errors..
if moko gt 0 then sigx[ok[oko]]=sl[ok[oko]]*x[ok[oko]]/(1.+exp(lx[ok[oko]]))

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,x
end
