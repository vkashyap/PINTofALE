function mk_poly1d,x,aa,pder,x0=x0,verbose=verbose, _extra=e
;+
;function	mk_poly1d
;	returns a polynomial, \Sum_{i=0}^{N} a_i (x-x0)^i
;
;syntax
;	py=mk_poly1d(x,aa,pder,x0=x0,verbose=verbose)
;
;parameters
;	x	[INPUT array; required] where G(X) must be computed
;	aa	[INPUT; required] array of coefficients for the polynomial
;		* the degree of the polynomial that is returned is N(AA)-1
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if 4 parameters are supplied in call.
;
;keywords
;	x0	[INPUT; default=0] zero-point of function
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[JUNK] here only to prevent crashing
;
;usage summary
;	* call as a function
;
;subroutines
;	NONE
;
;history
;	vinay kashyap (Aug01)
;	converted array notation to IDL 5; summation of powers now kosher;
;	  corrected sign of partial wrt X0 (VK; Apr02)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & na=n_elements(aa)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined' else $
  if na eq 0 then ok='Coefficiants A_i are not defined'
if ok ne 'ok' then begin
  print,'Usage: py=mk_poly1d(x,aa,pder,x0=x0,verbose=verbose)'
  print,'  returns a polynomial, \Sum_{i=0}^{N} a_i (x-x0)^i'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	inputs
xz=0. & if keyword_set(x0) then xz=0.0+x0[0]
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1

if vv gt 20 then message,'Generating polynomial of degree '+$
	strtrim(na-1,2),/info

;	compute function
py=0.*x
if na gt 1 then py=aa[na-1L]*(x-xz)
for i=na-2L,1L,-1L do py = (py + aa[i])*(x-xz) & py=py+aa[0]
;
;for i=0L,na-1L do py=py+aa[i]*(x-xz)^(i)

if vv gt 50 then plot,x,py,xtitle='X',ytitle='P!d'+strtrim(na-1,2)+'!n(X)'

;	compute partial derivatives
if np ge 3 then begin
  pder=fltarr(nx,na+1L)
  ;	partial wrt X0
  Pyx0=0.*x & for i=1L,na-1L do pyx0=pyx0-aa[i]*i*(x-xz)^(i-1L)
  pder[*,0]=Pyx0[*]
  ;	partial wrt A_i
  for i=0L,na-1L do pder[*,i+1L]=(x[*]-xz)^(i)
endif

return,Py
end
