function mk_efold,x,norm,scale,yoff,pder,verbose=verbose, _extra=e
;+
;function	mk_efold
;	returns an exponential, NORM*exp(-x*SCALE)+YOFF
;
;syntax
;	f=mk_efold(x,norm,scale,yoff,pder,verbose=verbose)
;
;parameters
;	x	[INPUT array; required] where F(X) must be computed
;	norm	[INPUT; required] normalization for exponential
;	scale	[INPUT; required] reciprocal of e-folding decay scale
;	yoff	[INPUT; required] offset from Y=0
;		* NORM,SCALE,YOFF must be scalars or single-element
;		  vectors -- extra elements, if present, are silently
;		  ignored
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if supplied in the calling sequence.
;
;keywords
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[JUNK] here only to prevent crashing
;
;example
;	x=findgen(100) & plot,x,mk_efold(x,100.,0.05,20.)
;
;subroutines
;	NONE
;
;history
;	vinay kashyap (Aug08)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) ;& na=n_elements(aa)
nn=n_elements(norm) & ns=n_elements(scale) & n0=n_elements(yoff)
if np lt 4 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined' else $
  if nn eq 0 then ok='NORM is not defined' else $
   if ns eq 0 then ok='SCALE is not defined' else $
    if n0 eq 0 then ok='YOFF is not defined'
if ok ne 'ok' then begin
  print,'Usage: f=mk_efold(x,norm,scale,yoff,pder,verbose=verbose)'
  print,'  returns an exponential, A*exp(-x*B)+C
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	inputs
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
anorm=norm[0]
ascale=scale[0]
ayoff=yoff[0]

if vv gt 20 then message,'Generating exponential',/informational

;	compute function
ff1=(x*ascale < 69) > (-69)
ff2=anorm*exp(-ff1)
ff=ff2+ayoff

if vv gt 50 then plot,x,ff,xtitle='X',ytitle='f(X)',/ylog

;	compute partial derivatives
if np ge 3 then begin
  pder=fltarr(nx,3)
  ;	partial wrt NORM
  pder[*,0]=ff2/anorm
  ;	partial wrt SCALE
  pder[*,1]=-x*ff2
  ;	partial wrt OFFSET
  pder[*,2]=0*x+1.
endif

return,ff
end
