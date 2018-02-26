function mk_ident,x,y0,x0,m0,pder,gigo=gigo,verbose=verbose, _extra=e
;+
;function	mk_ident
;	returns the input as the output (hence the name)
;	it is useful as an "identity translator" to make
;	a model out of whatever input array.
;
;syntax
;	gy=mk_ident(x,y0,x0,m0,pder,gigo=gigo,verbose=verbose)
;
;parameters
;	x	[INPUT array; required] the points at which the output is
;		computed.
;	y0	[INPUT] the normalization of the output function
;		* if not given, assumed to be 1.0
;	x0	[INPUT] the zero-point of the output, i.e., how much to shift
;		the input
;		* if not given, assumed to be 0.0
;	m0	[INPUT] a stretching factor
;		* if not given, assumed to be 1.0
;		* if Y0, X0, or M0 are arrays, only the first element is used
;		* NOTE: for X0.NE.0 and M0.NE.1, the new functions are simply
;		  interpolated from GIGO
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if 5 parameters are supplied in call.
;
;keywords
;	gigo	[INPUT] the array which goes into the output
;		* if not given, assumed to be 0.*X+1.
;		* if size does not match X, truncated or filled out with
;		  1st element
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
;	vinay kashyap (Apr02)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x)
if np lt 1 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined'
if ok ne 'ok' then begin
  print,'Usage: gy=mk_ident(x,y0,x0,m0,pder,gigo=gigo,verbose=verbose)'
  print,'  returns the input function as the output'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	inputs
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
ny0=n_elements(y0) & nx0=n_elements(x0) & nm0=n_elements(m0) & nf=n_elements(gigo)
if ny0 eq 0 then yz=1.0 else yz=y0[0]*1.0
if nx0 eq 0 then xz=0.0 else xz=x0[0]+0.0
if nm0 eq 0 then mz=1.0 else mz=m0[0]*1.0
if nf ne nx then begin
  fz=fltarr(nx)+x[0]*0.0
  if nf ne 0 then begin
    fz[*]=gigo[0]
    if nf le nx then fz[0L:nf-1L]=gigo[*] else fz[*]=gigo[0L:nx-1L]
  endif
endif else fz=gigo

;	compute function
gz=fz & xx=(x-xz)*mz
gy=interpol(gz,xx,x)

;	compute partial derivatives
if np ge 5 then begin
  pder=fltarr(nx,3)
  ;	partial wrt Y0
  pder[*,0]=gy
  ;	partial wrt X0
  delx=median(xx[1:*]-xx)/5.
  tmp1=interpol(gz,(x-xz-delx)*mz,x) & tmp2=interpol(gz,(x-xz+delx)*mz,x)
  pder[*,1]=(tmp2-tmp1)/delx/2.
  ;	partial wrt M0
  delm=0.01*mz
  tmp1=interpol(gz,(x-xz)*(mz-delm),x) & tmp2=interpol(gz,(x-xz)*(mz+delm),x)
  pder[*,2]=(tmp2-tmp1)/delm/2.
endif

gy=gy*yz

return,gy
end
