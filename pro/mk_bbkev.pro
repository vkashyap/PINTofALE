function mk_bbkev,x,norm,temp,pder,verbose=verbose, _extra=e
;+
;function	mk_bbkev
;	returns the Planck function B(E) [ph/s/cm^2/keV] as a
;	function of E [keV] at a given temperature T.
;
;	In general,
;	B(nu)d(nu) = (R/d)^2*(2*pi/c^2) * d(nu)*nu^2/(exp(h*nu/kT)-1)
;	for E=h*nu,
;	B(E)dE	= (R/d)^2 * (2*pi/c^2/h^3) * dE*E^2/(exp(E/kT)-1)
;		== A * (E^2/(exp(E/kT)-1)) * dE
;
;syntax
;	b=mk_bbkev(x,norm,temp,pder,verbose=verbose)
;
;parameters
;	x	[INPUT array; required] where B(X) must be computed
;	norm	[INPUT; required] normalization for B(X)
;	temp	[INPUT; required] temperature in [K]
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if 4 parameters are supplied in call.
;
;keywords
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing
;
;usage summary
;	* call as a function
;
;subroutines
;	NONE
;
;history
;	vinay kashyap (Aug08; based loosely on 1993 blackbody.pro)
;-

;	usage
ok='ok' & np=n_params()
nx=n_elements(x) & nn=n_elements(norm) & nt=n_elements(temp)
if np lt 3 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined' else $
  if nn eq 0 then ok='NORM is not defined' else $
   if nt eq 0 then ok='TEMPERATURE is not defined' else $
    if nn gt 1 then ok='NORM should not be an array' else $
     if nn ne nt then ok='TEMPERATURE incompatible with NORM'
if ok ne 'ok' then begin
  print,'Usage: B=mk_bbkev(x,norm,temperature,pder,verbose=verbose)'
  print,'  compute and return a Planck spectrum B(E) [ph/s/cm^2/keV]'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	inputs
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
;
aa=norm[0]
tt=temp[0]

if vv gt 20 then begin
  ct=strtrim(tt,2)
  cc='B(keV;T='+ct+')'
  message,'Generating a Planck curve '+cc,/informational
endif

;	compute function
hh=6.6261760d-27	;[erg s]
cc=2.9979246d+10	;[cm/s]
kB=1.3806620e-16	;[erg/K]
ee=1.6021892e-19	;[C]

onekeV=ee*1e7*1e3	;[erg/keV]
Eerg=x*onekeV		;[keV]*[erg/keV]=[erg]

xx=(Eerg/kB/tt) < 69.
tmp=(exp(xx)-1.) > 1d-30
;B(E)dE = (R/d)^2*(2*pi/c^2/h^3) * E^2/(exp(E/kT)-1) * dE
Bx = exp( alog(aa)+2.*alog(x)-alog(tmp) )

if vv gt 50 then plot,x,Bx,xtitle='E [keV]',ytitle='B(E)'

;	compute partial derivatives
if np ge 3 then begin
  pder=fltarr(nx,2)
  ;	partial wrt norm
  pder[*,0]=Bx/aa
  ;	partial wrt T
  tmp2=exp(alog(Bx)+alog(xx)-alog(tmp)+xx)
  pder[*,1]=tmp2
endif

return,Bx
end
