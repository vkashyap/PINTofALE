function mk_sine,x,ampl,freq,phase,pder,normflx=normflx,$
	missing=missing, _extra=e
;+
;function	mk_sinusoid
;	returns a sinusoid, A*sin(f*X+p)
;
;syntax
;	s=mk_sine(x,ampl,freq,phase,pder,/normflx,$
;	missing=missing)
;
;parameters
;	x	[INPUT array; required] where G(X) must be computed
;	ampl	[INPUT; default: 1] amplitude of the sinusoid
;	freq	[INPUT; default: 2*!pi/range(X)] frequency of sinusoid
;	phase	[INPUT; default: 0] phase of sinusoid [deg]
;	pder	[OUTPUT; optional] partial derivatives of model wrt
;		parameters at each X; calculated only if 5 parameters
;		are supplied in call.
;		* array of size [N(X),3], with columns containing the partial
;		  derivatives wrt AMPL, FREQ, and PHASE respectively
;
;keywords
;	normflx	[INPUT] not used yet
;	missing	[INPUT] 3 element array to populate missing values of
;		AMPL, FREQ, and OFFSET
;	_extra	[JUNK] here only to prevent crashing
;
;description
;	The sinusoid is
;		S(X) = AMPL * sin(FREQ*X+PHASE)
;
;usage summary
;	* call as a function
;	* generates sinusoid values only at specified points X
;	* needs AMPL, FREQ, and PHASE for complete specification
;
;examples
;
;subroutines
;	NONE
;
;history
;	vinay kashyap (Sep04)
;-

np=n_params()
if np lt 1 then begin
  print, 'Usage: s=mk_sine(x,ampl,freq,phase,pder,missing=m,/normflx)'
  print, '  generates a sinusoid s(x)=ampl*sin(freq*X+phase)'
  return,[-1L]
endif

;initialize
nx=n_elements(x) & x0=x[nx/2] & mxx=max(x,min=mnx)

;	figure out the defaults
if not keyword_set(missing) then missing=[1.,2*!pi/(mxx-mnx),0.]
if np lt 3 then faz=missing[2] else faz=phase[0]
if np lt 2 then nu=missing[1] else nu=freq[0]
if np lt 1 then aa=missing[0] else aa=ampl[0]
faz=faz*!pi/180.	;convert from [deg] to [rad]

;	compute function
ss=aa*sin(nu*x+faz)

;	compute partial derivatives
if np ge 5 then begin
  pder = fltarr(nx,3)

  ;z=(x-m+0.0)/c & zz=(1.+z^2)
  ;	partial wrt AMPL
  ssa=sin(nu*x+faz)
  pder[*,0] = ssa[*]
  ;	partial wrt FREQ
  ssf=aa*cos(nu*x+faz)*x
  pder[*,1] = ssf[*]
  ;	partial wrt PHASE
  ssp=aa*cos(nu*x+faz)
  pder[*,2]=ssp[*]
endif

return,ss
end
