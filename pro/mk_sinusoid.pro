function mk_sinusoid,x,ampl,freq,offset,pder,phase=phase,normflx=normflx,$
	missing=missing, _extra=e
;+
;function	mk_sinusoid
;	returns a sinusoid, A*sin(f*X+p)+c
;
;syntax
;	s=mk_sinusoid(x,ampl,freq,offset,pder,phase=phase,/normflx,$
;	missing=missing)
;
;parameters
;	x	[INPUT array; required] where G(X) must be computed
;	ampl	[INPUT; default: 1] amplitude of the sinusoid
;	freq	[INPUT; default: 2*!pi/range(X)] frequency of sinusoid
;	offset	[INPUT; default: 1] constant offset for sinusoid
;	pder	[OUTPUT; optional] partial derivatives of model wrt
;		parameters at each X; calculated only if 5 parameters
;		are supplied in call.
;		* array of size [N(X),4], with columns containing the partial
;		  derivatives wrt AMPL, FREQ, OFFSET, and PHASE respectively
;
;keywords
;	phase	[INPUT; default=0] the phase of the sinusoid in [deg]
;	normflx	[INPUT] not used yet
;	missing	[INPUT] 3 element array to populate missing values of
;		AMPL, FREQ, and OFFSET
;	_extra	[JUNK] here only to prevent crashing
;
;description
;	The sinusoid is
;		S(X) = AMPL * sin(FREQ*X+PHASE) + OFFSET
;
;usage summary
;	* call as a function
;	* generates sinusoid values only at specified points X
;	* needs AMPL, FREQ, OFFSET, and PHASE for complete specification
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
  print, 'Usage: s=mk_sinusoid(x,ampl,freq,offset,pder,phase=phase,missing=m,/normflx)'
  print, '  generates a sinusoid s(x)=ampl*sin(freq*X+phase)+offset'
  return,[-1L]
endif

;initialize
nx=n_elements(x) & x0=x[nx/2] & mxx=max(x,min=mnx)

;	figure out the defaults
if not keyword_set(missing) then missing=[1.,2*!pi/(mxx-mnx),1.]
if np lt 4 then coff=missing[2] else coff=offset[0]
if np lt 3 then nu=missing[1] else nu=freq[0]
if np lt 2 then aa=missing[0] else aa=ampl[0]
if not keyword_set(phase) then faz=0. else faz=float(phase[0])
faz=faz*!pi/180.	;convert from [deg] to [rad]

;	compute function
ss=aa*sin(nu*x+faz)+coff

;	compute partial derivatives
if np ge 5 then begin
  pder = fltarr(nx,4)

  ;z=(x-m+0.0)/c & zz=(1.+z^2)
  ;	partial wrt AMPL
  ssa=sin(nu*x+faz)
  pder[*,0] = ssa[*]
  ;	partial wrt FREQ
  ssf=aa*cos(nu*x+faz)*x
  pder[*,1] = ssf[*]
  ;	partial wrt OFFSET
  ssc=0*x+1.
  pder[*,2]=ssc[*]
  ;	partial wrt PHASE
  ssp=aa*cos(nu*x+faz)
  pder[*,3]=ssp[*]
endif

return,ss
end
