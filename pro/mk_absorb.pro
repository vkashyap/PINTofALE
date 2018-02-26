function mk_absorb,x,NH,pder, _extra=e
;+
;function	mk_absorb
;	returns a multiplicative transmission factor for ISM absorption
;
;syntax
;	ismcorr=mk_absorb(x,NH,pder, fH2=fH2,He1=He1,HeII=HeII,/Fano,/ikev,$
;	/wam,/bam,/mam,/noHeH,abund=abund,verbose=verbose)
;
;parameters
;	x	[INPUT array; required] where the function must be computed
;	NH	[INPUT; default: 1e18] H column
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameter
;		at each X; calculated only if 3 parameters are supplied in call.
;
;keywords
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		ISMTAU: FH2, HE1, HEII, FANO, IKEV, WAM, BAM, MAM, NOHEH, ABUND
;
;description
;	wrapper to ISMTAU
;
;subroutines
;	ISMTAU
;
;history
;	vinay kashyap (Jul08)
;-

ok='ok' & np=n_params() & nx=n_elements(x) & nnh=n_elements(NH)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is undefined' else $
  if nnh eq 0 then ok='NH is undefined' else $
   if nnh gt 1 then ok='NH must be a scalar'
if ok ne 'ok' then begin
  print,'Usage: ismcorr=mk_absorb(x,NH,pder, fH2=fH2,He1=He1,HeII=HeII,$'
  print,'       /Fano,/ikev,/wam,/bam,/mam,/noHeH,abund=abund,verbose=verbose)'
  print, '  returns a multiplicative transmission factor for ISM absorption'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;compute exponential part of gaussian
fcorr=exp(-ismtau(x,NH=NH[0], _extra=e))

;compute partial derivatives
if np ge 3 then begin
  pder = fltarr(nx,1)
  pder[*,0] = fcorr*alog(fcorr)/NH[0]
endif

return,fcorr
end
