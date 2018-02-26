function mid2bound,m,halfit=halfit,eps=eps,verbose=verbose, _extra=e
;+
;function	mid2bound
;	given mid-bin values, make an intelligent determination of
;	the bin boundaries and return said boundaries.  given input
;	array of N elements, returns array of N+1 elements.
;
;	(actually, it's not all that intelligent.  uses spline interpolation)
;
;syntax
;	dm=mid2bound(m,/halfit,eps=eps,verbose=verbose)
;
;parameters
;	m	[INPUT; required] the mid-bin values from which to
;		expand out
;
;keywords
;	halfit	[INPUT] if set, does the simpleminded thing, i.e.,
;		puts the bin half-way between 2 points, and extrapolates
;		by repeating previous bins
;	eps	[INPUT] a small number, default is 1e-6
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (OctMM)
;	added keyword HALFIT (VK; Jul01)
;	default EPS had not been defined (VK; Jan05)
;-

;	usage
ok='ok' & np=n_params() & nm=n_elements(m)
if np eq 0 then ok='Insufficient parameters' else $
 if nm eq 0 then ok='mid-bin values not defined'
if ok ne 'ok' then begin
  print,'Usage: b=mid2bound(m,/halfit,eps=eps,verbose=verbose)'
  print,'  return bin-boudaries estimated from mid-bin values'
  return,-1L
endif

;	special case of nm=1
if not keyword_set(eps) then eps=1e-6
if nm eq 1 then return,m+eps*[-1,1]

if keyword_set(HALFIT) then begin
  dM=m[1:*]-m
  b=[m[0]-dM[0]*0.5,m+0.5*dM,m[nm-1L]+dM[nM-2L]*0.5]
  return,b
endif

;	bin boundaries are likely directly in between the mid-bin values
x=lindgen(nm) & y2=spl_init(x,m)
ix=lindgen(nm+1L)-0.5 & b=spl_interp(x,m,y2,ix)

return,b
end
