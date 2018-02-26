function x3model_f,x,a,pder, _extra=e
;+
;function	x3model_f
;	function written to be compatible with MPFIT, and is
;	identical in all respects to X3MODEL, except that this
;	is a function and not a procedure.
;	returns F(X;A) and computes partial derivatives as needed
;	for a combination of 3-parameter family of curves
;
;syntax
;	f=x3model_f(x,a,pder, group=group,delp=delp,type=type,/asis,/allcomp,$
;	missing=missing,/fwhm,/norm,betap=betap,vrot=vrot)
;
;parameters
;	x	[INPUT; required] points at which to generate models
;	a	[INPUT; required] array of parameter values
;		* N_ELEMENTS(A) must be a multiple of 3
;		* expected to be sequential triples of (position, width, height)
;	pder	[OUTPUT; optional] partial derivatives for each parameter
;
;keywords
;	_extra	[INPUT] pass defined keywords to subroutines
;		MK_3MODEL: GROUP, DELP, TYPE, ASIS, ALLCOMP
;		MK_GAUSS: MISSING, FWHM, NORM
;		MK_LORENTZ: BETAP, MISSING, NORM
;		MK_SLANT: ANGLE, BETAP, MISSING, NORM
;		MK_ROGAUSS: VROT, FWHM, NORM, MISSING
;		MK_SINUSOID: PHASE, MISSING
;
;requires
;	MK_3MODEL
;	MK_GAUSS
;	MK_LORENTZ
;	MK_SLANT
;	MK_ROGAUSS
;	MK_POWLAM
;	MK_SINUSOID
;
;history
;	converted from X3MODEL (VK/LL; Oct02)
;	added MK_SINUSOID (VK; Sep04)
;-

np=n_params(0)
if np lt 2 then begin
  print,'Usage: f=x3model_f(x,a,pder, type=type,/asis,/allcomp,missing=missing,$'
  print,'       etc.)'
  print,'  evaluate F and partial(F)/partial(A) at X for parameters A'
  return,-1L
endif

forward_function mk_3model,mk_gauss,mk_lorentz,mk_slant,mk_rogauss,mk_powlam

;	how many components?
nc=n_elements(a)/3 & na=n_elements(a)
if na-3*nc ne 0 then message,'insufficient parameters!'

;	go ahead and deconstruct A
ic=lindgen(nc)
p=a(3*ic)				;positions
w=a(3*ic+1)			;widths
h=a(3*ic+2)			;heights

;	create the model
if np ge 3 then f=mk_3model(x,p,w,h,pder, _extra=e) else $
  f=mk_3model(x,p,w,h, _extra=e)

return,f
end
