pro x3model,x,a,f,pder,rmfstr=rmfstr,x_seg=x_seg,_extra=e
;+
;procedure	x3model
;	procedure written to be compatible with FIT_LEVMAR, and also
;	with IDL's CURVEFIT, GHRS-IDL's WFIT, etc.
;	computes F(X;A) and partial derivatives as needed for a
;	combination of 3-parameter family of curves
;
;syntax
;	x3model,x,a,f,pder, group=group,delp=delp,type=type,/asis,/allcomp,$
;	missing=missing,/fwhm,/norm,betap=betap,vrot=vrot
;
;parameters
;	x	[INPUT; required] points at which to generate models
;	a	[INPUT; required] array of parameter values
;		* N_ELEMENTS(A) must be a multiple of 3
;		* expected to be sequential triples of (position, width, height)
;	f	[OUTPUT; required] output f=f(x;a)
;	pder	[OUTPUT; optional] partial derivatives for each parameter
;
;keywords
;	rmfstr  [INPUT] structure w/ OGIP standard rmf with which to
;               convolve model. SEE: RD_OGIP_RMF.PRO
;	x_seg	[INPUT] x-axis wavelength grid on which to construct
;		the model spectrum; used when RMFSTR is set, and X is
;		then ignored
;       _extra	[INPUT] pass defined keywords to subroutines
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
;	SCRMF
;	REGROUP_RMF
;
;history
;	vinay kashyap (Oct98; modified from GMODEL.PRO)
;	added MK_SLANT (VK; MarMM)
;	added MK_ROGAUSS,MK_POWLAM (VK; JunMM)
;	added keywords RMFSTR and X_SEG and call to SCRMF (LL; Aug03)
;	added MK_SINUSOID (VK; Sep04)
;-

;help, rmfstr_cut, /str
np=n_params(0)
if np lt 3 then begin
  print,'Usage: x3model,x,a,f,pder, type=type,/asis,/allcomp,missing=missing,$'
  print,'       etc.'
  print,'  evaluate F and partial(F)/partial(A) at X for parameters A'
  return
endif

forward_function mk_3model,mk_gauss,mk_lorentz,mk_slant,mk_rogauss,mk_powlam,mk_sinusoid

;	how many components?
nc=n_elements(a)/3 & na=n_elements(a)
if na-3*nc ne 0 then message,'insufficient parameters!'

;	go ahead and deconstruct A
ic=lindgen(nc)
p=a(3*ic)				;positions
w=a(3*ic+1)			;widths
h=a(3*ic+2)			;heights

;	create the model
if keyword_set(x_seg) then xx = x_seg  else xx=x
if np ge 4 then f=mk_3model(xx,p,w,h,pder, _extra=e) else $
  f=mk_3model(xx,p,w,h, _extra=e)
if keyword_set(rmfstr) then begin
  scrmf,(12.398521/xx),f,chan,ff,rmfstr, _extra=e 
  f = ff 
endif

return
end
