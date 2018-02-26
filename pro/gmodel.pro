pro gmodel,x,a,f,pder
;+
;procedure	gmodel
;		procedure written to be compatible with IDL's CURVEFIT
;		GHRS-IDL's WFIT, etc.  computes F(X;A) and partial
;		derivatives as needed.
;
;syntax
;	gmodel,x,a,f,pder
;
;parameters	x	[INPUT; required] points at which to generate models
;		a	[INPUT; required] array of parameter values
;		f	[OUTPUT; required] output f=f(x;a)
;		pder	[OUTPUT; optional] partial derivatives for each
;			parameter
;
;requires
;	MK_3MODEL.PRO
;	MK_GAUSS.PRO
;
;history
;	vinay kashyap (Nov96)
;-

np=n_params(0)
if np lt 3 then begin
  print,'Usage: gmodel,x,a,f,pder'
  print,'  evaluate F and partial(F)/partial(A) at X for parameters A'
  return
endif

forward_function mk_3model,mk_gauss

;	the first element of A says how many components there are!
nc=long(a(0)) & ic=lindgen(nc)

;	the number of elements in A says if all parameters are available
n=n_elements(a)
if (n-1)/nc lt 3 then begin
  message,'No parameters supplied!',/info & stop ;& return
endif

;	go ahead and deconstruct A
p=a(5*ic+1)			;positions
w=a(5*ic+1+1)			;widths
h=a(5*ic+2+1)			;heights
if (n-1)/nc ge 5 then begin
  group=a(5*ic+3+1)		;grouping index
  delp=a(5*ic+4+1)		;delta_x's
endif else begin
  group=indgen(nc) & delp=fltarr(nc)
endelse

;	anything else appended to the end of A can be used as 
;	"keyword" specifier.  right now, none.
k=n-5*nc-1
if k gt 0 then begin
  print,'too many parameters?'
endif

if nc ge 4 then begin
  f=mk_3model(x,p,w,h,group,delp,pder,type='gauss')
endif else begin
  f=mk_3model(x,p,w,h,group,delp,type='gauss')
endelse

return
end
