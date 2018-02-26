function findscale,curve,dim,crunch=crunch,half=half,pick=pick,choice=choice,$
	eps=eps,_extra=e
;+
;function	findscale
;	returns the lengthscale in pixels at each point on the given curve
;
;syntax
;	ls=findscale(curve,dim,/crunch,/half,pick=pick,choice=choice,eps=eps)
;
;parameters
;	curve	[INPUT; required] regularly gridded array of function
;		values to be used to compute the length scales
;		* if scalar, returns 0
;		* if 2D,
;		  -- use DIM to specify primary dimension
;		  -- compute lengthscales separately along each projection
;		* if >2D, convert to 1D
;	dim	[INPUT; default=1] primary dimension in case of 2D array
;		(e.g., if CURVE=CURVE(NX,NY), DIM=2 returns SCALE=SCALE(NY))
;
;keywords
;	crunch	[INPUT] if set, and CURVE is 2D, collapses the array along
;		the secondary dimension to generate 1D curve
;	half	[INPUT] if set, returns the half-scale
;	pick	[INPUT; default=0] if 2D, specifies how to combine the
;		scales computed at the different cuts
;		0: pick the smallest scale
;		1: pick the largest scale
;		2: get the average
;	choice	[INPUT; default=0] what algorithm to use to find the scale?
;		0: MexicanHat wavelet
;		1: use inverse of 1st derivative
;		2: radius of curvature
;		3: stepped toggle
;	eps	[INPUT; default=1e-7] small number
;	_extra	[JUNK] ignore.  here only to prevent crashing program.
;
;subroutines
;	WVLT_SCALE [ROOFN]
;
;history
;	vinay kashyap (Apr97)
;	added CHOICE option 3 (VK; Feb03)
;-

;	usage
if n_elements(curve) eq 0 then begin
  print,'Usage: ls=findscale(curve,dim,/crunch,/half,pick=pick,choice=choice,eps=eps)'
  print,'  returns length scales at each point along curve'
  return,0L
endif

;	save inputs
f=curve & if keyword_set(dim) then d=fix(dim) else d=1

;	check dimensions
szf=size(f) & nszf=n_elements(szf)
if szf(0) eq 0 then return,[0L]		;scalar -- return 0
if szf(0) gt 2 then begin		;convert to 1D
  f=[temporary(f(*))] & szf=size(f) & nszf=n_elements(szf)
endif
if szf(0) ne 2 then d=1			;only 1 D, see?
nx=szf(1) & if szf(0) eq 1 then ny=1L else ny=szf(2)

;	if primary dimension is the 2nd D, then transpose the matrix
if d eq 2 then begin
  nx=szf(2) & ny=szf(1) & f=transpose(temporary(f))
endif

;	catch trivial errors
if nx lt 3 then return,lonarr(NX)		;scalar masquerading as array

;	collapse to 1D
if szf(0) eq 2 and keyword_set(crunch) then begin
  g=reform(f(*,0))
  for ix=0,nx-1 do g(ix)=total(f(ix,*))
  f=g & ny=1
endif

;	initialize
if not keyword_set(choice) then choice=0 ;how to make the scales
if not keyword_set(pick) then pick=0	;how to combine scales across dimensions
scale=lonarr(nx)+nx & if pick eq 2 then scale(*)=0	;the output
if not keyword_set(eps) then eps=1e-7			;"epsilon"
norm=lonarr(nx)						;for PICK=2

;	get length scale
scl=lonarr(nx)
for iy=0,ny-1 do begin			;{shtep through secondary dimensions
  g=reform(f(*,iy))		;la function
  dg=deriv(g)			;derivative
  d2g=deriv(dg)			;2nd derivative
  gb=intarr(nx)+1		;coverage function
  ok=where(g lt 0.01*max(g),mok) & if mok gt 0 then gb(ok)=0
  if not keyword_set(choice) then begin
    scl=wvlt_scale(g,_extra=e)
  endif else begin
    if choice(0) eq 1 then scl=ceil(abs(g)/(abs(dg)>eps)) else $
     if choice(0) eq 2 then scl=ceil((1.+dg^2)^(1.5)/(abs(d2g)>eps)) else $
      if choice(0) eq 3 then scl=scl+gb else $
       scl=wvlt_scale(g,_extra=e)
  endelse
  if choice(0) ne 3 then begin
    for ix=0,nx-1 do begin	;if 2D, we gotta pick
      s=scale(ix)
      if g(ix) gt eps*max(f) then begin
        if pick eq 0 then scale(ix)=scl(ix) < s
        if pick eq 1 then scale(ix)=scl(ix) > s
        if pick eq 2 then begin
	  scale(ix)=scl(ix) + s
	  norm(ix)=norm(ix)+1L
        endif
      endif
    endfor
  endif else scale=max(scl)-scl+1L
endfor					;IY=0,NY-1}
norm=norm>1 & if pick eq 2 then scale=ceil(float(scale)/float(norm))

;	return the half-scale if asked
if keyword_set(half) then scale=scale/2

return,scale
end
