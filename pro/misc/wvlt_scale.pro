function wvlt_scale,curve,coarse=coarse,_extra=e
;+
;function	findscale
;	returns the lengthscale in pixels at each point on the given curve
;	using the MexHat Wavelet Transforms as a guide
;
;syntax
;	ls=wvlt_scale(curve,/coarse)
;
;parameters
;	curve	[INPUT; required] regularly gridded 1D array of function
;		values to be used to compute the length scales
;		* if scalar, returns 0
;		* if >1D, convert to 1D
;
;keywords
;	coarse	[INPUT; default=2] explore scales at specified coarseness
;		e.g.: if COARSE=1, every scale from 1 to N(CURVE)-1 will
;		be tested
;	_extra	[INPUT] junk; ignore.  here only to prevent crashing program.
;
;history
;	vinay kashyap (Apr97)
;	fixed bug when N(CURVE) was a power of 2 (LiWei Lin/VK; Jan'03)
;-

;	usage
if n_elements(curve) eq 0 then begin
  print,'Usage: ls=wvlt_scale(curve,/coarse)'
  print,'  returns length scales at each point along curve'
  return,0L
endif

;	check inputs
f=curve & szf=size(f)
if szf(0) eq 0 then return,[0L]		;scalar -- return 0
if szf(0) gt 1 then begin		;convert to 1D
  f=[temporary(f(*))] & szf=size(f) & nszf=n_elements(szf)
endif
;
if not keyword_set(coarse) then coarse=2 else coarse=fix(coarse)>1

;	initialize
n=szf(1)				;number of elements in array
fmax=abs(f)				;array of maximum values of f
nn=roofn(n,2)				;change to power of 2 (for FFT)
n2=nn*2L				;double that, for padding (FFT)
x=findgen(n2)				;for defining the MexHat
if nn gt n then ff=[f,f(0:nn-n-1)] else ff=f	;fill out the reminder (periodic bc)...
g=[reverse(ff(0:nn/2-1)),ff,reverse(ff(nn/2:*))]	;...re-position...
f=g					;...rename
scale=lonarr(n)+long(n/sqrt(2))		;the output
	;the sqrt(2) factor because the maximum response of a MexHat to
	;a Gaussian of sigma=d occurs at a scale of d/sqrt(2)

;	catch trivial errors
if n eq 1 then return,[0L]		;scalar masquerading as an array

;	convolve with MexHat wavelet at different scales
for i=1,n-1,coarse do begin		;{shtep through scales
  scal=float(i)						;sigma of the MexHat
  g=(x-n2/2)/scal & g=g^2 & mh=(1.-g) & g=g<138
  norm=3./sqrt(2*scal)/!pi^(0.25)
  mh=mh*exp(-g/2) & mh=shift(mh,n2/2)*norm		;the MexHat
  ff=float(fft(fft(f,1)*fft(mh,1),-1))			;convolve
  ff=ff(n2/4:n2/4+n-1)					;unpad
  oo=where(abs(fmax) lt abs(ff),moo)
  if moo gt 0 then begin				;stronger signal?
    scale(oo)=long(scal/sqrt(2)+0.5) & fmax(oo)=abs(ff(oo))
  endif
endfor					;I=0,N-1,COARSE}

return,scale
end
