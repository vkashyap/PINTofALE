function morfo_scheletro,image,rmin,thresh=thresh,verbose=verbose,$
	_extra=e
;+
;function	morfo_scheletro
;	compute and return the skeleton representation of an
;	image A using a circular structuring element B.
;	a skeleton is defined as the union of the sets
;		Erosion(A,k*B) - Opening(Erosion(A,k*B),B)
;	for k=0,..K where K is the largest value of k before
;	the sets become empty.
;
;syntax
;	skel=morfo_scheletro(image,rmin,thresh=thresh,verbose=verbose)
;
;parameters
;	image	[INPUT; required] image whose skeleton must be
;		obtained.
;		* must be a 2-D array
;	rmin	[INPUT] radius of the circle that
;		acts as a structuring element
;		* if not given, the default is to use a 3x3 array
;		* if given, uses DIST() to define the array
;		  -- minimum possible radius is hardcoded as 1pix
;
;keywords
;	thresh	[INPUT; default=0] filter IMAGE such that only
;		values greater than THRESH are considered.
;		* if size matches that of IMAGE, a pixel-by-pixel
;		  filtering is carried out
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Jul2006)
;-

;	usage
ok='ok'
np=n_params() & nA=n_elements(image) & szA=size(image)
if np lt 1 then ok='Insufficient parameters' else $
 if nA eq 0 then ok='IMAGE is undefined' else $
  if szA[0] ne 2 then ok='IMAGE is not a 2-D array' else $
   if szA[1] eq 1 then ok='IMAGE axis 1 has collapsed' else $
    if szA[2] eq 1 then ok='IMAGE axis 2 has collapsed'
if ok ne 'ok' then begin
  print,'skel=morfo_scheletro(image,rmin,thresh=thresh)'
  print,'  compute and return the skeleton of an image'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	input image
AA=image & nx=szA[1] & ny=szA[1]

;	structuring element
BB=fltarr(3,3)+1
if n_elements(rmin) ne 0 then begin
  rr=abs(float(rmin[0]))>1
  rpix=long(rr+0.5)
  dd=shift(dist(2L*rpix+1L),rpix,rpix)
  BB = dd le rpix
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
thr=0.*AA & nth=n_elements(thresh)
if nth ne 0 then begin
  szt=size(thresh)
  if szt[0] eq 2 and szt[1] eq nx and szt[2] eq ny then $
	thr=thresh else thr[*]=thresh[0]
endif
AA = AA gt thr
;o0=where(AA le thr,mo0) & if mo0 ne 0 then AA[o0]=0
;o1=where(AA gt thr,mo1) & if mo1 ne 0 then AA[o1]=1

;	the smallest structuring element is a circle that fits into one pixel 
;		Erosion(A,k*B) - Opening(Erosion(A,k*B),B)
Sk=AA-morph_open(AA,BB)
skel=Sk & Ak=AA & kk=0L
while total(Sk) gt 0 do begin
  Ak=erode(Ak,BB)
  Sk=Ak-morph_open(Ak,BB)
  skel=skel+Sk
  ;
  kk=kk+1L
  if vv gt 0 then kilroy
  if vv gt 4 then print,kk
  tvscl,skel
endwhile

return,skel
end
