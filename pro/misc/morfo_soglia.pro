function morfo_soglia,img,snr=snr,gamma=gamma,nsigma=nsigma,$
	centrale=centrale,bitlev=bitlev,zeromin=zeromin,$
	mimg=mimg,thresh=thresh,verbose=verbose, _extra=e
;+
;function	morfo_soglia
;	threshold the image based on some cut derived from applying some
;	kind of filter, and return the image, ideally with the interesting
;	features enhanced
;
;syntax
;	fimg=morfo_soglia(img,gamma=gamma,nsigma=nsigma,centrale=centrale,$
;	bitlev=bitlev,/zeromin,mimg=mimg,thresh=thresh,verbose=verbose)
;
;parametes
;	img	[INPUT; required] input image to be filtered
;		* must be 2D
;
;keywords
;	snr	[INPUT] if given, computes a threshold based
;		on a S/N criterion
;		* if set, IMG is _not_ converted to a byte scale
;	gamma	[INPUT] if given and is +ve, calls GMASCL()
;		to transform the image first according to
;		the power-law scaling transformation
;	nsigma	[INPUT] if given, sets the threshold at
;		mean+NSIGMA*stddev
;		* default is 1
;		* if mean+NSIGMA*stddev < 0, this is not applied
;	centrale [INPUT] if set, carries out a local median
;		filtering with a box of size 2*fix(CENTRALE>1)+1
;		* NSIGMA thresholding is ignored unless NSIGMA
;		  is explicitly set
;		* the filtered pixels are zero'd out unless the keyword
;		  ZEROMIN is explicitly set to 0
;	bitlev	[INPUT] if set, carries out a bit-level filtering
;		down to level BITLEV, with 1=most significant and
;		7 or more=least significant
;		* both CENTRALE and NSIGMA are ignored
;	zeromin	[INPUT] if set, sets all pixels that do not pass
;		the threshold to 0; otherwise will be set to the
;		threshold value
;	mimg	[OUTPUT] the median image, constructed if CENTRALE is set
;	thresh	[OUTPUT] the threshold value used, if histogram-threshold
;		is applied
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	requires David Fanning's GMASCL() function
;
;history
;	vinay kashyap (Jul2006)
;	added _EXTRA keyword; bugfix when NSIGMA was very small (VK; May2007)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(img) & szi=size(img)
if np eq 0 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='IMG is undefined' else $
  if szi[0] ne 2 then ok='IMG must be 2D' else $
   if szi[1] le 3 then ok='IMG axis 1 is collapsed' else $
    if szi[2] le 3 then ok='IMG axis 2 is collapsed'
if ok ne 'ok' then begin
  print,'Usage: fimg=morfo_soglia(img,gamma=gamma,nsigma=nsigma,centrale=centrale,$'
  print,'       bitlev=bitlev,/zeromin,mimg=mimg,thresh=thresh,verbose=verbose)'
  print,'  filter the input image according to some threshold'
  if np ne 0 then message,ok,/informational
  return,-1L
endif
AA=img & nx=szi[1] & ny=szi[2]

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	S/N thresholding
if keyword_set(snr) then begin
  AA=img
  snrthr=3. & if snr[0] gt 0 then snrthr=float(snr[0])
  snrimg=sqrt(AA)
  o0=where(snrimg le snrthr,mo0)
  if mo0 gt 0 then begin
    if keyword_set(zeromin) then AA[o0]=0 else $
	AA[o0]=snrthr*sqrt(AA[o0])
  endif
  return,AA
endif

;	gamma scaling?
if keyword_set(gamma) then AA = gmascl(img,gamma=gamma) else AA=bytscl(img)

;	bit level scaling?
if keyword_set(bitlev) then begin
  ib=8-(long(bitlev[0])<8)
  o0=where(AA le 2^ib,mo0)
  if mo0 gt 0 then begin
    if keyword_set(zeromin) then AA[o0]=0 else AA[o0]=2^ib
  endif
  return,AA
endif

;	median filtering?
if keyword_set(centrale) then begin
  mimg=0*img+median(AA)	;a place to store the medians
  csiz=long(centrale[0])>1
  msiz=2*csiz+1
  for i=0L+csiz,nx-1L-csiz do $
	for j=0L+csiz,ny-1L-csiz do $
		mimg[i,j]=median(AA[i-csiz:i+csiz,j-csiz:j+csiz])
  o0=where(AA lt mimg,mo0)
  if mo0 gt 0 then begin
    ;if keyword_set(zeromin) then AA[o0]=0 else AA[o0]=mimg[o0]
    if n_elements(zeromin) gt 0 then begin
      if zeromin[0] eq 0 then AA[o0]=mimg[o0] else AA[o0]=0
    endif else AA[o0]=0
  endif
  if n_elements(nsigma) eq 0 then return,AA
endif

;	histogram filtering
Amin=min(AA,max=Amax)
if Amax eq Amin then begin
  ;	IDL seems to require this and it has no effect on anything anyway
  Amin=float(Amin) & AA=float(AA)
  Amax=Amin+1
endif
dA=(Amax-Amin)/255.
hAA=histogram(AA,min=Amin,max=Amax,binsize=dA) & xAA=findgen(256)*dA
avgAA=total(xAA*hAA)/total(hAA)
sigAA=sqrt(total(xAA^2*hAA)/total(hAA)-avgAA^2)
nsig=0 & if keyword_set(nsigma) then nsig=nsigma[0]
thresh=avgAA+nsig*sigAA
thresh=fix(thresh+0.5)	;make sure it is definitely greater than avgAA!
o0=where(AA lt thresh,mo0)
if mo0 gt 0 then begin
  if keyword_set(zeromin) then AA[o0]=0 else AA[o0]=thresh
endif

return,AA
end
