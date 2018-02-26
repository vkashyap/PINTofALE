function cellimg,img,hwidth=hwidth,expix=expix,median=median,mean=mean,sum=sum,$
	_extra=e
;+
;function	cellimg
;	compute pixel-wise scan statistics (sum, median, mean, sdev) for an
;	image by sweeping a square cell across it.
;
;syntax
;	cimg=cellimg(img,hwidth=hwidth,/expix,/sum,/median,/mean,/sdev)
;
;parameters
;	img	[INPUT; required] a 2-D array of numbers
;
;keywords
;	hwidth	[INPUT] half-width of the cell, leads to a (2*HWIDTH+1)x(2*HWIDTH+1) box
;		* default is 1, for a 3x3 cell
;	expix	[INPUT] set to exclude the pixel at the center
;		* default is to include it
;	sum	[INPUT] if set, computes the sum of the values within the cell
;	median	[INPUT] if set, computes the median of the values within the cell
;	mean	[INPUT] if set, computes the mean of the values within the cell
;	sdev	[INPUT] if set, computes the stddev of the values within the cell
;		* default is SUM
;		* if multiple flags are set, the hierarchy is the following:
;		  SUM > MEDIAN > MEAN > SDEV
;		  that is, if /MEDIAN and /MEAN are both set, only /MEDIAN is used
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (2013oct)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(img) & szi=size(img)
if np eq 0 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='IMG is not defined' else $
  if ni lt 9 then ok='IMG is not big enough' else $
   if szi[0] ne 2 then ok='IMG must be 2-D array' else $
    if szi[1] eq 1 or szi[2] eq 1 then ok='IMG is not really 2-D' else $
     if szi[n_elements(szi)-2] eq 6 then ok='IMG cannot be a string array' else $
      if szi[n_elements(szi)-2] eq 7 then ok='IMG cannot be complex array'
if ok ne 'ok' then begin
  print,'Usage: cimg=cellimg(img,hwidth=hwidth,/expix,/sum,/median,/mean,/sdev)'
  print,'  compute pixel-wise scan summary properties'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
hw=1 & if keyword_set(hwidth) then hw=abs(long(hwidth[0]))>1
cell=bytarr(2*hw+1,2*hw+1)+1
if keyword_set(expix) then cell[hw,hw]=0
;
isum=1 & imed=0 & iavg=0 & ivar=0
if keyword_set(sdev) then begin & ivar=1 & iavg=0 & imed=0 & isum=0 & endif
if keyword_set(mean) then begin & ivar=0 & iavg=1 & imed=0 & isum=0 & endif
if keyword_set(median) then begin & ivar=0 & iavg=0 & imed=1 & isum=0 & endif
if keyword_set(sum) then begin & ivar=0 & iavg=0 & imed=0 & isum=1 & endif
;
nx=szi[1] & ny=szi[2]

;	output
cimg=0.*img

for i=0L,nx-1L do begin
  i0=(i-hw) > 0
  i1=(i+hw) < (nx-1L)
  di0=(i0-i)+hw & di1=(i1-i)+hw
  for j=0L,ny-1L do begin
    j0=(j-hw) > 0
    j1=(j+hw) < (ny-1L)
    dj0=(j0-j)+hw & dj1=(j1-j)+hw
    subimg = img[i0:i1,j0:j1]
    subcell = cell[di0:di1,dj0:dj1] & ok=where(subcell gt 0,mok)
    if mok gt 0 then begin
      if keyword_set(isum) then cimg[i,j]=total(subimg[ok]) else $
       if keyword_set(imed) then cimg[i,j]=median(subimg[ok]) else $
        if keyword_set(iavg) then cimg[i,j]=mean(subimg[ok]) else $
	 if keyword_set(ivar) then cimg[i,j]=stddev(subimg[ok])
    endif
  endfor
endfor

return,cimg
end
