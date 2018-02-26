function morfo_circolare,image,arc,rzero,rstep,nstep,gray=gray,$
	mopen=mopen,mclose=mclose,$
	readd=readd,verbose=verbose, _extra=e
;+
;function        morfo_circolare
;	returns an image that is morphologically manipulated with a
;	rotated, annular arclike structure element
; 
;	if MOPEN is set then an opening operation is applied
;	if MCLOSE is set then a closing operation is applied
;
;syntax
;       result=morfo_circolare(image,arc,rzero,rstep,nstep,gray=gray,$
;       /mopen,/mclose,verbose=verbose)
;
;parameters
;	image	[INPUT; required] the image you want to smooth
;		* must be a 2-D array
;	arc	[INPUT; required] the angular size of the arc to be used
;		as a structuring element [degree]
;		* the full circle is completed by moving the arc by ARC/2
;		  and combining the results
;	rzero	[INPUT; required] radius of the innermost annulus to consider
;	rstep	[INPUT; required] the annular width of structure elements [pix]
;	nstep	[INPUT] number of annuli to use
;		* if not given, carries on until the outer radius equals
;		  the larger of the image width or height
;
;keywords
;	gray	[INPUT] explicitly set to 0 to carry out the morphological
;		operations in binary rather than grayscale; default is to
;		use grayscale
;       mopen	[INPUT] if set, runs the IDL function morph_open()
;       mclose	[INPUT] if set, runs the IDL function morph_close()
;		* NOTE: if more than one of MOPEN,MCLOSE
;		  are set, the operations will be daisy-chained in that
;		  particular order
;	readd	[INPUT] if set, adds up the output from the morphological
;		operations at one annular arc to the output from another
;		* the default is to record, at any pixel, the largest
;		  response obtained from all the bits
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (july 06; based on morfo_rottangoli)
;	added keyword _EXTRA (VK; may 07)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(image) & szi=size(image)
narc=n_elements(arc) & nrzero=n_elements(rzero) & nrstep=n_elements(rstep)
if np lt 3 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='IMAGE is undefined' else $
  if narc eq 0 then ok='ARC size is undefined' else $
   if nrzero eq 0 then ok='RZERO is undefined' else $
    if nrstep eq 0 then ok='RSTEP is undefined' else $
     if szi[0] ne 2 then ok='IMAGE is not a 2-D array' else $
      if szi[1] lt 4 then ok='IMAGE axis 1 is collapsed' else $
       if szi[2] lt 4 then ok='IMAGE axis 2 is collapsed'
if ok ne 'ok' then begin
  print,'Usage: result=morfo_circolare(image,arc,rzero,rstep,nstep,gray=gray,$'
  print,'      /mopen,/mclose,verbose=verbose)'
  print,'carries out morphological operations with annular arc-like structure elements'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
nx=szi[1] & ny=szi[2]
arclen=abs(arc[0])>1. & arcshift=arclen/2.
numarc=long(360./arcshift+0.5)
radmin=abs(rzero[0])>0.
delrad=abs(rstep[0])>1.
numann=long((nx>ny)/delrad+0.5)
if n_elements(nstep) ne 0 then numann=long(nstep[0])>1

;	keywords
igray=1
if n_elements(gray) gt 0 then igray=gray[0]
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	output
result=0*image[0]+fltarr(nx,ny)	;at least float

AA=image

for iann=0L,numann-1L do begin	;{for each annular element
  rmin=iann*delrad+radmin & rmax=rmin+delrad & pixmax=long(rmax+0.5)
  print,'annulus:',iann,rmin,rmax
  dd=shift(dist(2*pixmax),pixmax,pixmax)
  xx=findgen(2*pixmax)-pixmax & yy=xx
  pp=fltarr(2*pixmax,2*pixmax)
  for i=0L,2*pixmax-1L do for j=0L,2*pixmax-1L do pp[i,j]=atan(xx[i],yy[j])*180./!pi
  pp=(pp+360) mod 360
  karc=numarc & if rmin eq 0 then karc=1
  for iarc=0L,karc-1L do begin	;{for each arc
    arcmin=iarc*arcshift & arcmax=arcmin+arclen
    print,'arc:',iarc,arcmin,arcmax
    structElem = dd ge rmin and dd le rmax and pp ge arcmin and pp le arcmax
    if rmin eq 0 then structElem = dd ge rmin and dd le rmax

    AA=image
    if keyword_set(mopen) then AA=morph_open(AA,structElem,GRAY=iGRAY)
    if keyword_set(mclose) then AA=morph_close(AA,structElem,GRAY=iGRAY)

    if keyword_set(readd) then result=result+AA else $
      result=result>AA
    if vv gt 10 then tvscl,result
    tvscl,structElem,0

  endfor				;IARC=0,NUMARC-1}
endfor				;IANN=0,NUMANN-1}

return,result
end
