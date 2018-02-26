function morfo_rottangoli,image,xsize,ysize,angle,gray=gray,$
	mgrad=mgrad,mtop=mtop,mopen=mopen,mclose=mclose,mthin=mthin,$
	reuse=reuse,readd=readd,verbose=verbose, _extra=e
;+
;function        morfo_rottangoli
;	an image that is morphologically manipulated with a
;	rotated, rectangular structure element
; 
;	if MGRAD is set then a gradien operation is applied
;	if MTOP is set then a tophat operation is applied
;	if MOPEN is set then an opening operation is applied
;	if MCLOSE is set then a closing operation is applied
;	if MTHIN is set then a thinning operation is applied
;
;syntax
;       result=morfo_rottangoli(image,xsize,ysize,angle,gray=gray,$
;       /mgrad,/mtop,/mopen,/mclose,/mthin,/readd,/reuse,$
;	verbose=verbose)
;
;parameters
;	image	[INPUT; required] the image you want to smooth
;		* must be a 2-D array
;	xsize	[INPUT; required] the width of the structure element in pixels
;	ysize	[INPUT; required] the height of the structure element in pixels
;	angle	[INPUT; default=0] the angle in degrees by which to rotate the
;		rectangular structure element defined by (XSIZE,YSIZE)
;		* if a vector, the output will be the summed images
;		  from the application of the rectangular structure
;		  element at each of the angles
;
;keywords
;	gray	[INPUT] explicitly set to 0 to carry out the morphological
;		operations in binary rather than grayscale; default is to
;		use grayscale
;       mgrad	[INPUT] if set, runs the IDL function morph_grad()
;       mtop	[INPUT] if set, runs the IDL function morph_tophat()
;       mopen	[INPUT] if set, runs the IDL function morph_open()
;       mclose	[INPUT] if set, runs the IDL function morph_close()
;       mthin	[INPUT] if set, runs the IDL function morph_thin()
;		* NOTE: if more than one of MGRAD,MTOP,MOPEN,MCLOSE,MTHIN
;		  are set, the operations will be daisy-chained in that
;		  particular order
;	readd	[INPUT] if set, adds up the output from the morphological
;		operations at one ANGLE to the output from another ANGLE
;		* the default is to record, at any pixel, the largest
;		  response obtained from all the ANGLEs
;	reuse	[INPUT] if set, uses the output from the morphological
;		operations at one ANGLE as the input for the next ANGLE
;		* overrides both default and READD
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	julia sandell (july 06)
;	tightened up and made bulletproof and changed name from morph_rot
;	  to morfo_rottangoli; morphological analysis with rotating rectangles;
;	  added keywords GRAY, READD, REUSE, VERBOSE (Vinay Kashyap; july 06)
;	added _EXTRA keyword (VK; may 07)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(image) & szi=size(image)
nrx=n_elements(xsize) & nry=n_elements(ysize) & nang=n_elements(angle)
if np lt 3 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='IMAGE is undefined' else $
  if nrx eq 0 then ok='XSIZE is undefined' else $
   if nry eq 0 then ok='YSIZE is undefined' else $
    if szi[0] ne 2 then ok='IMAGE is not a 2-D array' else $
     if szi[1] eq 1 then ok='IMAGE axis 1 is collapsed' else $
      if szi[2] eq 1 then ok='IMAGE axis 2 is collapsed'
if ok ne 'ok' then begin
  print,'Usage: result=morfo_rottangoli(image,xsize,ysize,angle,gray=gray,$'
  print,'       /mgrad,/mtop,/mopen,/mclose,/mthin,/readd,/reuse,verbose=verbose)'
  print,'smooths with rotated, rectangular structure element'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
nx=szi[1] & ny=szi[2]
if nang gt 0 then ang=angle else begin
  ang=0. & nang=1L
endelse

;	keywords
igray=1
if n_elements(gray) gt 0 then igray=gray[0]
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	output
result=0*image[0]+fltarr(nx,ny)	;at least float

;    creating the rectangular structure element
xsz=xsize[0] > 2
ysz=ysize[0] > 2
nsz=ysz
if xsz gt ysz then begin
  message,'structure element is wider than it is taller?',/informational
  nsz=xsz
endif
structElem=fltarr(nsz,nsz)
;
x2sz=xsz/2L
y2sz=ysz/2L
n2sz=nsz/2L
imin=n2sz-x2sz > 0
imax=n2sz+x2sz < (n2sz-1L)
for i=imin,imax,1L do structElem[i,*]=1
if keyword_set(mthin) then begin
  missstructElem=1-structElem > 0
endif

;	apply the morphoplogical operations for each rotation
AA = image
for i=0L,nang-1L do begin		;{step through the angles

  ; rotating structure element by angle degrees
  ; for each element of the angle vector the function runs
  if ang[0] ge 0 then begin
    distances=shift(dist(nsz,nsz),n2sz,n2sz)
    rot_structelem=fix(rot(structelem,ang[i]))
    o0=where(distances gt n2sz,mo0)
    if mo0 gt 0 then rot_structelem[o0]=0
    if keyword_set(mthin) then begin
      rot_missstructElem=1-rot_structelem > 0
    endif
  endif
  print,i,ang[i]

  ; executing a smoothing function, specified by keyword, for information
  ; on what each function does, look at help manual

  if not keyword_set(reuse) then AA=image
  if keyword_set(mgrad) then AA=morph_gradient(AA,rot_structelem)
  if keyword_set(mtop) then AA=morph_tophat(AA,rot_structelem,GRAY=iGRAY)
  if keyword_set(mopen) then AA=morph_open(AA,rot_structelem,GRAY=iGRAY)
  if keyword_set(mclose) then AA=morph_close(AA,rot_structelem,GRAY=iGRAY)
  if keyword_set(mthin) then AA=morph_thin(AA,rot_structelem,rot_missstructElem)

  ;
  ; summing up the images generated by a smoothing function
  ;
  if keyword_set(reuse) then result=AA else $
   if keyword_set(readd) then result=result+AA else $
    result=result>AA
  if vv gt 10 then tvscl,result

endfor				;I=0,NANG-1}

return,result
END
