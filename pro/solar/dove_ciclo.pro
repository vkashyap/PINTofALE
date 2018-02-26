function dove_ciclo,inpimg,$
	xbgdat=xbgdat,icell=icell,$
	xrect=xrect,xsize=xsize,ysize=ysize,angle=angle,$
	xrclose=xrclose,$
	xbgsub=xbgsub,jcell=jcell,$
	xthresh=xthresh,thrthr=thrthr,nsigma=nsigma,mimg=mimg,$
	xtclose=xtclose,$
	xblob=xblob,areas=areas,blbthr=blbthr,arclev=arclev,deepimg=deepimg,$
	xroi=xroi,pixroi=pixroi,roithr=roithr,$
	xskel=xskel,rmin=rmin,skthr=skthr,$
	xprune=xprune,minpix=minpix,$
	outx=outx,outy=outy, sunrad=sunrad,pixsiz=pixsiz,$
	offsetx=offsetx,offsety=offsety,suncenx=suncenx,sunceny=sunceny,$
	loopsav=loopsav, verbose=verbose, _extra=e
;+
;function	dove_ciclo
;	this function acts as a wrapper to all the morphological
;	operations that are used to extract pixels that define a loop
;	from an image of the solar corona
;
;	the general workflow is:
;	- [CONVOL] - enhance contrast by subtracting background
;	- [MORFO_ROTTANGOLI] open with rectangular structure element
;	- [MORPH_CLOSE] close with small square kernel
;	- [CONVOL] enhance contrast again by subtracting background
;	- [MORFO_SOGLIA] threshold to selectively enhance structures
;	- [MORPH_CLOSE] close the image again
;	- [MORFO_SEGMENTO] group contiguous pixels into region blobs
;	  and throw out regions deemed too small
;	- [ROI_SELEZIONI] interactively select some blobs for further study
;	- [MORFO_SCHELETRO] create a skeleton for these blobs
;	- [MORFO_POTARE] prune the skeleton to make a cleaner image
;	- [] extract the pixels that form the pruned skeleton and
;	  convert them to heliographic coordinates
;
;	keywords control whether any of the calls to the subroutines
;	should be skipped or repeated.  keywords also set up the input
;	parameters to the subroutines.
;
;syntax
;	oimg=dove_ciclo(inpimg,$
;	/xbgdat,icell=icell,$
;	/xrect,xsize=xsize,ysize=ysize,angle=angle,gray=gray,/readd,/reuse,$
;	/xrclose,$
;	/xbgsub,jcell=jcell,$
;	/xthresh,thrthr=thrthr,nsigma=nsigma,mimg=mimg,gamma=gamma,$
;		centrale=centrale,bitlev=bitlev,/zeromin,$
;	/xtclose,$
;	/xblob,blbthr=blbthr,areas=areas,$
;		subidx=subidx,/hardcut,deepimg=deepimg,$
;	/xroi,pixroi=pixroi,roithr=roithr,jitter=jitter,/hadrian,$
;	/xskel,skthr=skthr,rmin=rmin,$
;	/xprune,minpix=minpix,$
;	outx=outx,outy=outy,sunrad=sunrad,pixsiz=pixsiz,$
;	offsetx=offsetx,offsety=offsety,suncenx=suncenx,sunceny=sunceny,$
;	loopsav=loopsav, verbose=verbose)
;
;parameters
;	inpimg	[INPUT; required] image of the corona
;
;keywords	For keywords with names beginning with "X", when set,
;		the corresponding call to the subroutine is skipped.
;		But if that keyword is set to a negative number -N,
;		the subroutine is called N times in sequence (unless
;		otherwise noted below).
;
;	xbgdat	[INPUT] if set, skips initial background subtraction
;	icell	[INPUT; default=1] defines the size of the cell to subtract
;		background - the "source" is measured in a square of size
;		2*ICELL+1 and the "background" in a surrounding layer of 1pix
;	xrect	[INPUT] if set, skips morphological open with rectangular
;		structure element
;	xsize	[INPUT; default=1] width of the structure element
;	ysize	[INPUT; default=10] height of the structure element
;	angle	[INPUT; default=findgen(36)*5] angle of tilt of the rectangle
;		defined by (XSIZE,YSIZE) in [degrees]
;	xrclose	[INPUT] if set, skips morphological close apres open
;	xbgsub	[INPUT] if set, skips background subtraction of modified image
;	jcell	[INPUT; default=ICELL] defines size of cell to subtract
;		background
;	xthresh	[INPUT] if set, skips image thresholding
;	thrthr	[OUTPUT] value used if histogram-threshold is applied
;	nsigma	[INPUT] if given, sets the threshold at mean+NSIGMA*stddev
;		* default is 1
;		* if mean+NSIGMA*stddev < 0, this is not applied
;	mimg	[OUTPUT] the median threshold image, constructed if
;		CENTRALE is set
;	xtclose	[INPUT] if set, skips morphological close apres thresholding
;	xblob	[INPUT] if set, skips grouping the image into distinct regions
;		* if set to -ve number, gets automatically unset
;	areas	[OUTPUT] number of pixels in each labeled region
;	arclev	[INPUT; default=0.68] fraction of the area of the largest
;		blob as the threshold below which to discard blobs
;	deepimg	[OUTPUT] an image that depicts the depth at which a
;		particular pixel was added to a blob
;	blbthr	[INPUT; default=0] threshold value at which to convert
;		grouped image to bitmap
;	xroi	[INPUT] if set, skips selecting a subset of regions
;		* if set to -ve number, gets automatically unset
;	pixroi	[INPUT] if set, assumed to be a set of pixels around
;		which to pick out the region of interest non-interactively
;	roithr	[INPUT; default=0] threshold value at which to ignore pixels
;	xskel	[INPUT] if set, skips making skeletons of regions
;		* if set to -ve number, gets automatically unset
;	rmin	[INPUT; default=1] radius of circle that acts as
;		structure element
;	skthr	[INPUT; default=0] threshold to filter region
;	xprune	[INPUT] if set, skips pruning skeletons
;		* if set to -ve number, gets automatically unset
;	minpix	[INPUT; default=4] number of pixels above which a given
;		branch must be kept while pruning
;	outx	[OUTPUT] the x-pixel indices of the pruned points
;	outy	[OUTPUT] the y-pixel indices of the pruned points
;		* if SUNRAD is set, then floating point values defining
;		  the locations in heliospheric coordinates are returned
;	sunrad	[INPUT] radius of the Sun as seen from the spacecraft,
;		in [arcsec]
;		* if this is set, then OUTX and OUTY are transformed
;		  to heliospheric coordinate system
;		* if set to 1, assumed to be (RSun/AU)*(180/!pi)*3600
;	pixsiz	[INPUT; default=0.5] pixel size, in [arcsec]
;	offsetx	[INPUT; default=0] offset to be applied to OUTX
;		to bring them in line with SUNCENX
;	offsety	[INPUT; default=0] offset to be applied to OUTY
;		to bring them in line with SUNCENY
;	suncenx	[INPUT; default=0] the X coordinate distance of the
;		center of INPIMG from Sun center in [arcsec]
;	sunceny	[INPUT; default=0] the Y coordinate distance of the
;		center of INPIMG from Sun center in [arcsec]
;	verbose	[INPUT] controls chatter
;	loopsav	[INPUT] if set to a filename, saves all the
;		intermediate arrays in the named IDL save file
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		MORFO_ROTTANGOLI: GRAY, READD, REUSE
;		MORFO_SOGLIA: GAMMA, CENTRALE, BITLEV, ZEROMIN
;		MORFO_SEGMENTO: SUBIDX, HARDCUT
;		ROI_SELEZIONI: JITTER, HADRIAN
;
;subroutines
;	MORFO_ROTTANGOLI
;	MORFO_SOGLIA [GMASCL,ERROR_MESSAGE,PATH_SEP,SCALE_VECTOR,CONVERT_TO_TYPE,FPUFIX]
;	MORFO_SEGMENTO
;	ROI_SELEZIONI
;	MORFO_SCHELETRO
;	MORFO_POTARE
;	INICON
;	KILROY
;
;history
;	vinay kashyap (May2007)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(inpimg) & szi=size(inpimg)
nx=szi[1] & ny=szi[2]
if np eq 0 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='INPIMG is undefined' else $
  if szi[0] ne 2 then ok='INPIMG must be a 2D array' else $
   if szi[n_elements(szi)-2] gt 5 then ok='INPIMG cannot be understood' else $
    if nx eq 1 then ok='INPIMG has a collapsed X-axis' else $
     if ny eq 1 then ok='INPIMG has a collapsed Y-axis'
if ok ne 'ok' then begin
  print,'oimg=dove_ciclo(inpimg, /xbgdat,icell=icell, /xrect,xsize=xsize,$
  print,'     ysize=ysize,angle=angle,gray=gray,/readd,/reuse, /xrclose,$'
  print,'     /xbgsub,jcell=jcell, /xthresh,thrthr=thrthr,nsigma=nsigma,$'
  print,'     gamma=gamma,mimg=mimg,centrale=centrale,bitlev=bitlev,/zeromin,$'
  print,'     /xtclose, /xblob,blbthr=blbthr,areas=areas,subidx=subidx,$'
  print,'     /hardcut,deepimg=deepimg,arclev=arclev, /xroi,pixroi=pixroi,$'
  print,'     roithr=roithr,jitter=jitter,/hadrian, /xskel,skthr=skthr,rmin=rmin,$'
  print,'     /xprune,minpix=minpix, outx=outx,outy=outy,sunrad=sunrad,pixsiz=pixsiz,$'
  print,'     offsetx=offsetx,offsety=offsety,suncenx=suncenx,sunceny=sunceny,$'
  print,'     loopsav=loopsav, verbose=verbose)'
  print,'  wrapper for morphological analysis based loop recognition'
  print,'  and identification subroutines'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
; verbose
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
; to run or not to run
ybgdat=1 & yrect=1 & yrclose=1 & ybgsub=1 & ythresh=1 & ytclose=1
yblob=1 & yarea=1 & yroi=1 & yskel=1 & yprune=1
if keyword_set(xbgdat) then begin & if xbgdat[0] lt 0 then ybgdat=fix(abs(xbgdat)) else ybgdat=0 & endif
if keyword_set(xrect) then begin & if xrect[0] lt 0 then yrect=fix(abs(xrect)) else yrect=0 & endif
if keyword_set(xrclose) then begin & if xrclose[0] lt 0 then yrclose=fix(abs(xrclose)) else yrclose=0 & endif
if keyword_set(xbgsub) then begin & if xbgsub[0] lt 0 then ybgsub=fix(abs(xbgsub)) else ybgsub=0 & endif
if keyword_set(xthresh) then begin & if xthresh[0] lt 0 then ythresh=fix(abs(xthresh)) else ythresh=0 & endif
if keyword_set(xtclose) then begin & if xtclose[0] lt 0 then ytclose=fix(abs(xtclose)) else ytclose=0 & endif
if keyword_set(xblob) then begin & if xblob[0] lt 0 then yblob=1 else yblob=0 & endif
if keyword_set(xroi) then begin & if xroi[0] lt 0 then yroi=1 else yroi=0 & endif
if keyword_set(xskel) then begin & if xskel[0] lt 0 then yskel=1 else yskel=0 & endif
if keyword_set(xprune) then begin & if xprune[0] lt 0 then yprune=fix(abs(xprune)) else yprune=0 & endif
;
; associated keywords
if not keyword_set(icell) then i_cell=1 else i_cell=abs(icell[0])>1
if not keyword_set(xsize) then x_size=1 else x_size=xsize[0]
if not keyword_set(ysize) then y_size=10 else y_size=ysize[0]
if not keyword_set(angle) then theta=findgen(36)*5 else theta=angle
if not keyword_set(jcell) then j_cell=i_cell else j_cell=abs(jcell[0])>1
if not keyword_set(thrthr) then thr_thr=0
if not keyword_set(nsigma) then msigma=1 else msigma=nsigma[0]
if not keyword_set(blbthr) then blb_thr=0
if n_elements(arclev) eq 0 then areaclev=0.68 else areaclev=arclev[0]
if not keyword_set(roithr) then roi_thr=0 else roi_thr=roithr[0]
if not keyword_set(skthr) then sk_thr=0 else sk_thr=skthr[0]
if not keyword_set(rmin) then r_min=1 else r_min=rmin[0]
if not keyword_set(minpix) then min_pix=4 else min_pix=minpix[0]
;
if keyword_set(sunrad) then begin
  if sunrad[0] ne 1 then radsun=sunrad[0] else begin
    inicon,fundae=fundae
    radsun=(fundae.RSun/fundae.AU)*(180./!pi)*3600.
  endelse
endif
if not keyword_set(pixsiz) then sizpix=0.5 else sizpix=pixsiz[0]
if not keyword_set(offsetx) then offsetx=0
if not keyword_set(offsety) then offsety=0
if not keyword_set(suncenx) then suncenx=0
if not keyword_set(sunceny) then sunceny=0
;
istep=0
if vv gt 0 then tvscl,inpimg

;	background subtraction
istep=istep+1
Bimg=inpimg
ib=2*abs(i_cell[0])+1 & bb=fltarr(ib,ib)+1
cc=fltarr(ib+2,ib+2)+1 & cc[1:ib,1:ib]=0
for i=0,ybgdat-1 do begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': CONVOL : background subtraction',/informational
  tmp=convol(1.0*Bimg,bb,/edge_wrap)/total(bb) - $
     convol(1.0*Bimg,cc,/edge_wrap)/total(cc)
  Bimg=tmp
endfor
if vv gt 0 then tvscl,Bimg

;	open with rectangle
istep=istep+1
OBimg=Bimg
for i=0,yrect-1 do begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': MORFO_ROTTANGOLI : morphological open',/informational
  tmp=morfo_rottangoli(OBimg>0,x_size,y_size,theta,/mopen,verbose=vv, _extra=e)
  OBimg=tmp
endfor
if vv gt 0 then tvscl,OBimg

;	close
istep=istep+1
COBimg=OBimg
for i=0,yrclose-1 do begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': MORPH_CLOSE : morphological close',/informational
  tmp=float(morph_close(COBimg,bb,/gray))
  COBimg=tmp
endfor
if vv gt 0 then tvscl,COBimg

;	contrast enhance
istep=istep+1
BCOBimg=COBimg
jb=2*abs(j_cell[0])+1 & bb2=fltarr(jb,jb)+1
cc2=fltarr(jb+2,jb+2)+1 & cc2[1:jb,1:jb]=0
for i=0,ybgsub-1 do begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': CONVOL : enhance contrast',/informational
  tmp=convol(1.0*BCOBimg,bb2,/edge_wrap)/total(bb2) - $
      convol(1.0*BCOBimg,cc2,/edge_wrap)/total(cc2)
  BCOBimg=tmp
endfor
if vv gt 0 then tvscl,BCOBimg

;	threshold
istep=istep+1
TBCOBimg=BCOBimg>0
for i=0,ythresh-1 do begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': MORFO_SOGLIA : thresholding',/informational
  tmp=morfo_soglia(TBCOBimg,thresh=thr_thr,nsigma=msigma, _extra=e)
  TBCOBimg=tmp
endfor
if vv gt 0 then tvscl,TBCOBimg

;	close
istep=istep+1
CTBCOBimg=fix(TBCOBimg) gt min(TBCOBimg)
for i=0,ytclose-1 do begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': MORPH_CLOSE : morphological close',/informational
  tmp=morph_close(CTBCOBimg,bb)
  CTBCOBimg=tmp
endfor
if vv gt 0 then tvscl,CTBCOBimg

;	group pixels into regions
istep=istep+1
GCTBCOBimg=CTBCOBimg
if keyword_set(yblob) then begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': MORFO_SEGMENTO : percolate into blobs',/informational
  GCTBCOBimg=morfo_segmento(CTBCOBimg,areas=areas,thresh=blb_thr,$
  arclev=areaclev[0],verbose=vv, _extra=e)
endif
if vv gt 0 then tvscl,GCTBCOBimg

;	pick out regions
istep=istep+1
RGCTBCOBimg=GCTBCOBimg
if keyword_set(yroi) then begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': ROI_SELEZIONI : select blobs',/informational
  print,'' & print,''
  if n_elements(pixroi) gt 0 then ccroi='Non-' else ccroi=''
  message,ccroi+'Interactively select regions of interest',/informational
  print,'' & print,''
  iroi=roi_selezioni(GCTBCOBimg,pixroi,thresh=roi_thr,verbose=vv, _extra=e)
  tmp=GCTBCOBimg
  if iroi[0] ge 0 then tmp[iroi]=0
  RGCTBCOBimg=GCTBCOBimg-tmp
endif
if vv gt 0 then tvscl,RGCTBCOBimg

;	make skeleton
istep=istep+1
SRGCTBCOBimg=RGCTBCOBimg
if keyword_set(yskel) then begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': MORFO_SCHELETRO : make skeleton',/informational
  SRGCTBCOBimg=morfo_scheletro(RGCTBCOBimg,r_min,thresh=sk_thr)
endif
if vv gt 0 then tvscl,SRGCTBCOBimg

;	prune
istep=istep+1
PSRGCTBCOBimg=SRGCTBCOBimg
if keyword_set(yprune) then begin
  if vv gt 0 then message,'STEP '+strtrim(istep,2)+$
	': MORFO_POTARE : prune skeleton',/informational
  PSRGCTBCOBimg=morfo_potare(SRGCTBCOBimg,min_pix,verbose=vv)
endif
if vv gt 0 then tvscl,PSRGCTBCOBimg

;	pick out the points
opix=where(PSRGCTBCOBimg gt 0,mopix) & outx=-1L & outy=-1L
if mopix gt 0 then begin
  outx = opix mod nx
  outy = opix/nx
endif

;	convert to heliospheric coordinates if possible
if mopix gt 0 and keyword_set(sunrad) then begin
  outx = ((suncenx-sizpix*nx/2) + sizpix*(outx+offsetx))/radsun
  outy = ((sunceny-sizpix*ny/2) + sizpix*(outy+offsety))/radsun
endif

;	stop if asked for
if vv gt 1 then tv,gmascl(inpimg+PSRGCTBCOBimg*max(inpimg)/2.,gamma=0.5)
if vv gt 100 then begin
  print,''
  help,inpimg,Bimg,OBimg,COBimg,BCOBimg,TBCOBimg,CTBCOBimg,GCTBCOBimg,RGCTBCOBimg,SRGCTBCOBimg,PSRGCTBCOBimg,opix,outx,outy
  stop,'HALTing; type .CON to continue'
endif

;	record the intermediate steps if asked for
if keyword_set(loopsav) then begin
  if size(loopsav,/type) eq 7 then begin
    if strpos(strlowcase(loopsav),'none') lt 0 then begin
      save,file=loopsav,$
        inpimg,Bimg,OBimg,COBimg,BCOBimg,$
        TBCOBimg,CTBCOBimg,GCTBCOBimg,$
        RGCTBCOBimg,SRGCTBCOBimg,$
        PSRGCTBCOBimg,$
        opix,outx,outy
    endif
  endif
endif

return,PSRGCTBCOBimg
end
