function morfo_segmento,image,areas,thresh=thresh,arclev=arclev,$
	subidx=subidx,hardcut=hardcut,deepimg=deepimg,verbose=verbose,$
	_extra=e
;+
;function	morfo_segmento
;	percolates along adjacent connected pixels in a bitmap
;	until all the pixels that are connected to each other get
;	grouped into contiguous regions and are so labeled; an image
;	containing the region numbers for each pixel is returned
;
;syntax
;	oimg=morfo_segmento(img,areas,thresh=thresh,arclev=arclev,$
;	subidx=subidx,/hardcut,deepimg=deepimg,verbose=verbose)
;
;parameters
;	image	[INPUT; required] 2D image in which to separate
;		out the connected blobs
;		* the blobs are defined as any pixel with intensity>0
;		  unless keyword THRESH is set
;	areas	[OUTPUT] number of pixels in each labeled region
;		* remember the IDL index shift! AREAS[0] is the number
;		  of pixels in the background, and AREAS[SOURCE] gives
;		  the number of pixels that constitute the blob SOURCE
;
;keywords
;	thresh	[INPUT] the threshold value at which to convert
;		the input image into a bitmap
;		* if not set, assumed to be 0
;	arclev	[INPUT] fraction of the maximum area size to keep
;		in the output -- all regions of size smaller than
;		this are zeroed out
;		* default is 0, i.e., no filtering is done
;		* if .le. -1, taken to be an absolute value cut
;		  if .gt. -1 && .lt. 0, absolute value is choosen
;		  if .ge. 1 && .lt. 100, assumed to be percentage
;		  if .ge. 100, assumed to be 1-1/ARCLEV
;	subidx	[INPUT] index array defining a subset of IMAGE
;		to which one must confine attention
;		* if not set, assumed to be LINDGEN(N_ELEMENTS(IMAGE))
;	hardcut	[INPUT] if set, places a hard cut at the boundaries
;		defined by SUBIDX -- the default is to track the region
;		out of SUBIDX if any part of it overlaps
;	deepimg	[OUTPUT] an image that depicts the depth at which a
;		particular pixel was added to a blob
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	this is an updated version of morfo_filtri(), and is much faster
;	because it uses a few algorithmic and IDL tricks to speed things
;	up.  for example, rather than do a linear search going through
;	every pixel once, it does array operations.  each pixel is initially
;	assigned a unique number, and at each iteration, the numbers of
;	each adjacent pixels are compared and the smallest of the numbers
;	is assigned to all of them.  thus, all the pixels that belong to
;	each blob end up having the same number.  these numbers are then
;	recast into a gap-free sequence thereafter.
;	note: there will be no significant improvements in speed in some
;	special cases, such as if there is a very large connected region;
;	this is never *slower* than morfo_filtri() however.
;
;subroutines
;	KILROY
;
;history
;	vinay kashyap (Sep2006; based on morfo_filtri.pro)
;	added keywords ARCLEV and _EXTRA (VK; May2007)
;	was failing when there was only one blob (VK; Jun2007)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(image) & szi=size(image)
if np eq 0 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='IMAGE is undefined' else $
  if szi[1] eq 1 then ok='IMAGE X-axis is collapsed' else $
   if szi[2] eq 1 then ok='IMAGE Y-axis is collapsed'
if ok ne 'ok' then begin
  print,'Usage: oimage=morfo_segmento(image,areas,thresh=thresh,$'
  print,'       subidx=subidx,verbose=verbose)'
  print,'  percolates to find connected blobs of pixels'
  if np ne 0 then message,ok,/informational
  return,-1L
endif
nx=szi[1] & ny=szi[2]

;	keywords
thr=0. & if keyword_set(thresh) then thr=thresh[0]
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
subimg=0*image+1 & nsub=n_elements(subidx)
if nsub gt 0 then begin
  subimg[*]=0 & subimg[subidx]=1
endif
arthr=0.
if keyword_set(arclev) then begin
  tmp=0.0+arclev[0]
  if tmp lt -1 then arthr=abs(tmp) else begin
    if tmp lt 0 then arthr=abs(tmp)
    if tmp ge 0 and tmp lt 1 then arthr=tmp
    if tmp ge 1 and tmp lt 100 then arthr=tmp/100.
    if tmp ge 100 then arthr=(1.D - 1./tmp)
  endelse
endif

;	convert input image to a bitmap
img = image gt thr

;	initialize
idx=where(img ne 0,midx)
if midx eq 0 then begin
  message,'No blobs found in image',/informational
  areas=[midx]
  return,img
endif
;
ix = idx mod nx
iy = idx / nx
idximg=0*long(image) & idximg[idx]=lindgen(midx)+1
;
sidx=where(idximg gt 0 and subimg eq 0,msidx)
if msidx gt 0 then begin
  if keyword_set(hardcut) then begin
    idximg[sidx]=0
    idx=where(idximg ne 0,midx)
    ix = idx mod nx
    iy = idx / nx
    if midx gt 0 then idximg[idx]=lindgen(midx)+1
  endif else begin
    idximg[sidx]=midx+1L
  endelse
endif

if vv gt 3000 then stop,'HALTING; type .CON to continue'

;	percolate and identify all segments
nuncorr=total(idximg) & depth=0L & deepimg=long(idximg gt 0)
while nuncorr gt 0 do begin
  tmp=idximg
  jx=ix+1L & jy=iy+0L & idximg[jx,jy]=idximg[jx,jy] < (idximg[ix,iy])
  jx=ix-1L & jy=iy+0L & idximg[jx,jy]=idximg[jx,jy] < (idximg[ix,iy])
  jx=ix+0L & jy=iy+1L & idximg[jx,jy]=idximg[jx,jy] < (idximg[ix,iy])
  jx=ix+0L & jy=iy-1L & idximg[jx,jy]=idximg[jx,jy] < (idximg[ix,iy])
  jx=ix+1L & jy=iy+1L & idximg[jx,jy]=idximg[jx,jy] < (idximg[ix,iy])
  jx=ix-1L & jy=iy+1L & idximg[jx,jy]=idximg[jx,jy] < (idximg[ix,iy])
  jx=ix+1L & jy=iy-1L & idximg[jx,jy]=idximg[jx,jy] < (idximg[ix,iy])
  jx=ix-1L & jy=iy-1L & idximg[jx,jy]=idximg[jx,jy] < (idximg[ix,iy])
  ;
  nuncorr=total(tmp-idximg)
  depth=depth+1L
  deepimg=deepimg + ((tmp-idximg) gt 0)
  ;
  if vv gt 0 then kilroy
  if vv gt 10 then tvscl,idximg
  ;if vv gt 10 then tvscl,deepimg
endwhile

if vv gt 2000 then stop,'HALTING; type .CON to continue'

;	output
oimg=long(0*idximg)
;	this is what it would look like for readability's sake
;ilab=idximg[uniq(idximg,sort(idximg))] & nlabels=n_elements(ilab)
;areas=fltarr(nlabels+1) & areas[0]=nx*ny-midx	;how many pixels in the background?
;for i=0L,nlabels-1L do begin
;  ok=where(idximg eq ilab[i] and ilab[i] lt midx,mok)
;  areas[i+1L]=mok & if mok gt 0 then oimg[ok]=i+1L
;endfor
;	and this is what it looks like for speed's sake
h=histogram(idximg,min=0,reverse_indices=ri)
nh=n_elements(h)	;generally equals max(idximg)+1
rvec=ri[0L:nh] & drvec=rvec[1:*]-rvec & ilab=where(drvec ne 0,nlabels)
areas=fltarr(nlabels)
;
for i=0L,nlabels-1L do begin
  j=ilab[i]
  areas[i]=drvec[j]
  oimg[ri[ri[j]:ri[j+1L]-1L]]=i
endfor
;
if keyword_set(arthr) then begin
  if arthr gt 0 and arthr lt 1 then begin
    tmpar=areas[1:*]
    os=sort(tmpar) & careas=total(tmpar[os],/cumul)
    if nlabels gt 1 then begin
      cf=findgen(nlabels-1L)/(nlabels-1L)
      if nlabels gt 2 then areacut=interpol(tmpar[os],cf,arthr) else $
	areacut=0.
    endif else areacut=areas[0]*0.9
  endif else areacut=arthr
  for i=0L,nlabels-1L do begin
    j=ilab[i]
    if areas[i] lt areacut then oimg[ri[ri[j]:ri[j+1L]-1L]]=0
  endfor
endif

if vv gt 1000 then stop,'HALTING; type .CON to continue'

return,oimg
end
