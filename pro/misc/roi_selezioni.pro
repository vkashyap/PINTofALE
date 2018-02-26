function roi_selezioni,image,iseed,thresh=thresh,jitter=jitter,hadrian=hadrian,$
	verbose=verbose, _extra=e
;+
;function	roi_selezioni
;	simple tool to pick out a set of pixels defining a region of
;	interest in a 2D image;	returns selection as 1D index array
;
;syntax
;	iroi=roi_selezioni(image,iseed,thresh=thresh,jitter=jitter,$
;	hadrian=hadrian,verbose=verbose)
;
;parameters
;	image	[INPUT; required] image from which to select region of interest
;	iseed	[INPUT] pixel indices to act as a seed selection
;		* if given, will use these points as the selected
;		  pixels and work in batch processing mode, i.e.,
;		  will NOT be interactive
;
;keywords
;	thresh	[INPUT] if given, then filters the input image such
;		that IMAGE > THRESH
;		* if not given, assumed to be 0
;	jitter	[INPUT] sensitivity in distinguishing between click and drag
;		* default is 1.5 pixels
;	hadrian	[INPUT] thus far and no further.  if set, places a fixed
;		boundary around the area selected with the cursor, i.e.,
;		connected regions are not followed out of this area.  set
;		this to pick out individual pixels, for instance.
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Aug06)
;	added parameter ISEED and keywords THRESH, HADRIAN _EXTRA;
;	  allowed selection extending outside cursor selected areas;
;	  changed name from pickoutroi to roi_selezioni (VK; May07)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(image) & szi=size(image)
nx=szi[1] & ny=szi[2]
if np eq 0 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='IMAGE is undefined' else $
  if szi[0] ne 2 then ok='IMAGE must be a 2-D array' else $
   if nx eq 1 then ok='IMAGE x-axis is collapsed' else $
    if ny eq 1 then ok='IMAGE y-axis is collapsed'
if ok ne 'ok' then begin
  print,'Usage: iroi=roi_selezioni(image,iseed,jitter=jitter,verbose=verbose)'
  print,'  returns selected area as an array of indices' 
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	interactive or not?
ibatch=0
if n_elements(iseed) gt 0 then begin
  if iseed[0] ne -1 then ibatch=1
endif
if keyword_set(ibatch) then begin
  isx=iseed mod nx
  isy=iseed/nx
  xmin0=min(isx,max=xmax0)
  ymin0=min(isy,max=ymax0)
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
dthresh=1.5 & if n_elements(jitter) gt 0 then dthresh=jitter[0]
thrlev=0 & if keyword_set(thresh) then thrlev=thresh[0]

;	initialize
baseimg=image & imgmax=max(image)
xpix=findgen(nx) # (fltarr(ny)+1)
ypix=(fltarr(nx)+1) # findgen(ny)

;	what is the largest screen size we can have?
nxre=nx & nyre=ny
if ibatch eq 0 then begin
  tmp=obj_new('IDLgrWindow')
  tmp->GetProperty, screen_dimensions=ss
  obj_destroy,tmp
  xsizemax=ss[0]-40 & ysizemax=ss[1]
  rescale=1
  if nx gt xsizemax or ny gt ysizemax then begin
    xrescale=fix(nx/xsizemax)+1
    yrescale=fix(ny/ysizemax)+1
    rescale=xrescale > yrescale
    if xsizemax-(nx/rescale) ne 0 or ysizemax-(ny/rescale) ne 0 then $
	usecongrid=1
  endif
  nxre=nx/rescale & nyre=nx/rescale
  window,xsize=nxre,ysize=nyre
endif else rescale=1

;	how to use it
cc=[	'',$
	'	Interactively pick out a region of interest',$
	'LEFT DRAG : select rectangular area and associated pixels',$
	'MIDDLE DRAG : de-select rectangular area and associated pixels',$
	'RIGHT CLICK or DRAG : exit program',$
	'LEFT/MIDDLE CLICK : selects/deselects pixels connected to single pixel',$
	'']
if ibatch eq 0 then for i=0,n_elements(cc)-1 do print,cc[i]

;	use the cursor to select appropriate region
go_on=1 & miroi=0L & iroi=-1L
while go_on do begin		;{keep on selectin

  ;	display image
  tmporimg=baseimg
  if miroi gt 0 then tmporimg[iroi]=max(tmporimg)+5
  if rescale eq 1 then tmpimg=tmporimg else begin
    if keyword_set(usecongrid) then tmpimg=congrid(tmporimg,nxre,nyre) else $
      tmpimg=rebin(tmporimg,nxre,nyre)
  endelse
  ;	if any regions have already been selected, show it
  tvscl,tmpimg

  ;	what to do next?
  if ibatch eq 0 then begin
    cursor,xd,yd,/down,/dev & cursor,xu,yu,/up,/dev
    ;
    dx=xu-xd & dy=yu-yd & dd=sqrt(float(dx)^2+float(dy)^2)
    xmin=min([xd,xu],max=xmax) & ymin=min([yd,yu],max=ymax)
    if dd lt dthresh then drag=0 else drag=1
    if !MOUSE.BUTTON eq 1 then click='left' else $
     if !MOUSE.BUTTON eq 2 then click='middle' else $
      if !MOUSE.BUTTON eq 4 then click='right'
    xmin=xmin*rescale & xmax=xmax*rescale
    ymin=ymin*rescale & ymax=ymax*rescale
  endif else begin
    drag=0 & click=''
    xmin=xmin0 & xmax=xmax0
    ymin=ymin0 & ymax=ymax0
  endelse

  ;	cases
  if vv gt 5 then begin
    ok=click
    if drag eq 0 then ok=ok+'+click' else ok=ok+'+drag'
    kilroy,dot=ok+' '
  endif

  if ibatch eq 1 then click='left'

  if click eq 'left' then begin			;(left click or drag
    ;	select all pixels in the box defined by the click and drag
    if keyword_set(ibatch) then begin
     ok=iseed & mok=n_elements(iseed)
    endif else $
      ok=where(xpix ge xmin and xpix le xmax and $
      ypix ge ymin and ypix le ymax and baseimg gt thrlev,mok)
    if mok gt 0 then begin
      if keyword_set(hadrian) then ok2=ok else begin
        tmp=baseimg[ok] & ut=tmp[uniq(tmp,sort(tmp))] & mut=n_elements(ut)
        ok2=where(baseimg eq ut[0])
        for jj=1,mut-1L do ok2=[ok2,where(baseimg eq ut[jj])]
      endelse
      if miroi eq 0 then begin
	iroi=ok2 & miroi=n_elements(ok2)
      endif else begin
	iroi=[iroi,ok2] & iroi=iroi[uniq(iroi,sort(iroi))]
	miroi=n_elements(iroi)
      endelse
    endif
    if ibatch eq 0 and drag eq 0 and vv gt 2 then $
	for i=0,n_elements(cc)-1 do print,cc[i]
  endif						;left click or drag)

  if click eq 'middle' then begin		;(middle click or drag
    ;	deselect all pixels in the box defined by the click and drag
    ok=where(xpix ge xmin and xpix le xmax and $
	ypix ge ymin and ypix le ymax and baseimg gt thrlev,mok)
    if mok gt 0 then begin
      if keyword_set(hadrian) then ok2=ok else begin
        tmp=baseimg[ok] & ut=tmp[uniq(tmp,sort(tmp))] & mut=n_elements(ut)
        ok2=where(baseimg eq ut[0])
        for jj=1,mut-1L do ok2=[ok2,where(baseimg eq ut[jj])]
      endelse
      tmp=0*fix(baseimg) & if miroi gt 0 then tmp[iroi]=1
      tmp[ok2]=0
      iroi=where(tmp ne 0,miroi)
    endif
    if drag eq 0 and vv gt 2 then for i=0,n_elements(cc)-1 do print,cc[i]
  endif						;middle drag)

  if click eq 'right' then begin		;(right click or drag
    ;	quit
    go_on=0
  endif						;right click or drag)
  if ibatch eq 1 then go_on=0

endwhile			;GO_ON}

return,iroi
end
