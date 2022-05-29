function areapoly,xpt,ypt,xcen=xcen,ycen=ycen,nosort=nosort,$
    outx=outx,outy=outy,outar=outar,verbose=verbose, _extra=e
;+
;function	areapoly
;	compute and return area of an irregular polygon
;
;syntax
;	area=areapoly(xpt,ypt,xcen=xcen,ycen=ycen,/nosort,$
;	outx=outx,outy=outy,outar=outar,verbose=verbose)
;
;parameters
;	xpt	[INPUT; required] x-coordinates of points that make up polygonal shape
;	ypt	[INPUT; required] y-coordinates of points that make up polygonal shape
;		* sizes must match
;
;keywords
;	xcen	[I/O] x-coordinate of a point inside the polygon
;	ycen	[I/O] y-coordinate of a point inside the polygon
;		* if not given, *CEN are computed as mean(*PT) and returned on output
;		* this point _must_ be inside the polygon, otherwise the output
;		  will be meaningless
;		* if set to 0 (I*4, not 0L == I*8 or 0b ==  I*2), will get recomputed and overwritten
;	nosort	[INPUT] if set, does not reorder the points for increasing angular
;		coordinate
;	outx	[OUTPUT] the XPT in translated (and resorted) coordinates
;	outy	[OUTPUT] the YPT in translated (and resorted) coordinates
;	outar	[OUTPUT] the area of each triangle calculated sitting on the base
;		of (XCEN,YCEN)->(XPT,YPT)
;	verbose	[INPUT] controls chatter
;
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	the area of the polygon is the sum of the swept triangles that all
;	have vertices at some point inside the polygon.
;
;example
;	.run areapoly
;	which approximates a 16-edge polygon out of a (3,4) ellipse
;
;history
;	vinay kashyap (2013oct)
;	bug fix for polygons with indentations (VK; 2015sep)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(xpt) & ny=n_elements(ypt)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='XPT is undefined' else $
  if ny eq 0 then ok='YPT is undefined' else $
   if nx lt 3 then ok='too few XPT' else $
    if ny lt 3 then ok='too few YPT' else $
     if nx ne ny then ok='XPT and YPT are incompatible'
if ok ne 'ok' then begin
  print,'Usage: area=areapoly(xpt,ypt,xcen=xcen,ycen=ycen,/nosort,$
  print,'       outx=outx,outy=outy,outar=outar,verbose=verbose)
  print,'  compute and return area of an irregular polygon'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
cenx=mean(xpt,/double) & ceny=mean(ypt,/double)
nxc=n_elements(xcen) & nyc=n_elements(ycen)
if nxc gt 1 and vv gt 1 then message,'XCEN is an array; only first element will be used',/informational
if nyc gt 1 and vv gt 1 then message,'YCEN is an array; only first element will be used',/informational
if nxc gt 0 then begin
  if size(xcen[0],/type) eq 2 and xcen[0] eq 0 then xcen=cenx else cenx=xcen[0]
endif
if nyc gt 0 then begin
  if size(ycen[0],/type) eq 2 and ycen[0] eq 0 then ycen=ceny else ceny=ycen[0]
endif
if nxc eq 0 then xcen=cenx
if nyc eq 0 then ycen=ceny
;
xx=[xpt[*]]-cenx & yy=[ypt[*]]-ceny
tht=(atan(yy,xx)+2.*!pi) mod (2.*!pi)
if not keyword_set(nosort) then begin
  os=sort(tht) & xx=xx[os] & yy=yy[os] & tht=tht[os]
endif
outx=xx & outy=yy

;	now compute the areas of each triangle sitting on the base of each point
d1=sqrt(xx^2+yy^2) & d2=shift(d1,1) & dtheta=tht-shift(tht,1)
outar=d1*d2*sin(dtheta)/2.
area=total(outar,/nan,/double)

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,area
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example calling sequence

;	usage
area=areapoly()

;	test case
aa=4. & bb=3.
rth=randomu(seed,16)*2*!pi & xpt=aa*cos(rth) & ypt=bb*sin(rth) & plot,xpt,ypt,psym=1,title=!pi*aa*bb,xr=[-6,6],yr=[-4,4],/xs,/ys
xcen=0 & ycen=0
area=areapoly(xpt,ypt,xcen=xcen,ycen=ycen,outx=outx,outy=outy,outar=outar,verbose=verbose)
plot,outx,outy,psym=-1,title=!pi*4*3 & oplot,[xcen],[ycen],psym=4,symsize=2
print,'calculated area, xcen, ycen:',area,xcen,ycen
print,'true area of full ellipse is '+strtrim(!pi*aa*bb,2)
print,'calculated area is always smaller because points are connected by straight line segments'

end
