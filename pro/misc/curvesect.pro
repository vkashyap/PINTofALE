function curvesect,y1,y2,xpt, verbose=verbose, _extra=e
;+
;function	curvesect
;	finds and returns the point(s) of intersection of two curves
;
;syntax
;	xy=curvesect(y1,y2,xpt,verbose=verbose)
;
;parameters
;	y1	[INPUT; required] points of curve Y1(X)
;	y2	[INPUT; required] points of curve Y2(X)
;		* Y1 and Y2 must be defined on the same grid
;	xpt	[INPUT] x-coordinates of Y1 and Y2
;		* if not supplied, uses the array index
;
;keywords
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	only returns the intersections within the domain of XPT
;	inputs must be sorted in X
;
;description
;	a very simple algorithm that checks whether two curves
;	have any crossings, and if they do, finds the point of
;	intersection using the line segments that connect across
;	the crossings.  Don't push it.
;
;example
;	.run curvesect
;
;history
;	Vinay Kashyap (2012dec)
;-

;	usage
ok='ok' & np=n_params() & n1=n_elements(y1) & n2=n_elements(y2)
if np lt 2 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='Y1 is not defined' else $
  if n2 eq 0 then ok='Y2 is not defined' else $
   if n1 ne n2 then ok='Y1 and Y2 must be on same grid' else $
    if n1 lt 2 then ok='curves must have at least 2 points'
if ok ne 'ok' then begin
  print,'Usage: xy=curvesect(y1,y2,xpt,verbose=verbose)'
  print,'  finds and returns points of intersection between Y1 and Y2'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	initialize
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
xx=findgen(n1) & nx=n_elements(xpt)
if nx eq n1 then xx=xpt else if vv gt 0 then message,$
	'XPT incompatible with Y1 and Y2; using array indices',/informational

;	are there any intersections at all?
i1=y1 gt y2
ti1=total(i1)
if ti1 eq 0 or ti1 eq n1 then begin	;(one curve lies entirely above the other
  message,'No intersections in this range',/informational
  return,[!values.F_NAN,!values.F_NAN]
endif					;ti1=0 or ti1=N1)

;	how many intersections?
di1=i1[1:*]-i1
od1=where(di1 ne 0,nsec)
;if mod1 ne nsec then message,'BUG!'

;for each intersection, find the adjacent points and find the intersect of the segment
;	y-y10 = m1 * (x-x0) ==> y = m1*x + (y10-m1*x0)
;	similarly, y = m2*x + (y20-m2*x0)
;	so (m1-m2)*x_int = (y20-m2*x0)-(y10-m1*x0) ==> x_int = ((y20-y10)-x0*(m2-m1))/(m1-m2)
;	and y_int = m1*x_int + (y10-m1*x0)
xy=fltarr(2,nsec)
for i=0L,nsec-1L do begin
  k0=od1[i] & k1=k0+1
  dx=xx[k1]-xx[k0] & dy1=y1[k1]-y1[k0] & dy2=y2[k1]-y2[k0]
  m1=dy1/dx & m2=dy2/dx
  if m1 eq m2 then begin
    xy[*,i]=[xx[k0],y1[k0]]
    continue	;skip to next point
  endif
  xint=((y2[k0]-y1[k0])-xx[k0]*(m2-m1))/(m1-m2)
  yint=0.5*((m1*xint+(y1[k0]-m1*xx[k0])) + (m2*xint+(y2[k0]-m2*xx[k0])))
  xy[*,i]=[xint,yint]
endfor

if vv gt 10000 then stop,'HALTing; type .CON to continue'

return,xy
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;example usage

if not keyword_set(verbose) then verbose=1
xpt=findgen(20) & y1=xpt^2-20*x+100. & y2=1*xpt+40
plot,xpt,y1,psym=-1,xtitle='X',ytitle='Y',title='[PINTofALE] curvesect.pro',line=1 & oplot,xpt,y2,psym=-1,line=1

xy=curvesect(y1,y2,xpt,verbose=verbose)

oplot,xy[0,*],xy[1,*],psym=4,symsize=2

end
