function get_continuum,x,y
;+
;function	get_continuum
;		returns the values of the continuum as determined
;		by fitting a spline curve to a set of interactively
;		selected points.
;
;parameters	x	[INPUT; required] points at which "reference curve"
;			is defined
;		y	[INPUT] y(x) (if absent, X is taken to be Y and X
;			is set to the index positions)
;
;restrictions
;	requires X display and mouse interaction with plot
;
;description
;	select a bunch of points on a plot, get a spline fit with the
;	default "tension", use differences between computed curve and
;	input curve to get new value of "tension", recompute the curve.
;
;history
;	vinay kashyap (Oct96)
;-

message,'OBSOLETE!',/informational

np=n_params(0)
if np lt 1 then begin
  print,'Usage: yc=get_continuum(x,y)'
  print,'  return a "continuum", spline fit to points selected interactively'
  return,[-1L]
endif
xx=x							;save
if np eq 2 then yy=y else yy=xx				;save/set
ny=n_elements(yy) & if np eq 1 then xx=lindgen(ny)	;reset

;plot
plot,xx,yy,psym=10

;get list of points
xp=[0.] & yp=[0.] & go_on=1
while go_on eq 1 do begin
  cursor,x0,y0,/down,/data
  mbutton=!err
  if mbutton eq 2 then begin			;middle button -- delete
    nn=n_elements(xp)
    if nn gt 1 then begin
      oplot,[xp(nn-1)],[yp(nn-1)],psym=1,col=0
      xp=xp(0:nn-2) & yp=yp(0:nn-2)
    endif
  endif else begin
    if mbutton eq 4 then begin			;right button -- quit
      go_on=0
    endif else begin				;left button -- pick
      xp=[xp,x0] & yp=[yp,y0]
      ;oplot,[x0],[y0],psym=1
      plot,xx,yy,psym=10 & oplot,xp(1:*),yp(1:*),psym=1
    endelse
  endelse
  nn=n_elements(xp)
  if nn eq 2 then oplot,xx,yy*0.+yp(1),linestyle=1
  if nn eq 3 then begin
    if xp(2) ne xp(1) then slope=(yp(2)-yp(1))/(xp(2)-xp(1)) else slope=1e10
    yc=yp(1)+slope*(xx-xp(1)) & oplot,xx,yc,linestyle=1
  endif
  if nn gt 3 then begin
    px=xp(1:*) & py=yp(1:*) & oo=sort(px)
    yc=spline(px(oo),py(oo),xx) & oplot,xx,yc,linestyle=1
  endif
endwhile
nn=n_elements(xp)-1
if nn lt 1 then begin & message,'interact, dammit!',/info & return,yy & endif
if nn eq 1 then return,0*y+yp(1)		;constant
if nn eq 2 then begin
  if xp(1) ne xp(0) then slope=(yp(1)-yp(0))/(xp(1)-xp(0)) else slope=1e10
  yc=yp(0)+slope*(xx-xp(0))
  plot,xx,yy,psym=10 & oplot,xx,yc,linestyle=1
  return,yc
endif
xp=xp(1:*) & yp=yp(1:*)
oo=sort(xp) & xp=xp(oo) & yp=yp(oo)

;do a spline fit
yc=spline(xp,yp,xx)
plot,xx,yy,psym=10 & oplot,xx,yc,linestyle=1

;residuals & sigma
res=yc-yy & sig=sqrt((moment(res))(1)) & print,sig

;refit
yc=spline(xp,yp,xx,sig)
plot,xx,yy,psym=10 & oplot,xx,yc

return,yc
end
