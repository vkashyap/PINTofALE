pro locus_ellipse,xout,yout,xcen=xcen,ycen=ycen,amajor=amajor,aminor=aminor,$
	theta=theta,npt=npt,verbose=verbose, _extra=e
;+
;procedure	locus_ellipse
;	make a curve describing the locus of an ellipse 
;
;syntax
;	locus_ellipse,xout,yout,xcen=xcen,ycen=ycen,amajor=amajor,$
;	aminor=aminor,theta=theta,npt=npt,verbose=verbose
;
;parameters
;	xout	[OUTPUT; required] x-coordinates of the locus
;	yout	[OUTPUT; required] y-coordinates of the locus
;
;keywords
;	xcen	[INPUT; default=0] x-coordinate of the center
;	ycen	[INPUT; default=0] y-coordinate of the center
;	amajor	[INPUT; default=1] length of the semi-major axis
;	aminor	[INPUT; default=1] length of the semi-minor axis
;		* if AMINOR is not defined, it is set to AMAJOR
;		  to make a circle
;		* if AMINOR>AMAJOR, no harm done, program will
;		  run regardless of mathematical purity
;	theta	[INPUT; default=0 degrees] anticlockwise angle made
;		by AMAJOR with x-axis
;	npt	[INPUT; default=101] number of points along the locus,
;		spaced equally along THETA
;		* note that the first and the last points are identical
;		  unless NPT < 0, in which case abs(NPT) points are
;		  returned without the overlap
;	verbose	[INPUT; default=0] controls chatter
;	
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Mar2008)
;-

;	usage
ok='ok' & np=n_params()
if np eq 0 then ok='Insufficient parameters' else $
 if np eq 1 then ok='YOUT is not defined'
if ok ne 'ok' then begin
  print,'Usage: locus_ellipse,xout,yout,xcen=xcen,ycen=ycen,amajor=amajor,$'
  print,'       aminor=aminor,theta=theta,npt=npt,verbose=verbose'
  print,'  make a curve describing the locus of an ellipse'
  if np ne 0 then message,ok,/informational
  return
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	initialize
nx0=n_elements(xcen) & ny0=n_elements(ycen)
naa=n_elements(amajor) & nbb=n_elements(aminor) & nth=n_elements(theta)
nm=n_elements(npt)
;
if nx0 eq 0 then x0=0. else x0=float(xcen[0])
if ny0 eq 0 then y0=0. else y0=float(ycen[0])
if naa eq 0 then aa=1. else aa=abs(float(amajor[0]))
if aa eq 0 then begin & if vv gt 0 then message,'AMAJOR cannot be 0; setting to 1',/informational & aa=1. & endif
if nbb eq 0 then bb=aa else bb=abs(float(aminor[0]))
if bb eq 0 then begin & if vv gt 0 then message,'AMINOR cannot be 0; setting to 1',/informational & bb=1. & endif
if nth eq 0 then th=0. else th=float(theta[0])
foldover=1
if nm eq 0 then nn=101 else begin
  nn=long(npt[0])
  if nn eq 0 then nn=101
  if nn lt 0 then begin & foldover=0 & nn=abs(nn) & endif
  if nn eq 1 then foldover=0
endelse

;	make the ellipse
tt=findgen(nn)*360./nn & if keyword_set(foldover) then tt=tt*float(nn)/float(nn-1.)
rr=tt*!pi/180.
xx=aa*cos(rr) & yy=bb*sin(rr)

;	rotate
thr=th*!pi/180.
ax=[0.,aa] & ay=[0.,0.]
if th ne 0 then begin
  xout=xx*cos(thr)-yy*sin(thr)
  yout=xx*sin(thr)+yy*cos(thr)
  xa=ax*cos(thr)-ay*sin(thr) & ya=ax*sin(thr)-ay*cos(thr) & ax=xa & ay=ya
endif else begin
  xout=xx & yout=yy
endelse

;	translate
xout=xout+x0 & yout=yout+y0
ax=ax+x0 & ay=ay+y0

;	plot
if vv gt 10 then begin
  xmin=min([xx,xout,yy,yout],max=xmax) & ymin=xmin & ymax=xmax
  plot,[0],xtitle='X',ytitle='Y',xrange=[xmin,xmax],yrange=[ymin,ymax],/nodata
  peasecolr
  oplot,[xmin-abs(xmax),xmax+abs(xmin)],[ycen,ycen],line=1 & oplot,[xcen,xcen],[ymin-abs(ymax),ymax+abs(ymin)],line=1
  oplot,xx,yy,line=1,col=1
  oplot,xout,yout,col=2,thick=2
  ;oplot,xx,tan(th*!pi/180.)*xx,line=1,col=3
  oplot,ax,ay,line=1,col=3
endif

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return
end
