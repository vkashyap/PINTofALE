function mk_lceclipse,x,midpt,hwidth,depth,pder,angfall=angfall,angrise=angrise,$
	steep=steep,ybase=ybase,yoff=yoff,eps=eps, _extra=e
;+
;function	mk_lceclipse
;	returns the light curve of a hard eclipse of an almost point source
;
;syntax
;	lc=mk_lceclipse(x,midpt,hwidth,depth,pder,angfall=angfall,angrise=angrise,$
;	ybase=ybase,yoff=yoff,/steep)
;
;parameters
;	x	[INPUT array; required] where LC(X) must be computed
;	midpt	[INPUT; default: mid_point(X)] position of central point of eclipse
;	hwidth	[INPUT; default: 0.1*range(X)] half-width of eclipse
;	depth	[INPUT; default: 1] depth of eclipse
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if 5 parameters are supplied in call.
;		* array of size [N(X),7], with columns containing the partial
;		  derivatives wrt MIDPT,WIDTH,DEPTH,ANGFALL,ANGRISE,YBASE,YOFF
;
;keywords
;	angfall	[INPUT; default=3*!pi/2] angle of lightcurve descent in radians
;	angrise	[INPUT; default=!pi/2] angle of lightcurve ascent in radians
;		* NOTE: if STEEP is set, ANGFALL and ANGRISE are assumed to be slopes
;		  and the default values are set to -1000 and +1000 respectively
;	steep	[INPUT] if set, assumes that ANGFALL and ANGRISE are given as slopes
;	ybase	[INPUT; default=0] baseline emission
;	yoff	[INPUT; default=0] offset on the other side of eclipse
;	eps	[INPUT; default=1e-6] a small number
;	_extra	[JUNK] here only to prevent crashing
;
;description
;	The eclipse curve is:
;	slp1=tan(ANGFALL) & slp2=tan(ANGRISE)
;	LC(X)	= YBASE+DEPTH	for X < MIDPT-HWIDTH+DEPTH/slp1
;		= slp1*X+(YBASE-slp1*MIDPT+slp1*HWIDTH)
;			for MIDPT-HWIDTH-DEPTH/slp1 .LE. X .LE. MIDPT-HWIDTH
;		= YBASE	for MIDPT-HWIDTH < X < MIDPT+HWIDTH
;		= slp2*X+(YBASE-slp2*MIDPT-slp2*HWIDTH)
;			for MIDPT+HWIDTH .LE. X .LE. MIDPT+HWIDTH+(DEPTH+YOFF)/slp2
;		= YBASE+DEPTH+YOFF	for X > MIDPT+HWIDTH+(DEPTH+YOFF)/slp2
;	
;usage summary
;	* call as a function
;	* generates hard eclipse light curve at specified points
;
;subroutines
;	NONE
;
;example
;	x=findgen(100)
;	y=mk_lceclipse(x,60,20,0.77,pder,ybase=0.4,yoff=0.11,$
;	angfall=178*!pi/180.,angrise=10.*!pi/180.) & plot,x,y
;
;history
;	vinay kashyap (MayMMVII; based on MK_LORENTZ)
;	added keyword STEEP, also corrected bug of pderiv wrt ANG (VK; AugMMVII)
;-

np=n_params()
if np lt 1 then begin
  print,'Usage: lc=MK_LCECLIPSE(X,midpt,hwidth,depth,pder,angfall=angfall,angrise=angrise,$'
  print,'       ybase=ybase,yoff=yoff,/steep)'
  print,'  returns the light curve of a hard eclipse of an almost point source'
  return,-1L
endif

;initialize
nx=n_elements(x) & x0=x[nx/2] & mxx=max(x,min=mnx)

;figure out the defaults
if np lt 4 then ydepth=1. else ydepth=depth[0]
if np lt 3 then xhwidth=0.1*(mxx-mnx) else xhwidth=hwidth[0]
if np lt 2 then xmidpt=x0 else xmidpt=midpt[0]
if keyword_set(steep) then begin
  if keyword_set(angfall) then afall=angfall[0] else afall=-1000.
  if keyword_set(angrise) then arise=angrise[0] else arise=1000.
endif else begin
  if keyword_set(angfall) then afall=(angfall[0] mod (2*!pi)) else afall=3*!pi/2.
  if keyword_set(angrise) then arise=(angrise[0] mod (2*!pi)) else arise=!pi/2.
endelse
if keyword_set(ybase) then basey=ybase[0] else basey=0.
if keyword_set(yoff) then yoffset=yoff[0] else yoffset=0.

;	compute function
if keyword_set(steep)  then begin
  slp1=afall & slp2=arise
endif else begin
  slp1=tan(afall) & slp2=tan(arise)
endelse
p1=xmidpt-xhwidth+ydepth/slp1 & o1=where(x le p1,mo1)
p2=xmidpt-xhwidth & o2=where(x ge p1 and x le p2,mo2)
p3=xmidpt+xhwidth & o3=where(x ge p2 and x le p3,mo3)
p4=xmidpt+xhwidth+(ydepth+yoffset)/slp2 & o4=where(x ge p3 and x le p4,mo4) & o5=where(x ge p4,mo5)
;
y=0.*x 
if mo1 gt 0 then y[o1]=basey+ydepth
if mo2 gt 0 then y[o2]=slp1*X[o2]+(basey-slp1*xmidpt+slp1*xhwidth)
if mo3 gt 0 then y[o3]=basey
if mo4 gt 0 then y[o4]=slp2*X[o4]+(basey-slp2*xmidpt-slp2*xhwidth)
if mo5 gt 0 then y[o5]=basey+ydepth+yoffset

;	compute partial derivatives
if np ge 5 then begin
  pder = fltarr(nx,7)
  ;	partial wrt MIDPT
  if mo2 gt 0 then pder[o2,0]=-slp1
  if mo4 gt 0 then pder[o4,0]=-slp2
  ;	partial wrt WIDTH
  if mo2 gt 0 then pder[o2,1]=slp1
  if mo4 gt 0 then pder[o4,1]=-slp2
  ;	partial wrt DEPTH
  if mo1 gt 0 then pder[o1,2]=1.
  if mo5 gt 0 then pder[o5,2]=1.
  ;	partial wrt ANGFALL
  if keyword_set(steep) then begin
    if mo2 gt 0 then pder[o2,3]=X[o2]-xmidpt+xhwidth
  endif else begin
    if mo2 gt 0 then pder[o2,3]=(1.+slp1^2)*(x[o2]+xmidpt+xhwidth) ;bug correction -- originally the () enclosed only the X
  endelse
  ;	partial wrt ANGRISE
  if keyword_set(steep) then begin
    if mo4 gt 0 then pder[o4,3]=X[o2]-xmidpt-xhwidth
  endif else begin
    if mo4 gt 0 then pder[o4,4]=(1.+slp2^2)*(X[o4]+xmidpt+xhwidth) ;bug correction -- originally the () enclosed only the X
  endelse
  ;	partial wrt YBASE
  pder[*,5]=1.
  ;	partial wrt YOFF
  if mo5 gt 0 then pder[o5,6]=1.
endif

return,y
end
