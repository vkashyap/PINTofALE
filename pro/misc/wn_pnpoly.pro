function wn_pnpoly,x0,y0,xp,yp,eps=eps, _extra=e
;+
;function	wn_pnpoly
;	returns the winding number test for whether a point P0 is
;	inside (wn>0) or outside (wn=0) of polygon Pp
;
;	translated from C++ code of same name by Dan Sunday
;	see http://geomalgorithms.com/a03-_inclusion.html#wn_PnPoly()
;
;syntax
;	wn=wn_pnpoly(x0,y0,xp,yp,eps=eps)
;
;parameters	[ALL INPUT, ALL REQUIRED]
;	x0	x-coordinate of point to be tested
;	y0	y-coordinate of point to be tested
;	xp	array of x-coordinates of vertices of polygon
;	yp	array of y-coordinates of vertices of polygon
;		* XP and YP must have at least 3 points, else what's the point?
;		* the order is taken as is to be the order of vertices
;		* the algorithm assumes that the last point is the same as
;		  the first, so if it is not, the array gets inflated to
;		  close the polygon
;
;keywords
;	eps	[INPUT; default=1e-6] a small number
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2016.08)
;-

;	usage
ok='ok' & np=n_params()
nx0=n_elements(x0) & ny0=n_elements(y0)
nxp=n_elements(xp) & nyp=n_elements(yp)
if np lt 4 then ok='Insufficient parameters' else $
 if nx0 eq 0 then ok='X0 is not defined' else $
  if ny0 eq 0 then ok='Y0 is not defined' else $
   if nxp eq 0 then ok='Xp is not defined' else $
    if nyp eq 0 then ok='Yp is not defined' else $
     if nx0 ne ny0 then ok='X0 and Y0 are not compatible' else $
      if nxp lt 3 then ok='polygon has insufficient vertices' else $
       if nxp ne nyp then ok='Xp and Yp are not compatible'
if ok ne 'ok' then begin
  print,'Usage: wn=wn_pnpoly(x0,y0,xp,yp,eps=eps)'
  print,'  returns the winding number test for whether (x0,y0) is'
  print,'  inside (>0) or outside (=0) of polygon (xp,yp)'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	initialize winding number counter
wn=lonarr(nx0)
if nx0 eq 1 then P0=[x0,y0] else P0=make_array(2,nx0,value=0*x0[0])

;	what's a small number?
ee=1e-7 & if keyword_set(eps) then ee=double(abs(eps[0]))

;	close the polygon
d0n=sqrt((xp[0]-xp[nx0-1L])^2+(yp[0]-yp[ny0-1L])^2)
if d0n gt ee then begin
  xz=[xp,xp[0]] & yz=[yp,yp[0]] & nn=nxp+1L
endif else begin
  xz=xp & yz=yp & nn=nxp
endelse

;	loop through all edges of the polygon
for i=0L,nn-2L do begin	;{edge from P(i) to P(i+1)

  ;	figure out whether P0 is left or right of this edge
  PA=[xz[i],yz[i]] & PB=[xz[i+1L],yz[i+1L]]
  pp=isLeft(PA,PB,P0)
  	;note that P0 here is P2 inside of isLeft()

  for j=0L,nx0-1L do begin	;{for each X0,Y0

    if yz[i] le y0[j] then begin	;(start y .LE. Py
      if yz[i+1L] gt y0[j] then begin	;(an upward crossing
        if pp[j] gt 0 then wn[j]=wn[j]+1L	;P0 is left of edge and has a valid up intersect
      endif				;)
    endif else begin		;)(start y > Py
      if yz[i+1L] le y0[j] then begin	;(a downward crossing
        if pp[j] lt 0 then wn[j]=wn[j]-1L	;P0 is right of edge and has a valid down intersect
      endif				;)
    endelse			;)

  endfor			;J=0,NX0-1}

endfor			;I=0,NX}

return,wn
end

;	original C++ code pasted below, in case reference URL goes stale
;	http://geomalgorithms.com/a03-_inclusion.html#wn_PnPoly()
;
;// Copyright 2000 softSurfer, 2012 Dan Sunday
;// This code may be freely used and modified for any purpose
;// providing that this copyright notice is included with it.
;// SoftSurfer makes no warranty for this code, and cannot be held
;// liable for any real or imagined damage resulting from its use.
;// Users of this code must verify correctness for their application.
;
;// wn_PnPoly(): winding number test for a point in a polygon
;//      Input:   P = a point,
;//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
;//      Return:  wn = the winding number (=0 only when P is outside)
;int
;wn_PnPoly( Point P, Point* V, int n )
;{
;    int    wn = 0;    // the  winding number counter
;
;    // loop through all edges of the polygon
;    for (int i=0; i<n; i++) {   // edge from V[i] to  V[i+1]
;        if (V[i].y <= P.y) {          // start y <= P.y
;            if (V[i+1].y  > P.y)      // an upward crossing
;                 if (isLeft( V[i], V[i+1], P) > 0)  // P left of  edge
;                     ++wn;            // have  a valid up intersect
;        }
;        else {                        // start y > P.y (no test needed)
;            if (V[i+1].y  <= P.y)     // a downward crossing
;                 if (isLeft( V[i], V[i+1], P) < 0)  // P right of  edge
;                     --wn;            // have  a valid down intersect
;        }
;    }
;    return wn;
;}
;//===================================================================
