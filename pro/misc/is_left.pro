function is_left,P0,P1,P2
;+
;function	is_left
;	tests if a point left, on, or right of a 2D line and returns
;	>0 if left, =0 if on, and <0 if right
;
;	translated from C++ code of same name by Dan Sunday
;	see http://geomalgorithms.com/a01-_area.html
;
;syntax
;	pp=is_left(P0,P1,P2)
;
;parameters	[ALL INPUT, ALL REQUIRED]
;	P0	2-element array defining one of the points on line
;	P1	2-element array defining the other point on line
;	P2	2,N element array of points to test for whether it is to
;		left (>0), on (=0), or right (<0) of line P0-P1
;		* in all cases, PX[0,..] is assumed to be the x-coordinate
;		  and PX[1,..] is assumed to be the y-coordinate
;
;keywords
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2016.08)
;	changed name from isLeft to is_left because got tired of .comp'ing endlessly (VK; 2021.06)
;-

;	usage
ok='ok' & np=n_params()
n0=n_elements(P0) & n1=n_elements(P1)
n2=n_elements(P2) & sz2=size(P2)
if np lt 3 then ok='Insufficient parameters' else $
 if n0 eq 0 then ok='P0 is not defined' else $
  if n1 eq 0 then ok='P1 is not defined' else $
   if n2 eq 0 then ok='P2 is not defined' else $
    if n0 ne 2 then ok='P0 must be (x,y) pair' else $
     if n1 ne 2 then ok='P1 must be (x,y) pair' else $
      if sz2[0] eq 1 and n2 ne 2 then ok='P2 must be (x,y) pair' else $
       if sz2[0] gt 2 then ok='P2 cannot be understood' else $
        if sz2[0] eq 2 and sz2[1] ne 2 then ok='P2 must be 2xN array of (x,y)'
if ok ne 'ok' then begin
  print,'Usage: pp=is_left(P0,P1,P2)'
  print,'  tests if P2 is left (>0), on (=0), or right (<0) of line P0-P1'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	initialize
nn=sz2[2] & if sz2[0] eq 1 then nn=1L
if sz2[0] eq 1 then begin
  zx=P2[0] & zy=P2[1]
endif else begin
  zx=reform(P2[0,*]) & zy=reform(P2[1,*])
endelse
pp=(P1[0]-P0[0])*(zy-P0[1]) - (zx-P0[0])*(P1[1]-P0[1])

return,pp
end

;	original C++ code pasted below, in case reference URL goes stale
;	http://geomalgorithms.com/a01-_area.html
;// Copyright 2000 softSurfer, 2012 Dan Sunday
;// This code may be freely used and modified for any purpose
;// providing that this copyright notice is included with it.
;// iSurfer.org makes no warranty for this code, and cannot be held
;// liable for any real or imagined damage resulting from its use.
;// Users of this code must verify correctness for their application.
; 
;// isLeft(): test if a point is Left|On|Right of an infinite 2D line.
;//    Input:  three points P0, P1, and P2
;//    Return: >0 for P2 left of the line through P0 to P1
;//          =0 for P2 on the line
;//          <0 for P2 right of the line
;inline int
;isLeft( Point P0, Point P1, Point P2 )
;{
;    return ( (P1.x - P0.x) * (P2.y - P0.y)
;           - (P2.x - P0.x) * (P1.y - P0.y) );
;}
;//===================================================================


