function circle_intersect,rad1,rad2,roff, verbose=verbose, _extra=e
;+
;function	circle_intersect
;	computes and returns the area of overlap between two circles of different radii, at a specified offset location
;
;syntax
;	area=circle_intersect(rad1,rad2,roff,verbose=verbose)
;
;parameters	ALL REQUIRED
;	rad1	[INPUT] radius of circle assumed wlog to be at origin
;	rad2	[INPUT] radius of overlapping circle
;	roff	[INPUT] distance of center of circle 2 from origin
;		* ROFF may be an array, in which case areas are calculated for all of them
;
;keywords
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to preventc crashing the program
;
;history
;	Vinay Kashyap (2019oct)
;-

;	usage
ok='ok' & np=n_params() & n1=n_elements(rad1) & n2=n_elements(rad2) & nd=n_elements(roff)
if np lt 3 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='RAD1 is undefined' else $
  if n2 eq 0 then ok='RAD2 is undefined' else $
   if nd eq 0 then ok='ROFF is undefined'
if ok ne 'ok' then begin
  print,'Usage: area=circle_intersect(rad1,rad2,roff,verbose=verbose)'
  print,'  compute and return area of overlap between two circles'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	gather inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
if vv gt 0 and n1 gt 1 then message,'RAD1 must be a scalar; ignoring RAD1[1:*]',/informational
if vv gt 0 and n2 gt 1 then message,'RAD2 must be a scalar; ignoring RAD2[1:*]',/informational
r1=rad1[0] & r2=rad2[0]
doff=roff[0] & if nd gt 1 then doff=roff

;	output
rmin=min([r1,r2],max=rmax) & maxarea=!pi*rmin^2
area=0.*doff-1

;	deal with simple cases
o1=where(doff gt r1+r2,mo1) & if mo1 gt 0 then area[o1]=0.
o2=where(rmax ge doff+rmin,mo2) & if mo2 gt 0 then area[o2]=maxarea

;	the full calculation
o3=where(area lt 0,mo3)
if mo3 gt 0 then begin
  alp1=acos((doff[o3]^2+r1^2-r2^2)/(2.*doff[o3]*r1))
  alp2=acos((doff[o3]^2+r2^2-r1^2)/(2.*doff[o3]*r2))
  area[o3]=r1^2*(alp1-cos(alp1)*sin(alp1)) + r2^2*(alp2-cos(alp2)*sin(alp2))
endif

if vv gt 1000 then stop,'halting; type .CON to continue'

return,area
end

;================================================================================
;example
; case 1: a smaller circle transits right through the center of a larger one
; case 2: a larger circle transits through the edge of a smaller one
; case 3: a smaller circle transits through the edge of a larger one
; case 4: a small circle transits through and just inside of a much larger one
;================================================================================

;	print usage
jnk=circle_intersect()
print,''

if not keyword_set(verbose) then verbose=1

;	case 1
rad1=2. & rad2=1. & xoff1=(findgen(37)-18)*0.2 & roff1=abs(xoff1)
area1=circle_intersect(rad1,rad2,roff1,verbose=verbose)/!pi
help,rad1,rad2,xoff1,roff1,area1

;	case 2
rad1=1. & rad2=2. & xoff2=(findgen(37)-18)*0.2 & yoff2=rad1 & roff2=sqrt(xoff2^2+yoff2^2)
area2=circle_intersect(rad1,rad2,roff2,verbose=verbose)/!pi
help,rad1,rad2,xoff2,roff2,yoff2,area2

;	case 3
rad1=2. & rad2=1. & xoff3=(findgen(37)-18)*0.2 & yoff3=rad1 & roff3=sqrt(xoff3^2+yoff3^2)
area3=circle_intersect(rad1,rad2,roff3,verbose=verbose)/!pi
help,rad1,rad2,xoff3,roff3,yoff3,area3

;	case 4
rad1=5. & rad2=1. & xoff4=(findgen(37)-18)*0.2 & yoff4=rad1-rad2 & roff4=sqrt(xoff4^2+yoff4^2)
area4=circle_intersect(rad1,rad2,roff4,verbose=verbose)/!pi
help,rad1,rad2,xoff4,roff4,yoff4,area4

!p.multi=[0,1,4]
if !d.name eq 'X' then window,0,xsize=1200,ysize=1000
plot,xoff1,area1,xtitle='X',title='center of R2=1 going through center of R1=2',xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,charsize=3
plot,xoff2,area2,xtitle='X',title='center of R2=2 going through edge of R1=1',xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,charsize=3,ytitle='fractional area of R=1 covered'
plot,xoff3,area3,xtitle='X',title='center of R2=1 going through edge of R1=2',xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,charsize=3
plot,xoff4,area4,xtitle='X',title='R2=1 going through just inside of R1=5',xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,charsize=3
!p.multi=0

end
