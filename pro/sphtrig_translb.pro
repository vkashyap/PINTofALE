function sphtrig_translb,theta,phi,tilt=tilt,twist=twist,eps=eps,$
	verbose=verbose, _extra=e
;+
;function	sphtrig_translb
;	transform the given {theta,phi} coordinates on a sphere that
;	has been tilted and twisted (in that order) to {theta',phi'}
;	on an untwisted and untilted pov reference frame, and
;	returns them in a structure.
;
;syntax
;	tpstr=sphtrig_translb(theta,phi,tilt=tilt,twist=twist,$
;	eps=eps,verbose=vv)
;
;parameters
;	theta	[INPUT; required] angular coordinate from Z-axis [deg]
;	phi	[INPUT; required] azimuthal coordinate around Z-axis [deg]
;		* THETA,PHI may be vectors or scalars
;		* if vectors, assumed to be paired sets, with
;		  values filled in with 0's if one falls short
;
;keywords
;	tilt	[INPUT; default=0] inclination angle of sphere
;		relative to pov coordinate frame [deg]
;	twist	[INPUT; default=0] rotation angle by which the
;		vector going through the pole is rotated around
;		the Z-axis of the pov reference frame [deg]
;		* note that TWIST is given in the negative sense --
;		  the angle by which the frame must be rotated to
;		  return to 0.  the program takes the -ve of the
;		  supplied TWIST.
;	eps	[INPUT; default=1e-6] a small number
;		* There are some singularities when THETA and PHI are zero.
;		  This is added on to make sure at least something
;		  close to the limiting value is returned in that case
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	.run sphtrig_translb
;	makes lines of constant latitude, constant longitude, and
;	draws the figure used to derive the transformation used here.
;	(the ultimate in self-documenting!)
;
;subroutines
;	example run requires PEASECOLR
;
;history
;	vinay kashyap (2013feb)
;	added keyword EPS (2013mar)
;	corrected sign error that was making PHI go clockwise instead of
;	  counter-clockwise (2015sep)
;-

;	usage
ok='ok' & npar=n_params() & nt=n_elements(theta) & np=n_elements(phi)
if npar lt 2 then ok='Insufficient parameters' else $
 if nt eq 0 then ok='THETA not defined' else $
  if np eq 0 then ok='PHI not defined'
if ok ne 'ok' then begin
  print,'Usage: tpstr=sphtrig_translb(theta,phi,tilt=tilt,twist=twist,eps=eps,verbose=vv)'
  print,'  compute theta,phi transformed from tilted and twisted sphere'
  print,'  to line up with point-of-view coordinate frame'
  if npar gt 0 then message,ok,/informational
  return,-1L
endif

;	parse input
nn=nt > np
zth=dblarr(nn) & zth[0L:np-1L]=theta[*]
zph=dblarr(nn) & zph[0L:nt-1L]=phi[*]

;	parse keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
if not keyword_set(eps) then eps=1d-6
tht0=0.D + eps & if keyword_set(tilt) then tht0=tilt[0]
while tht0 lt 0 do tht0=tht0+180.D
if tht0 gt 180. then tht0=(tht0 mod 180.D)
phi0=0.D + eps & if keyword_set(twist) then phi0=twist[0]
while phi0 lt 0 do phi0=phi0+360.D
if phi0 gt 360. then phi0=(phi0 mod 360.D)
phi0=-phi0

;	convert to radians
zth=zth*!dpi/180.D & zph=zph*!dpi/180.D & tht0=tht0*!dpi/180.D & phi0=phi0*!dpi/180.D

;	points, arcs, and angles refer to the documentation figure, which
;	can be viewed by running the example script to drive this program.

;	first, get angle CSB
;arcBC=2.*!DPI-zph	;this changes sign of phi
arcBC=zph
arcSC=!DPI/2.
arcSB=!DPI/2.
;	cos(arcBC) = cos(arcSC)*cos(arcSB) + sin(arcSC)*sin(arcSB)*cos(angCSB)
;	==> cos(arcBC) = 0*0 + 1*1*cos(angCSB)
;	==> angCSB = arcBC
angCSB=arcBC

;	now, get arc ZX
arcSX=zth
arcZS=tht0
angZSX=!pi-angCSB
;	cos(arcZX) = cos(arcZS)*cos(arcSX) + sin(arcZS)*sin(arcSX)*cos(angZSX)
arcZX = acos(cos(arcZS)*cos(arcSX) + sin(arcZS)*sin(arcSX)*cos(angZSX))
sgnZSO = intarr(nn)+1 & o0=where(angZSX lt 0,mo0) & if mo0 gt 0 then sgnZSO[o0]=-1

;	get angle XZS
;	cos(arcSX) = cos(arcZS)*cos(arcZX) + sin(arcZS)*sin(arcZX)*cos(angXZS)
;	cos(angXZS) = (cos(arcSX)-cos(arcZS)*cos(arcZX))/(sin(arcZS)*sin(arcZX))
angXZS = sgnZSO*acos((cos(arcSX)-cos(arcZS)*cos(arcZX))/(sin(arcZS)*sin(arcZX)))

;	angle UZV
arcZU=!DPI/2.
arcZV=!DPI/2.
arcUV=phi0
;	cos(arcUV) = cos(arcZU)*cos(arcZV) + sin(arcZU)*sin(arcZV)*cos(angUZV)
;	==> cos(arcUV) = 0*0 + 1*1*cos(angUZV)
;	==> angUZV = arcUV
angUZV=arcUV

;	angle WZU
angWZU = angXZS - angUZV

;	arc UW
arcZW=!DPI/2.
;	cos(arcUW) = cos(arcZW)*cos(arcZU) + sin(arcZW)*sin(arcZU)*cos(angWZU)
;	==> cos(arcUW) = 0*0 + 1*1*cos(angWZU)
arcUW=angWZU

;	output
;	convert to [deg]
tpstr=create_struct('THETA',arcZX*180./!pi,'PHI',arcUW*180./!pi)

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,tpstr
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example
;	set variables TILT, TWIST, STHETA, SPHI
;
;	made slightly more verbose (VK; 2015sep)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

peasecolr & loadct,3 & peasecolr

;	calling sequence
print,'' & print,'CALLING SEQUENCE:' & print,''
jnk=sphtrig_translb()
print,''

;	initialize
if n_elements(tilt) ne 1 then tilt=35. else tilt=tilt[0]
if n_elements(twist) ne 1 then twist=50. else twist=twist[0]
if not keyword_set(verbose) then verbose=1

print,'SET THE INCLINATION OF THE STAR USING VARIABLES --'
help,tilt,twist

if tilt lt 0 then begin
  message,'what does negative TILT mean anyway?',/informational
  tilt=abs(tilt)
  help,tilt,twist
endif
if tilt gt 180 then begin
  message,'nothing can TILT beyond a vertical dive',/informational
  tilt=(tilt mod 180.)
  help,tilt,twist
endif

wait,1
print,'FIRST DISPLAY A LATITUDE GRID'
if !d.name eq 'X' then window,0,xsize=800,ysize=800

;	set up observer's frame
surface,fltarr(2,2),/nodata,xrange=1.2*[-1,1],yrange=1.2*[-1,1],zrange=1.2*[-1,1],xstyle=5,ystyle=5,zstyle=5,ax=20,az=-100,/save
plots,[0,1],[0,0],[0,0],col=3,/t3d,line=1
plots,[0,0],[0,1],[0,0],col=3,/t3d,line=1
plots,[0,0],[0,0],[0,1],col=3,/t3d,line=1

;	now set up the source sphere frame
plots,[0,sin(tilt*!pi/180.)*cos(twist*!pi/180.)],[0,sin(tilt*!pi/180.)*sin(twist*!pi/180.)],[0,cos(tilt*!pi/180.)],col=1,/t3d,thick=2

;	draw latitudes
ang2=findgen(361)
for ilat=0,180,10 do begin
  ang1=0*ang2+ilat
  tpstr=sphtrig_translb(ang1,ang2,tilt=tilt,twist=twist)
  zz=cos(tpstr.THETA*!pi/180.) & xx=sin(tpstr.THETA*!pi/180.)*cos(tpstr.PHI*!pi/180.) & yy=sin(tpstr.THETA*!pi/180.)*sin(tpstr.PHI*!pi/180.)
  plots,xx,yy,zz,/t3d,psym=3,col=1
  ox=where(xx gt 0,mox) & if mox gt 0 then plots,xx[ox],yy[ox],zz[ox],/t3d,psym=4,col=1
  if ilat eq 90 then plots,xx,yy,zz,/t3d,col=1,thick=2
endfor

wait,1
print,'NEXT DISPLAY A LONGITUDE GRID'
if !d.name eq 'X' then window,2,xsize=800,ysize=800

;	set up observer's frame
surface,fltarr(2,2),/nodata,xrange=1.2*[-1,1],yrange=1.2*[-1,1],zrange=1.2*[-1,1],xstyle=5,ystyle=5,zstyle=5,ax=20,az=-100,/save
plots,[0,1],[0,0],[0,0],col=3,/t3d,line=1
plots,[0,0],[0,1],[0,0],col=3,/t3d,line=1
plots,[0,0],[0,0],[0,1],col=3,/t3d,line=1

;	now set up the source sphere frame
plots,[0,sin(tilt*!pi/180.)*cos(twist*!pi/180.)],[0,sin(tilt*!pi/180.)*sin(twist*!pi/180.)],[0,cos(tilt*!pi/180.)],col=1,/t3d,thick=2

;	draw longitudes
ang1=findgen(181)
for ilon=0,340,20 do begin
  ang2=0*ang1+ilon
  tpstr=sphtrig_translb(ang1,ang2,tilt=tilt,twist=twist)
  zz=cos(tpstr.THETA*!pi/180.) & xx=sin(tpstr.THETA*!pi/180.)*cos(tpstr.PHI*!pi/180.) & yy=sin(tpstr.THETA*!pi/180.)*sin(tpstr.PHI*!pi/180.)
  plots,xx,yy,zz,/t3d,psym=3,col=1
  ox=where(xx gt 0,mox) & if mox gt 0 then plots,xx[ox],yy[ox],zz[ox],/t3d,psym=4,col=1
  if ilon eq 0 then plots,xx,yy,zz,/t3d,col=1,thick=2
endfor

;	pick a point, any point, on the source sphere
if n_elements(stheta) ne 1 then stheta=44.
if n_elements(sphi) ne 1 then sphi=111.
while stheta lt 0 do stheta=stheta+180.
if stheta gt 180. then stheta=(stheta mod 180.)
while sphi lt 0 do sphi=sphi+360.
if sphi gt 360. then sphi=(sphi mod 360.)

wait,1
print,'SELECT AN ARBITRARY POINT ON THE STAR USING VARIABLES --'
help,stheta,sphi

print,'SHOW WHERE THE POINT LIES, ALONG WITH EXPLICATORY SPHERICAL TRIANGLES --'
if !d.name eq 'X' then window,1,xsize=800,ysize=800

;	set up the observer's frame
surface,fltarr(2,2),/nodata,xrange=1.2*[-1,1],yrange=1.2*[-1,1],zrange=1.2*[-1,1],xstyle=5,ystyle=5,zstyle=5,ax=20,az=-100,/save
plots,[0,1],[0,0],[0,0],col=3,/t3d,line=1
plots,[0,0],[0,1],[0,0],col=3,/t3d,line=1
plots,[0,0],[0,0],[0,1],col=3,/t3d,line=1
xyouts,0.05,0.05,z=0.05,text_axes=2,charthick=2,charsize=2,col=3,/t3d,'O'
xyouts,1.05,0,z=-0.05,text_axes=2,charthick=2,charsize=2,col=3,/t3d,'U'
xyouts,0.05,0,z=1.05,text_axes=2,charthick=2,charsize=2,col=3,/t3d,'Z'

;	draw some orienting circles
ang=findgen(360)*!pi/180.
xx=cos(ang) & yy=sin(ang) & zz=0*ang & plots,xx,yy,zz,col=3,line=1,/t3d
xx=0*ang & yy=sin(ang) & zz=cos(ang) & plots,xx,yy,zz,col=3,line=1,/t3d

;	mark the tilt
ang=[findgen(fix(tilt)+1),tilt]*!pi/180. & xx=sin(ang)*cos(twist*!pi/180.) & yy=sin(ang)*sin(twist*!pi/180.) & zz=cos(ang)
plots,xx,yy,zz,col=3,/t3d
nang=n_elements(ang) & xyouts,xx[nang-1],yy[nang-1],z=zz[nang-1]+0.05,text_axes=2,charthick=2,charsize=2,col=1,/t3d,'S'

;	now set up the source sphere frame
plots,[0,sin(tilt*!pi/180.)*cos(twist*!pi/180.)],[0,sin(tilt*!pi/180.)*sin(twist*!pi/180.)],[0,cos(tilt*!pi/180.)],col=1,/t3d

;	draw some orienting circles in source frame
ang2=findgen(360) & ang1=0*ang2+90
tpstr=sphtrig_translb(ang1,ang2,tilt=tilt,twist=twist,verbose=verbose)
zz=cos(tpstr.THETA*!pi/180.) & xx=sin(tpstr.THETA*!pi/180.)*cos(tpstr.PHI*!pi/180.) & yy=sin(tpstr.THETA*!pi/180.)*sin(tpstr.PHI*!pi/180.)
plots,xx,yy,zz,psym=3,col=1,/t3d

ang1=findgen(91) & ang2=0*ang1
tpstr=sphtrig_translb(ang1,ang2,tilt=tilt,twist=twist,verbose=verbose)
zz=cos(tpstr.THETA*!pi/180.) & xx=sin(tpstr.THETA*!pi/180.)*cos(tpstr.PHI*!pi/180.) & yy=sin(tpstr.THETA*!pi/180.)*sin(tpstr.PHI*!pi/180.)
plots,xx,yy,zz,line=2,col=1,/t3d
nang=n_elements(ang1) & xyouts,xx[nang-1]-0.05,yy[nang-1]-0.05,z=zz[nang-1]-0.1,text_axes=2,charthick=2,charsize=2,col=1,/t3d,'B'
jnk=min(abs(zz),imn) & plots,[0,xx[imn]],[0,yy[imn]],[0,zz[imn]],line=1,col=1,/t3d
xyouts,xx[imn]-0.05,yy[imn]-0.05,z=zz[imn]-0.1,text_axes=2,charthick=2,charsize=2,col=1,/t3d,'V'
plots,[0,xx[nang-1]],[0,yy[nang-1]],[0,zz[nang-1]],line=1,col=2,/t3d

;	mark the point
slb=sphtrig_translb(stheta,sphi,tilt=tilt,twist=twist,verbose=verbose)
szz=cos(slb.THETA*!pi/180.) & sxx=sin(slb.THETA*!pi/180.)*cos(slb.PHI*!pi/180.) & syy=sin(slb.THETA*!pi/180.)*sin(slb.PHI*!pi/180.)
plots,[0,sxx],[0,syy],[0,szz],col=2,/t3d,symsize=2,thick=2
xyouts,sxx-0.05,syy-0.05,z=szz+0.05,text_axes=2,charthick=2,charsize=2,col=2,/t3d,'X'
print,'Cartesian (x,y,z) coordinates of X in upright frame:',sxx[0],syy[0],szz[0]
print,'Polar (theta,phi) coordinates of X in upright frame:',(slb.THETA)[0],(slb.PHI)[0]

;	mark the triangle in the source frame
ang1=findgen(91) & ang2=0*ang1+sphi
tpstr=sphtrig_translb(ang1,ang2,tilt=tilt,twist=twist,verbose=verbose)
zz=cos(tpstr.THETA*!pi/180.) & xx=sin(tpstr.THETA*!pi/180.)*cos(tpstr.PHI*!pi/180.) & yy=sin(tpstr.THETA*!pi/180.)*sin(tpstr.PHI*!pi/180.)
plots,xx,yy,zz,line=2,col=1,/t3d
nang=n_elements(ang1) & xyouts,xx[nang-1]-0.05,yy[nang-1]-0.05,z=zz[nang-1]-0.1,text_axes=2,charthick=2,charsize=2,col=1,/t3d,'C'
plots,[0,xx[nang-1]],[0,yy[nang-1]],[0,zz[nang-1]],line=1,col=2,/t3d
;
ang1=[findgen(fix(stheta)+1),stheta] & ang2=0*ang1+sphi
tpstr=sphtrig_translb(ang1,ang2,tilt=tilt,twist=twist,verbose=verbose)
zz=cos(tpstr.THETA*!pi/180.) & xx=sin(tpstr.THETA*!pi/180.)*cos(tpstr.PHI*!pi/180.) & yy=sin(tpstr.THETA*!pi/180.)*sin(tpstr.PHI*!pi/180.)
plots,xx,yy,zz,col=1,/t3d
;
ang2=[findgen(fix(sphi)+1),sphi] & ang1=0*ang2+90.
tpstr=sphtrig_translb(ang1,ang2,tilt=tilt,twist=twist,verbose=verbose)
zz=cos(tpstr.THETA*!pi/180.) & xx=sin(tpstr.THETA*!pi/180.)*cos(tpstr.PHI*!pi/180.) & yy=sin(tpstr.THETA*!pi/180.)*sin(tpstr.PHI*!pi/180.)
plots,xx,yy,zz,line=2,col=1,/t3d

;	mark the triangle in the observer's frame
ang1=findgen(91) & ang2=0*ang1+(slb.PHI)[0]
zz=cos(ang1*!pi/180.) & xx=sin(ang1*!pi/180.)*cos(ang2*!pi/180.) & yy=sin(ang1*!pi/180.)*sin(ang2*!pi/180.)
plots,xx,yy,zz,line=2,col=4,/t3d
nang=n_elements(ang1) & xyouts,xx[nang-1]-0.05,yy[nang-1]-0.05,z=zz[nang-1]-0.1,text_axes=2,charthick=2,charsize=2,col=4,/t3d,'W'
plots,[0,xx[nang-1]],[0,yy[nang-1]],[0,0],line=1,col=4,/t3d
;
ang1=[findgen(fix((slb.THETA)[0])+1),(slb.THETA)[0]] & ang2=0*ang1+(slb.PHI)[0]
zz=cos(ang1*!pi/180.) & xx=sin(ang1*!pi/180.)*cos(ang2*!pi/180.) & yy=sin(ang1*!pi/180.)*sin(ang2*!pi/180.)
plots,xx,yy,zz,col=4,/t3d
;
ang1=findgen(91) & ang2=0*ang1
zz=cos(ang1*!pi/180.) & xx=sin(ang1*!pi/180.)*cos(ang2*!pi/180.) & yy=sin(ang1*!pi/180.)*sin(ang2*!pi/180.)
plots,xx,yy,zz,line=2,col=4,/t3d

;	some descriptive labels
xyouts,0.05,0.95-0.*0.03,/norm,charthick=1.5,charsize=2,'ZS = TILT = '+string(tilt,'(f7.2)')
xyouts,0.05,0.95-1.*0.03,/norm,charthick=1.5,charsize=2,'UOV = TWIST = '+string(twist,'(f7.2)')
xyouts,0.05,0.95-2.*0.03,/norm,charthick=1.5,charsize=2,'SX = STHETA = '+string(stheta,'(f7.2)')
xyouts,0.05,0.95-3.*0.03,/norm,charthick=1.5,charsize=2,'COB = SPHI = '+string(sphi,'(f7.2)')
xyouts,0.05,0.95-5.*0.03,/norm,charthick=1.5,charsize=2,'WX = LAT = '+string(90.-(slb.theta)[0],'(f7.2)')
xyouts,0.05,0.95-6.*0.03,/norm,charthick=1.5,charsize=2,'UOW = LON = '+string((slb.phi)[0],'(f7.2)')

end
