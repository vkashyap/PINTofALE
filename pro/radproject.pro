function radproject,lat,lon,inrad=inrad,maxrad=maxrad,rrng=rrng,rdel=rdel,$
	opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0,$
	xout=xout,yout=yout,zout=zout,rproj=rproj,rout=rout,volout=volout,$
	verbose=verbose, _extra=e
;+
;function	radproject
;	Projects a ray originating from an arbitrary point on an arbitrarily
;	oriented sphere onto plane at \infty.  Returns a structure containing
;	the radial coordinate, and the corresponding x-coords, y-,z-offsets,
;	the magnitude of the projection, the corresponding volume of the
;	pillbox, as well as all the input values.
;
;syntax
;	rstr=radproject(lat,lon,inrad=inrad,maxrad=maxrad,rrng=rrng,rdel=rdel,$
;	opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0,$
;	xout=xout,yout=yout,zout=zout,rproj=rproj,rout=rout,volout=volout,$
;	verbose=verbose)
;
;parameters
;	lat	[INPUT; required] latitude on sphere at which ray is anchored
;	lon	[INPUT; required] longitude on sphere at which ray is anchored
;		* LAT,LON must be scalar
;		* (LAT,LON)=(0,0) projects to (Y,Z)=(0,0)
;
;keywords
;	inrad	[INPUT; default=1] minimum r value of ray to consider
;	maxrad	[INPUT; default=10] maximum r value of ray to consider
;	rrng	[INPUT; default=[INRAD,MAXRAD]] filter the output by limiting
;		the radii to consider
;		* minmax(RRNG) is forced to not spill beyond [INRAD,MAXRAD]
;		* use this only if you want everything outside this range to be
;		  explicitly thrown away. discarded.  not zeroed, but removed
;		  from the output.
;	rdel	[INPUT; default=0.1] steps in r
;	opaq	[INPUT; default=0] decides how to treat those portions of rays
;		that are occluded by the sphere defined by radius INRAD
;		- by default, this is assumed to be completely transparent
;		- set OPAQ to negative number to set the opacity as exp(-OPAQ)
;		- any value OPAQ>0 results in a completely opaque inner sphere
;		  (so set /OPAQ to get an opaque inner sphere)
;	arcone	[INPUT; default=1] area of the opening cone at radius=1 [deg^2]
;	theta0	[INPUT; default=0] tilt of sphere [deg]
;	phi0	[INPUT; default=0] twist of sphere [deg]
;		* theta0 and phi0 are defined in a coordinate frame where the
;		  projection plane is at x=\infty.
;		* it is assumed that to get to the inclined sphere from the
;		  observer's frame, you first apply THETA0 and then PHI0.
;	rout	[OUTPUT] radii at which projections are computed
;	xout	[OUTPUT] computed X
;	yout	[OUTPUT] projected Y
;	zout	[OUTPUT] projected Z
;	rproj	[OUTPUT] projection of the radial unit vector towards the Y-Z
;		plane, weighted by the opacity, at each ROUT
;	volout	[OUTPUT] volume of the pillbox at each ROUT
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	.run radproject
;
;history
;	vinay kashyap (2012dec; supersedes pathrofile)
;	added call to SPHTRIG_TRANSLB (VK; 2013feb)
;	bug fixes (theta0,phi0=0) (VK; 2013mar)
;-

;	usage
ok='ok' & np=n_params() & nl=n_elements(lat) & nb=n_elements(lon)
if np lt 2 then ok='Insufficient parameters' else $
 if nl eq 0 then ok='LAT is not defined' else $
  if nb eq 0 then ok='LON is not defined' else $
   if nl gt 1 then ok='LAT must be a scalar' else $
    if nb gt 1 then ok='LON must be a scalar'
if ok ne 'ok' then begin
  print,'Usage: rstr=radproject(lat,lon,inrad=inrad,maxrad=maxrad,rrng=rrng,rdel=rdel,$'
  print,'       opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0,$'
  print,'       xout=xout,yout=yout,zout=zout,rproj=rproj,rout=rout,$
  print,'       verbose=verbose)'
  print,'  projects ray from sphere onto Y-Z plane at \infty'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	figure out input keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
r0=1.D & if keyword_set(inrad) then r0=abs(double(inrad[0]))
r1=10.D & if keyword_set(maxrad) then r1=abs(double(maxrad[0]))
radrng=[r0,r1] & nrrng=n_elements(rrng)
if nrrng gt 0 then begin
  if nrrng eq 1 then begin
    if abs(rrng[0]-r0) lt abs(r1-rrng[0]) then $
      radrng[0]=rrng[0]>r0 else radrng[1]=rrng[0]<r1
  endif else begin
    rr0=min(rrng,max=rr1)
    radrng[0]=rr0>r0
    radrng[1]=rr1<r1
  endelse
endif
dr=0.1D & if keyword_set(rdel) then dr=abs(double(rdel[0]))
ocor=1. & if keyword_set(opaq) then begin
  if opaq[0] lt 0 then ocor=exp(opaq[0])
  if opaq[0] gt 0 then ocor=0.
endif
acone=1.D & if keyword_set(arcone) then acone=abs(double(arcone[0]))
if acone eq 0 then acone=1.D
th0=0. & if keyword_set(theta0) then th0=abs(theta0[0] mod 180)
ph0=0. & if keyword_set(phi0) then ph0=abs(phi0[0] mod 360)

;	figure out Rout
nr=long((r1-r0)/dr+0.5) & rout=dindgen(nr)*dr+r0
if nrrng gt 0 then begin
  ok=where(rout ge radrng[0] and rout le radrng[1],mok)
  if mok eq 0 then begin
    message,'No points in selected radius range: '+$
    	'probably a digitization issue -- expand range',/informational
    return,-1L
  endif
  rout=rout[ok]
endif

;	rotate LAT and LON as seen from projected plane
;
;THESE ARE SMALL ANGLE CORRECTIONS.  FULL SPHERICAL CALC TBD
;lon2=lon-ph0
;lat2=lat-th0
;
tht=90-lat
tbstr=sphtrig_translb(tht,lon,tilt=th0,twist=ph0,_extra=e)
lat2=90.-tbstr.THETA
lon2=tbstr.PHI

;	compute projected (X,Y,Z)
zout=rout*cos((90.D -lat2[0])*(!pi/180.))
xout=rout*sin((90.D -lat2[0])*(!pi/180.))*cos(lon2[0]*(!pi/180.))
yout=rout*sin((90.D -lat2[0])*(!pi/180.))*sin(lon2[0]*(!pi/180.))
roff=sqrt(yout^2+zout^2)
volout=acone*rout^2*dr
rproj=sin((90.-lat2[0])*(!pi/180.))*cos(lon2[0]*(!pi/180.))+0.*rout
o0=where(roff le r0 and rproj le 0,mo0) & if mo0 gt 0 then volout[o0]=volout[o0]*ocor

if vv gt 1000 then stop,'HALTing; type .CON to continue'

parstr=create_struct('LAT',lat,'LON',lon,'INRAD',r0,'MAXRAD',r1,'RDEL',dr,$
	'OPAQ',ocor,'ARCONE',acone,'THETA0',th0,'PHI0',ph0)
outstr=create_struct('R',rout,'X',xout,'Y',yout,'Z',zout,'P',rproj,'VOL',volout,$
	'PAR',parstr)

return,outstr
end

;	example script to illustrate usage

if not keyword_set(npt) then npt=1000L	;number of (LAT,LON) pairs to look at
rb=randomu(seed,npt)*180-90 & rl=randomu(seed,npt)*360
if not keyword_set(inrad) then inrad=4.
if not keyword_set(maxrad) then maxrad=10.
if not keyword_set(rrng) then rrng=[5.,8.]
if not keyword_set(rdel) then rdel=0.01
if not keyword_set(opaq) then opaq=1
if not keyword_set(arcone) then arcone=1.
theta0=0. & phi0=0.

plot,[0],/nodata,xr=maxrad*[-1,1],yr=maxrad*[-1,1],xtitle='Y',ytitle='Z' & peasecolr
for i=0L,npt-1L do begin
  ;rstr=radproject(rb[i],rl[i],inrad=inrad,maxrad=maxrad,rrng=rrng,rdel=rdel,opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0,verbose=verbose)
  rstr=radproject(rb[i],rl[i],inrad=inrad,maxrad=maxrad,rdel=rdel,opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0,verbose=verbose)
  op=where(rstr.P gt 0 and rstr.VOL gt 0,mop) & om=where(rstr.P lt 0 and rstr.VOL gt 0,mom)
  if mop gt 0 then oplot,rstr.Y[op],rstr.Z[op],col=1
  if mom gt 0 then oplot,rstr.Y[om],rstr.Z[om],col=2
  if mop gt 0 then oplot,rstr.Y[op[[0]]],rstr.Z[op[[0]]],col=3,psym=4
  if mom gt 0 then oplot,rstr.Y[om[[0]]],rstr.Z[om[[0]]],col=4,psym=4
  oplot,rstr.Y,rstr.Z,line=2,col=5
endfor

end
