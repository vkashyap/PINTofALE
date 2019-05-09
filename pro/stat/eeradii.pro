function eeradii,xx,yy,eelev,eree=eree,bkgscal=bkgscal,bkgct=bkgct,$
	cenX=cenX,cenY=cenY,nmin=nmin,verbose=verbose, _extra=e
;+
;function	eeradii
;	computes the enclosed energy radii and corresponding error at the
;	specified levels for a list of events, accounting for background
;	contamination
;
;syntax
;	ree=eeradii(xx,yy,eelev,eree=eree,nmin=nmin,bkgct=bkgct,bkgscal=bkgscal,$
;	    cenX=cenX,cenY=cenY,nmin=nmin,verbose=verbose)
;
;parameters
;	xx	[INPUT; required] X positions of events
;	yy	[INPUT; required] Y positions of events
;	eelev	[INPUT] enclosed energy levels at which to compute radii
;		* if not given, computes the radii corresponding to EE=85%
;
;keywords
;	eree	[OUTPUT] error bars on REE computed assuming a symmetric binomial error
;		(will not work well for EELEV close to 0 or 1, and for small numbers of events)
;	bkgct	[INPUT; default=0] number of counts in the background region
;	bkgscal	[INPUT; default=1] ratio of the background to source areas
;	cenX	[INPUT] if given overrides the central X location determined by centroiding
;	cenY	[INPUT] if given overrides the central Y location determined by centroiding
;	nmin	[INPUT; default=100] minimum number of events before any calculations are done
;		* hardcoded minimum is 50
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2018nov)
;	added keywords CENX,CENY (VK; 2019may)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(xx) & ny=n_elements(yy) & nl=n_elements(eelev)
minph=100L & if keyword_set(nmin) then minph=long(nmin[0])>50L
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X positions of events are not given' else $
  if ny eq 0 then ok='Y positions of events are not given' else $
   if nx ne ny then ok='X and Y positions are incompatible' else $
    if nx lt minph then ok='Too few events, reset NMIN (which cannot be <50)'
if ok ne 'ok' then begin
  print,'Usage: ree=eeradii(xx,yy,eelev,eree=eree,bkgct=bkgct,bkgscal=bkgscal,nmin=nmin,verbose=verbose)'
  print,'  returns EE radii of list of events at specified EE levels, accounting for background'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
bgct=0L & if keyword_set(bkgct) then bgct=long(abs(bkgct[0]))>1L
backscal=1.D & if keyword_set(bkgscal) then backscal=double(bkgscal[0])
clev=0.85 & if nl gt 0 then clev=eelev[*] & nl=n_elements(clev)

;	outputs
ree=fltarr(nl) & eree=ree

;	compute
xcen=mean(xx,/double) & ycen=mean(yy,/double)
if keyword_set(cenX) then xcen=cenX[0]
if keyword_set(cenY) then ycen=cenY[0]
dd=sqrt((xx-xcen)^2+(yy-ycen)^2)
os=sort(dd) & dd=dd[os] & xs=xx[os]-xcen & ys=yy[os]-ycen
areacirc=!dpi*dd^2 & areas=areacirc
go_on=1 & i=nx-1L
while go_on do begin
  j=(i-minph/2)>0L
  xs2=xs[j:i] & ys2=ys[j:i]
  qhull,xs2,ys2,tr & jj=reform(tr[0,*])
  areas[i]=areapoly(xs2[jj],ys2[jj])
  if abs((areas[i]-areacirc[i])/areacirc[i]) lt 0.05 then go_on=0
  if i eq 3 then go_on=0
  i=i-1L
endwhile
maxarea=max(areas)

;	this is the slower, but more exact version
;for i=nx-1L,3L,-1L do begin
  ;j=(i-minph/2)>0L & xs2=xs[j:i] & ys2=ys[j:i]
  ;qhull,xs2,ys2,tr & jj=reform(tr[0,*])
  ;areas[i]=areapoly(xs2[jj],ys2[jj])
  ;;
  ;;xs=xs[0L:i-1L] & ys=ys[0L:i-1L]
;endfor

cdf=dindgen(nx)
cdfb=(cdf/(nx-1L))*(areas/maxarea)*(float(bgct)/backscal)

cdfs=cdf-cdfb & cdfs=cdfs/max(cdfs)
o1=where(cdfs ge 1,mo1) & cdfs[o1[0]:*]=1.

ree=interpol(dd,cdfs,clev)
cdfsig=sqrt(clev*(1.-clev)/(nx+bgct/backscal^2))
uree=interpol(dd,cdfs,(clev+cdfsig)<1)
lree=interpol(dd,cdfs,(clev-cdfsig)>0)
eree=(uree-lree)/2.

if vv gt 1000 then stop,'halting; type .CON to continue'

return,ree
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example case, a Gaussian sitting on a pedestal
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

peasecolr
nsrc=10000L & nbkg=1000L
xs=randomn(seed,nsrc) & ys=randomn(seed,nsrc) & osrc=where(abs(xs) lt 6 and abs(ys) lt 6,mosrc) & xs=xs[osrc] & ys=ys[osrc]
xb=randomu(seed,nbkg)*12.-6. & yb=randomu(seed,nbkg)*12.-6.
xx=[xs,xb] & yy=[ys,yb] & xcen=mean(xx,/double) & ycen=mean(yy,/double) & dd=sqrt((xx-xcen)^2+(yy-ycen)^2) & os=sort(dd)
nph=n_elements(dd)

window,0,ysize=1000
!p.multi=[0,1,2]
plot,xx,yy,psym=3,/xs,/ys,xr=[-6,6],yr=[-6,6] & oplot,xb,yb,psym=4,col=2
plot,dd[os],dindgen(n_elements(dd))/float(nph-1L),/xs,/ys

if not keyword_set(eelev) then eelev=[0.39,0.64,0.85,0.9]
if not keyword_set(verbose) then verbose=1
nl=n_elements(eelev)

ree=eeradii(xx,yy,eelev,eree=eree,bkgct=nbkg,bkgscal=1.,nmin=nmin,verbose=verbose)

for i=0L,nl-1L do begin
  oplot,[!x.crange[0],ree[i]],eelev[i]*[1,1],col=1,line=1,thick=2
  polyfill,ree[i]+eree[i]*[-1,1,1,-1,-1],eelev[i]*[0,0,1,1,0],col=3
endfor
;	NOTE: the overshoot is because background is removed from the original cdf
!p.multi=0

snl=strtrim(nl,2)
print,eelev,ree,eree,form='("EE ="'+snl+'(f5.2,3x),/,"rEE="'+snl+'(f6.3,2x),/," +- " '+snl+'(f7.4,1x))'

end
