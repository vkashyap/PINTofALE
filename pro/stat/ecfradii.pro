function ecfradii,xx,yy,ecflev,ecfmax=ecfmax,erecf=erecf,bkgscal=bkgscal,bkgct=bkgct,$
	cenX=cenX,cenY=cenY,nmin=nmin,verbose=verbose, _extra=e
;+
;function	ecfradii
;	computes the enclosed energy radii and corresponding error at the
;	specified levels for a list of events, accounting for background
;	contamination
;
;syntax
;	recf=ecfradii(xx,yy,ecflev,ecfmax=ecfmax,erecf=erecf,nmin=nmin,bkgct=bkgct,bkgscal=bkgscal,$
;	    cenX=cenX,cenY=cenY,nmin=nmin,verbose=verbose)
;
;parameters
;	xx	[INPUT; required] X positions of events
;	yy	[INPUT; required] Y positions of events
;	ecflev	[INPUT] enclosed energy levels at which to compute radii
;		* if not given, computes the radii corresponding to ecf=85%
;
;keywords
;	ecfmax	[INPUT] sometimes (xx,yy) does not cover the full set, and only a subset
;		is included.  If given, this factor corrects where the max of the cdf
;		falls.
;		* default (and hardcoded maximum) is 1
;	erecf	[OUTPUT] error bars on Recf computed assuming a symmetric binomial error
;		(will not work well for ecfLEV close to 0 or 1, and for small numbers of events)
;	bkgct	[INPUT; default=0] number of counts in the background region
;	bkgscal	[INPUT; default=1] ratio of the background to source areas
;	cenX	[INPUT] if given overrides the central X location determined by centroiding
;	cenY	[INPUT] if given overrides the central Y location determined by centroiding
;	nmin	[INPUT; default=100] minimum number of events before any calculations are done
;		* hardcoded minimum is 10
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2018nov)
;	added keywords CENX,CENY (VK; 2019may)
;	added keyword EEMAX; changed hardcoded min for NMIN from 50 to 10 (VK; 2019oct)
;	changed name from EERADII to ECFRADII, changed all instances of EE to ECF (VK; 2024jul)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(xx) & ny=n_elements(yy) & nl=n_elements(ecflev)
minph=100L & if keyword_set(nmin) then minph=long(nmin[0])>10L
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X positions of events are not given' else $
  if ny eq 0 then ok='Y positions of events are not given' else $
   if nx ne ny then ok='X and Y positions are incompatible' else $
    if nx lt minph then ok='Too few events, reset NMIN (which cannot be <10)'
if ok ne 'ok' then begin
  print,'Usage: recf=ecfradii(xx,yy,ecflev,ecfmax=ecfmax,erecf=erecf,bkgct=bkgct,bkgscal=bkgscal,nmin=nmin,verbose=verbose)'
  print,'  returns ecf radii of list of events at specified ecf levels, accounting for background'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
bgct=0L & if keyword_set(bkgct) then bgct=long(abs(bkgct[0]))>1L
backscal=1.D & if keyword_set(bkgscal) then backscal=double(bkgscal[0])
clev=0.85 & if nl gt 0 then clev=ecflev[*] & nl=n_elements(clev)
ecfcorr=1.0D & if n_elements(ecfmax) gt 0 then ecfcorr=abs(ecfmax[0])<1.0D

;	outputs
recf=fltarr(nl) & erecf=recf

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
cdfs=cdfs*ecfcorr
o1=where(cdfs ge 1,mo1) & cdfs[o1[0]:*]=1.

recf=interpol(dd,cdfs,clev)
cdfsig=sqrt(clev*(1.-clev)/(nx+bgct/backscal^2))
urecf=interpol(dd,cdfs,(clev+cdfsig)<1)
lrecf=interpol(dd,cdfs,(clev-cdfsig)>0)
erecf=(urecf-lrecf)/2.

if vv gt 1000 then stop,'halting; type .CON to continue'

return,recf
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

if not keyword_set(ecflev) then ecflev=[0.39,0.64,0.85,0.9]
if not keyword_set(verbose) then verbose=1
nl=n_elements(ecflev)

recf=ecfradii(xx,yy,ecflev,erecf=erecf,bkgct=nbkg,bkgscal=1.,nmin=nmin,verbose=verbose)

for i=0L,nl-1L do begin
  oplot,[!x.crange[0],recf[i]],ecflev[i]*[1,1],col=1,line=1,thick=2
  polyfill,recf[i]+erecf[i]*[-1,1,1,-1,-1],ecflev[i]*[0,0,1,1,0],col=3
endfor
;	NOTE: the overshoot is because background is removed from the original cdf
!p.multi=0

snl=strtrim(nl,2)
print,ecflev,recf,erecf,form='("ecf ="'+snl+'(f5.2,3x),/,"recf="'+snl+'(f6.3,2x),/," +- " '+snl+'(f7.4,1x))'

end
