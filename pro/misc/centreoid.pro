function centreoid,xarr,yarr,bkgval,$
	xbound=xbound,ybound=ybound,pbound=pbound,nxgrid=nxgrid,nygrid=nygrid,$
	xcen=xcen,ycen=ycen,xcerr=xcerr,ycerr=ycerr,poserr=poserr,poswdt=poswdt,$
	verbose=verbose, _extra=e
;+
;function	leastasymcentre
;	compute the least asymmetric center (LAC) for a given 2D array, using
;	the algorithm of Lust et al. 2014 (PASP 126, 946, 1092)
;	returns the (x,y) location of the LAC
;
;syntax
;	least_asymmetric_center=centreoid(xarr,yarr,bkgval,$
;	xbound=xbound,ybound=ybound,pbound=pbound,nxgrid=nxgrid,nygrid=nygrid,$
;	xcen=xcen,ycen=ycen,xcerr=xcerr,ycerr=ycerr,poserr=poserr,poswdt=poswdt,$
;	verbose=verbose)
;
;parameters
;	xarr	[INPUT; required] array of X positions
;	yarr	[INPUT; required] array of Y positions
;		* XARR and YARR sizes must match
;	bkgval	[INPUT; default=0] expected background per unit area
;		(see [XY]BIN below for how big area should be)
;
;keywords
;	xbin	[INPUT; default=1] size of a bin along X-axis
;	ybin	[INPUT; default=XBIN] size of a bin along Y-axis
;	xbound	[INPUT] if set, limits the X range over which the LAC is
;		calculated to this
;		* default is to use central PBOUND% of the mass of (XARR,YARR)
;		* if given as scalar, then set to +-[XY]CEN
;	ybound	[INPUT] as XBOUND, for Y axis
;		* default is equal to XBOUND
;	pbound	[INPUT; default=39%] percentage of mass of (XARR,YARR) to use
;		to determine the space over which to check for LAC positions
;		* PBOUND% from centroid position, note
;		* ignored if XBOUND is set
;		* if between -100:-1, counts off from 100 (e.g., -10 ==> include 90%)
;		* if between -1:+1, assumes it is given as fraction
;		* if >100, uses that many elements of (XARR,YARR)
;		* if <-100, excludes abs(PBOUND) many of (XARR,YARR) farthest from (XCEN,YCEN)
;	nxgrid	[INPUT; default=10] the number of grid points used to make the search grid
;		for LAC along the X-axis is 2*NXGRID+1
;	nygrid	[INPUT; default=NXGRID] like NXGRID, for Y-axis
;	xcen	[OUTPUT] simple centroid of XARR
;	ycen	[OUTPUT] simple centroid of YARR
;	xcerr	[OUTPUT] 1-sigma error on XCEN
;	ycerr	[OUTPUT] 1-sigma error on YCEN
;	poserr	[OUTPUT] 1-sigma standard error on LAC position
;	poswdt	[OUTPUT] 1-sigma width of the (XARR,YARR) distribution, computed
;		for all points within the 90% background-corrected EE radius
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (MMXIX.V)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(xarr) & ny=n_elements(yarr)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='XARR is undefined' else $
  if ny eq 0 then ok='YARR is undefined' else $
   if nx lt 2 then ok='XARR has too few elements to bother' else $
    if ny lt 2 then ok='YARR has too few elements to bother' else $
     if nx ne ny then ok='XARR and YARR are incompatible'
if ok ne 'ok' then begin
  print,'Usage: least_asymmetric_center=centreoid(xarr,yarr,bkgval,$'
  print,'       xbound=xbound,ybound=ybound,nxgrid=nxgrid,nygrid=nygrid,$'
  print,'       xcen=xcen,ycen=ycen,xcerr=xcerr,ycerr=ycerr,poserr=poserr,poswdt=poswdt,$'
  print,'       verbose=verbose)
  print,'  compute and return least asymmetric center (Lust et al 2014)'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	first compute the plain centroid
xcen=mean(xarr,/nan,/double)
ycen=mean(yarr,/nan,/double)
;	and errors thereon
xcerr=sqrt(mean(xarr^2,/nan,/double)-xcen^2)/sqrt(nx)
ycerr=sqrt(mean(yarr^2,/nan,/double)-ycen^2)/sqrt(ny)

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
mxgrid=10L & if keyword_set(nxgrid) then mxgrid=nxgrid
mygrid=mxgrid & if keyword_set(nygrid) then mygrid=nygrid
nlim=0 & plim=39.
if keyword_set(pbound) then begin
  if pbound[0] lt 0 then begin
    if pbound[0] lt -100 then nlim=-1 else $
     if pbound[0] ge -100 and pbound[0] lt -1 then plim=100+pbound[0] else $
      if pbound[0] ge -1 then plim=100.*(1.+pbound[0])
  endif else begin
    if pbound[0] gt 0 and pbound[0] lt 1 then plim=100.*pbound[0] else $
     if pbound[0] ge 1 and pbound[0] lt 100 then plim=pbound[0] else $
      if pbound[0] ge 100 then nlim=1
  endelse
endif
xlim=0 & if keyword_set(xbound) then xlim=xbound
ylim=xlim & if keyword_set(ybound) then ylim=ybound
;
if keyword_set(xlim) and n_elements(xlim) ne 2 then xlim=xcen+xbound[0]*[-1,1]
if keyword_set(ylim) and n_elements(ylim) ne 2 then ylim=ycen+ybound[0]*[-1,1]
;
if not keyword_set(xlim) or not keyword_set(ylim) then begin
  ;	compute cdf to get a sense of the range
  dd=sqrt((xarr-xcen)^2+(yarr-ycen)^2) & os=sort(dd) & dd=dd[os] & cct=(dindgen(nx)+1.)	;total counts within DD
  cbg=bkgval*!pi*dd^2	;expected background counts within DD
  cdf=(cct-cbg)>0	;estimated source counts within DD
  ;	the >0 to avoid getting into numerical difficulties with voids near the centroid
  cdf=cdf/max(cdf) & o1=where(cdf eq 1) & cdf[o1[0]:*]=1.
  if nlim eq 0 then dlim=interpol(dd,cdf,plim/100.)>0 else begin
    if nlim lt 0 then mlim=(((nx-1L)-abs(pbound[0]))>0L)<(nx-1L) else mlim=((pbound[0])>0L)<(nx-1L)
    dlim=dd[mlim]
  endelse
  if dlim eq 0 then dlim=stddev([xarr-xcen,yarr-ycen])
  if not keyword_set(xlim) then xlim=xcen+dlim*[-1,1]
  if not keyword_set(ylim) then ylim=ycen+dlim*[-1,1]
endif
nxlim=n_elements(xlim) & nylim=n_elements(ylim)
if nxlim ne 2 or nylim ne 2 then message,'BUG!'

;	compute variance of radial profiles for each grid point
nxbin=2L*mxgrid+1L & nybin=2L*mygrid+1L
binx=(xlim[1]-xlim[0])/(nxbin-1L) & biny=(ylim[1]-ylim[0])/(nybin-1L)
var=dblarr(nxbin,nybin) & xvar=fltarr(nxbin) & yvar=fltarr(nybin)
zr=shift(dist(4L*nxbin+1L,4L*nybin+1L),2*nxbin,2*nybin) & zur=zr[uniq(zr,sort(zr))] & nzur=n_elements(zur)
for i=0L,nxbin-1L do begin
  if vv gt 0 then kilroy,dot=strtrim(i,2)+'.. '
  xvar[i]=i*binx+xlim[0]
  for j=0L,nybin-1L do begin
    if vv gt 0 then kilroy,dot=strtrim(j,2)+' '
    yvar[j]=j*biny+ylim[0]
    zx=(xarr-xvar[i])/binx & zy=(yarr-yvar[j])/biny
    img=hist_2d(zx,zy,min1=-2*nxbin,min2=-2*nybin,max1=2*nxbin,max2=2*nybin,bin1=1,bin2=1)
    for k=0L,nzur-1L do begin
      oz=where(zr eq zur[k],moz)
      if moz gt 2 then var[i,j]=var[i,j]+variance(img[oz])*moz
    endfor
    if vv gt 5 then print,xvar[i],yvar[j],var[i,j]
  endfor
endfor

;	quick first order estimate until interpolation is implemented
jnk=min(var,imn) & ii=array_indices(var,imn) & pos0=[xvar[ii[0]],yvar[ii[1]]]
;	including quadratic interpolation
ipos=gcentroid(max(var)-var)	;returns the indices corresponding to max of input array, inverted here so we pick minimum of variance
xpos=interpol(xvar,findgen(nxbin),ipos[0])
ypos=interpol(yvar,findgen(nybin),ipos[1])
pos=[xpos,ypos]

;	now compute stderr simply as sqrt of mean of squared residuals, divided by sqrt of number of events
dd=sqrt((xarr-xpos)^2+(yarr-ypos)^2) & dd2=dd & os=sort(dd) & dd=dd[os] & cct=(dindgen(nx)+1.)
cbg=bkgval*!pi*dd^2
cdf=(cct-cbg)>0
cdf=cdf/max(cdf) & o1=where(cdf eq 1) & cdf[o1[0]:*]=1.
dlim=interpol(dd,cdf,0.9) & ok=where(dd2 le dlim,mok)
xposerr=sqrt(total((xarr[ok]-xpos)^2)/mok)/sqrt(mok)
yposerr=sqrt(total((yarr[ok]-ypos)^2)/mok)/sqrt(mok)
poserr=[xposerr,yposerr]	;std error on LAC position
poswdt=poserr*sqrt(mok)		;width of PSF

if vv gt 1000 then stop,'halting; type .CON to continue'

return,pos
end
