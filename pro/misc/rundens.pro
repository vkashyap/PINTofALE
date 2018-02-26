;+
;function	rundens
;	compute and return the running density of a set of measurements.
;
;syntax
;	yrun=rundens(xpt,hwidth,/usex,/optwdt,outhwdt=outhwdt,$
;	yerr=yerr,xerr=xerr,nsim=nsim,ysim=ysim,yesim=yesim,xtol=xtol,$
;	wts=wts,dwts=dwts,lbound=lbound,rbound=rbound,cvmse=cvmse,$
;	verbose=verbose)
;
;parameters
;	xpt	[INPUT; required] array of points for which the running density
;		must be calculated
;		* a running density is calculated at each value of XPT
;	hwidth	[INPUT] the half-width on either side of a point to include in
;		the density
;		* this is assumed to be the number of array indices
;		  unless USEX is set
;		* default is to use 5% of n_elements(XPT)
;		  or 5% of the range of XPT if USEX is set
;		* if OPTWDT is set, computes the optimum width to use
;		  and treats this as just a default
;		* actual value used is returned via keyword OUTHWDT
;
;keywords
;	usex	[INPUT] if set, assumes HWIDTH is in same units as XPT
;	optwdt	[INPUT] if set, figures out the optimum HWIDTH to use by
;		comparing the width of the distribution of the running
;		density estimates to the expected statistical stddev for
;		a flat density
;		* if this is set, HWIDTH will be overridden
;		* the optimum HWIDTH is calculated based on the array
;		  indices -- i.e., USEX is ignored during the calculation,
;		  but if set, is used afterwards to scale the optimum
;		  HWIDTH to the units in XPT.  So beware that it may not
;		  really be the optimum in this case.
;		* if this module fails, the smallest width checked will
;		  be returned
;	outhwdt	[OUTPUT] will contain the value of HWIDTH actually used in
;		the calculation
;	yerr	[OUTPUT] error bar on the density at each point
;		* WARNING: adjacent values of YERR are not independent.
;		  while the error estimate at each point is fine, to
;		  get a truly independent estimate, one has to skip by
;		  2*HWIDTH
;	xerr	[INPUT] error bars on XPT
;		* used to compute a Monte Carlo uncertainty band on the
;		  output running density
;	nsim	[INPUT] number of simulations to carry out
;		* default is none
;	ysim	[OUTPUT] a NSIMxN(XPT) array to contain the running
;		density for each simulation
;	yesim	[OUTPUT] point-wise stddev based on YSIM
;	xtol	[INPUT] the tolerance for considering two points to be
;		identical
;		* essentially, if any two adjacent points are found to
;		  be identical, then all points are jittered by +-0.5*XTOL
;	wts	[INPUT] if given, computes the total of the WTS and returns
;		that rather than the running density
;		* MUST match the size of XPT, otherwise ignored
;	dwts	[INPUT] same as WTS, but the total is computed and placed
;		in the denominator 
;		* MUST match the size of XPT, otherwise ignored
;		* if both WTS and DWTS are given, the ratio of the totals
;		  is computed, not the average of the ratios
;	lbound	[OUTPUT] array indices of leftside bounds for all XPT
;	rbound	[OUTPUT] array indices of rightside bounds for all XPT
;	cvmse	[OUTPUT] cross-validation mean square error between estimate
;		and leave-one-out estimate.  the extra calculations are
;		performed only if this keyword is present and are not done
;		for the simulated data.
;	verbose	[INPUT] controls chatter
;	
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	RUNDENS_OPTWDT (included here)
;	CURVESECT
;	EQT_INTERVAL
;
;history
;	vinay kashyap (MMXI.VII)
;	added keywords WTS and DWTS and corrected bug (VK; MMXI.VIII)
;	added keywords LBOUND, RBOUND, and CVMSE (VK; MMXII.V)
;	added keyword OPTWDT, added function RUNDENS_OPTWDT, and added
;	  calls to RUNDENS_OPTWDT(), CURVESECT(), and EQT_INTERVAL()
;	  (VK; MMXII.XII)
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function rundens_optwdt,xpt,verbose=verbose, _extra=e
;{
;function	rundens_optwdt
;	repeatedly calls rundens() with different values of HWIDTH
;	to find the optimal HWIDTH and returns that value
;
;	the optimum is defined as that value of HWIDTH where the
;	scatter in the running means is comparable to the expected
;	scatter for the given number of points in a flat density curve.
;
;dedicated function, expected to called only from RUNDENS()
;
;vinay kashyap (MMXII.XII)
;}

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L

nph=n_elements(xpt)
if nph eq 0 then begin
  message,'XPT is not defined; quitting',/informational & return,-1L
endif
if nph lt 3 then begin
  message,'XPT is too small; quitting',/informational & return,-1L
endif
if nph lt 16 then begin
  message,'XPT too small; half the distance to the goal line',/informational
  return,nph/2
endif

kph2=fix(alog(nph)/alog(2))
kk=round(sqrt(2)^(findgen(kph2*2-6)+6))
nkk=n_elements(kk)
dfvar=fltarr(nkk) & sfvar=fltarr(nkk)
for k=0L,nkk-1L do begin
  kph=kk[k]
  f=rundens(xpt,kph, _extra=e)
  df=f[1:*]-f
  ;dfvar[k]=stddev(df)	;this is too conservative
  tmp=eqt_interval(df,/fsample,clev=0.68) & dfvar[k]=(tmp[1]-tmp[0])/2.
  sfvar[k]=stddev(f)/sqrt(nph)
endfor

xmax=max(xpt,min=xmin) & dx=xmax-xmin
tmp=sqrt(2*kk+1L)/dx & tmp=sqrt(tmp^2+sfvar^2)
;xy=curvesect(dfvar,sqrt(2*kk+1L)/dx,kk)
xy=curvesect(dfvar,tmp,kk)
if finite(xy[1]) eq 0 then begin
  if vv gt 0 then message,'No intersections found; trying conservative match',$
  	/informational
  tmp=sqrt(2*kk+1L)/dx
  xy=curvesect(dfvar,tmp,kk)
  if finite(xy[1]) eq 0 then begin
    message,'No intersections found; choosing the smallest value',$
    	/informational
    xy[0]=min(kk)
  endif
endif
hw=long(xy[0])

if vv gt 2000 then stop,'HALTing; type .CON to continue'

return,hw
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function rundens,xpt,hwidth,usex=usex,optwdt=optwdt,outhwdt=outhwdt,$
	yerr=yerr,xerr=xerr,nsim=nsim,ysim=ysim,yesim=yesim,xtol=xtol,$
	wts=wts,dwts=dwts,lbound=lbound,rbound=rbound,cvmse=cvmse,$
	verbose=verbose, _extra=e

;	usage
ok='ok' & np=n_params() & nx=n_elements(xpt) & nw=n_elements(hwidth)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='XPT is undefined' else $
  if nx eq 1 then ok='XPT is not an array'
if ok ne 'ok' then begin
  print,'Usage: yrun=rundens(xpt,hwidth,/usex,/optwdt,outhwdt=outhwdt,$
  print,'       yerr=yerr,xerr=xerr,nsim=nsim,ysim=ysim,yesim=yesim,xtol=xtol,$
  print,'       wts=wts,dwts=dwts,lbound=lbound,rbound=rbound,cvmse=cvmse,$
  print,'       verbose=verbose)
  print,'  compute and return running density'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
xmin=min(xpt,max=xmax)
hw=long(0.05*nx+0.5) & if keyword_set(usex) then hw=0.05*(xmax-xmin)
if nw ne 0 then hw=1.*abs(hwidth[0]) & if not keyword_set(usex) then hw=long(hw)
nw=n_elements(hw)
;
msim=0L & if keyword_set(nsim) then msim=long(abs(nsim[0]))>1
;
xsig=fltarr(nx) & nxe=n_elements(xerr)
if nxe gt 0 then xsig[0L:nxe<(nx-1L)]=abs(xerr[0L:nxe<(nx-1L)])
;
tolx=1d-10 & if keyword_set(xtol) then begin
  if xtol[0] ne 1 then tolx=abs(xtol[0])
endif
;
numwt=fltarr(nw)+1. & denwt=numwt
nnwt=n_elements(wts) & if nnwt ne nx then nnwt=0L else numwt=1.0*wts
ndwt=n_elements(dwts) & if ndwt ne nx then ndwt=0L else denwt=1.0*dwts
;
if keyword_set(optwdt) then begin
  ;opthw=rundens_optwdt(xpt,xtol=tolx,wts=numwt,dwts=denwt,verbose=vv) ;UNCOMMENT THIS ONLY AFTER IMPLEMENTING YERR FUNCTIONALITY IN RUNDENS_OPT
  opthw=rundens_optwdt(xpt,xtol=tolx,verbose=vv)
  if opthw[0] ne -1 then begin
    if keyword_set(usex) then hw=opthw*(xmax-xmin)/float(nx) else hw=opthw
    if vv gt 5 then message,$
    	'using HWIDTH='+strtrim(hwidth,2),/informational
  endif else if vv gt 0 then message,$
    	'OPTWDT did not work; continuing with inputs and/or defaults',$
	/informational
endif

;	outputs
outhwdt=hw
yrun=fltarr(nx)*xpt[0] & yerr=yrun
if msim gt 0 then begin
  ysim=xpt[0]*fltarr(msim,nx) & yesim=xpt[0]*fltarr(nx)
endif
lbound=lonarr(nx)-1L & rbound=lbound
if arg_present(cvmse) then begin
  crossval=1
  yrun2=fltarr(nx)*xpt[0] & yerr2=yrun2
endif else crossval=0

;	step through the array and compute the running density
for isim=0L,msim do begin	;{isim=0 is standard run, rest are sims

  xx=xpt & if isim gt 0 then xx=xx+randomn(seed,nx)*xsig
  zn=numwt & zd=denwt
  os=sort(xx) & xx=xx[os] & zn=zn[os] & zd=zd[os]
  delx=min(xx[1:*]-xx)
  if delx eq 0 then begin
    xx=xx+randomu(seed)*(tolx/2.)-tolx/2.
    os=sort(xx) & xx=xx[os] & zn=zn[os] & zd=zd[os]
  endif
  if isim eq ((100-vv)>1)*long(isim/((100-vv)>1)) then kilroy

  for i=0L,nx-1L do begin		;{step through the input array

    if keyword_set(usex) then begin	;(find array range over which to compute density
      ok=where(abs(xx-xx[i]) le hw,mok)
      i0=ok[0] & i1=ok[mok-1L]
      if mok gt 0 then begin
        yy=float(mok)
        if isim eq 0 then ye=sqrt(mok) & yne=0 & yde=0
        ;
        if nnwt gt 0 then begin & yn=total(zn[ok]) & if mok gt 1 and isim eq 0 then yne=stddev(zn[ok]) & endif
        if ndwt gt 0 then begin & yd=total(zd[ok]) & if mok gt 1 and isim eq 0 then yde=stddev(zd[ok]) & endif
        if nnwt gt 0 and ndwt eq 0 then begin
          yy=yn & if isim eq 0 then ye=sqrt(yne^2+mok*(yy/mok)^2)	;E(var)+var(E), with E(var)~stddev and var(E)~N*<y>^2
        endif
        if nnwt eq 0 and ndwt gt 0 then begin
          yy=1./yd & if isim eq 0 then ye=yy*sqrt(yde^2/yd^4+mok/(yy/mok)^2)
        endif
        if nnwt gt 0 and ndwt gt 0 then begin
          yy=yn/yd & if isim eq 0 then ye=yy*sqrt((yne/yd)^2+(yn*yde/yd^2)^2+mok*(yy/mok)^2)
        endif
        ;
        yy=yy/hw/2.
        if isim eq 0 then ye=ye/hw/2.
	;
	if isim eq 0 and keyword_set(crossval) then begin	;(compute array for cross-validation
	  ok0=where((xx-xx[i]) ge -hw and (xx-xx[i]) le 0,mok0)
	  ok1=where((xx-xx[i]) ge 2*hw and (xx-xx[i]) le 3*hw,mok1)
	  j00=ok0[0] & j01=ok0[mok0-1L] & j10=ok1[0] & j1=ok1[mok1-1L]
	  if mok0 gt 0 and mok1 gt 0 then ok2=[ok0,ok1] else $
	   if mok0 gt 0 and mok1 eq 0 then ok2=ok0 else $
	    if mok0 eq 0 and mok1 gt 0 then ok2=ok1 else $
	     ok2=[-1L] & mok2=mok0+mok1
	  if mok2 gt 0 then begin
	    yy2=float(mok2)
	    ye2=sqrt(mok2) & yne2=0 & yde2=0
	    ;
	    if nnwt gt 0 then begin & yn2=total(zn[ok2]) & if mok2 gt 1 then yne2=stddev(zn[ok2]) & endif
	    if ndwt gt 0 then begin & yd2=total(zd[ok2]) & if mok2 gt 1 then yde2=stddev(zd[ok2]) & endif
	    if nnwt gt 0 and ndwt eq 0 then begin
	      yy2=yn2 & ye2=sqrt(yne2^2*mok2*(yy2/mok2)^2)
	    endif
            if nnwt eq 0 and ndwt gt 0 then begin
              yy2=1./yd2 & ye2=sqrt(yde2^2/yd2^4+mok2/(yy2/mok2)^2)
            endif
            if nnwt gt 0 and ndwt gt 0 then begin
              yy2=yn2/yd2 & ye2=sqrt((yne2/yd2)^2+(yn2*yde2/yd2^2)^2+mok2*(yy2/mok2)^2)
            endif
	    ;
	    yy2=yy2/hw/2.
	    ye2=ye2/hw/2.
	  endif
	endif							;ISIM=0 && CROSSVAL1)
      endif
    endif else begin			;USEX)(do not USE X
      i0=(i-hw)>0
      i1=(i+hw)<(nx-1L)
      mok=i1-i0+1L
      dx=xx[i1]-xx[i0]
      yy=float(mok) & if isim eq 0 then ye=sqrt(mok) & yne=0 & yde=0
      ;
      if nnwt gt 0 then begin & yn=total(zn[i0:i1]) & if mok gt 1 and isim eq 0 then yne=stddev(zn[i0:i1]) & endif
      if ndwt gt 0 then begin & yd=total(zd[i0:i1]) & if mok gt 1 and isim eq 0 then yde=stddev(zd[i0:i1]) & endif
      if nnwt gt 0 and ndwt eq 0 then begin
        yy=yn & if isim eq 0 then ye=sqrt(yne^2+mok*(yy/mok)^2)
      endif
      if nnwt eq 0 and ndwt gt 0 then begin
        yy=1./yd & if isim eq 0 then ye=yy*sqrt(yde^2/yd^4+mok/(yy/mok)^2)
      endif
      if nnwt gt 0 and ndwt gt 0 then begin
        yy=yn/yd & if isim eq 0 then ye=yy*sqrt((yne/yd)^2+(yn*yde/yd^2)^2+mok*(yy/mok)^2)
      endif
      ;
      if dx gt 0 then yy=yy/dx
      if dx gt 0 and isim eq 0 then ye=ye/dx
      ;
      if isim eq 0 and keyword_set(crossval) then begin	;(compute array for cross-validation
        j00=i0 & j01=i & j10=(i+2*hw)<(nx-1) & j11=(i+3*hw)<(nx-1)
	mok0=j01-j00+1L & mok1=j11-j10+1L & mok2=mok0+mok1
	dx0=xx[j01]-xx[j00] & dx1=xx[j11]-xx[j10] & dx2=dx0+dx1
	yy2=float(mok2) & ye2=sqrt(mok2) & yne2=0 & yde2=0
	;
	if mok0 gt 0 and mok1 gt 0 then zn2=[zn[j00:j01],zn[j10:j11]] else $
	 if mok0 gt 0 and mok1 eq 0 then zn2=zn[j00:j01] else $
	  if mok0 eq 0 and mok1 gt 0 then zn2=zn[j10:j11]
	if mok0 gt 0 and mok1 gt 0 then zd2=[zd[j00:j01],zd[j10:j11]] else $
	 if mok0 gt 0 and mok1 eq 0 then zd2=zd[j00:j01] else $
	  if mok0 eq 0 and mok1 gt 0 then zd2=zd[j10:j11]
        if nnwt gt 0 then begin & yn2=total(zn2) & if mok2 gt 1 then yne2=stddev(zn2) & endif
        if ndwt gt 0 then begin & yd2=total(zd2) & if mok2 gt 1 then yde2=stddev(zd2) & endif
        if nnwt gt 0 and ndwt eq 0 then begin
          yy2=yn2 & ye2=sqrt(yne2^2+mok2*(yy2/mok2)^2)
        endif
        if nnwt eq 0 and ndwt gt 0 then begin
          yy2=1./yd2 & ye2=yy2*sqrt(yde2^2/yd2^4+mok2/(yy2/mok2)^2)
        endif
        if nnwt gt 0 and ndwt gt 0 then begin
          yy2=yn2/yd2 & ye2=yy2*sqrt((yne2/yn2)^2+(yn2*yde2/yd2^2)^2+mok2*(yy2/mok2)^2)
        endif
        ;
        if dx gt 0 then yy2=yy2/dx2
        if dx gt 0 then ye2=ye2/dx2
      endif						;ISIM=0 && CROSSVAL)
    endelse				;I0,I1)
    lbound[i]=i0 & rbound[i]=i1

    if isim eq 0 then begin
      yrun[i]=yy
      yerr[i]=ye
      if keyword_set(crossval) then yrun2[i]=yy2
      if keyword_set(crossval) then yerr2[i]=ye2
    endif else ysim[isim-1L,i]=yy

  endfor				;I=0,NX-1}
endfor				;ISIM=0,MSIM}
if msim gt 0 then for i=0L,nx-1L do yesim[i]=stddev(ysim[*,i])

;	now compute the cross-validation mean-square error
;	for the running mean, there are K=NX/(2*HW) independent points.
;	so pick these K points and compute mean of squared difference
;	between estimate vs gapped estimate.  The gapped estimates are
;	simply the ones where the Kth bin is removed and the running
;	mean is computed by conflating the next bin with the previous one.
if keyword_set(crossval) then begin
  if keyword_set(usex) then kk=long((max(xx)-min(xx))/hw/2.+1)>1 else kk=(1+nx/hw/2)>1
  jkk=lindgen(kk+1)*float(nx)/float(kk+1)
  ;ikk=[(jkk-1L)>0,jkk,(jkk+1L)<(nx-1L)] & ikk=ikk[uniq(ikk,sort(ikk))]
  ikk=jkk
  yk=yrun[ikk] & yk2=yrun2[ikk] & yek=yerr[ikk]
  ;cvmse=mean(((yk-yk2)/yek)^2,/nan)	;don't do this, it'll "correct" for a crucial factor!
  cvmse=mean(((yk-yk2))^2,/nan)
endif

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,yrun
end
