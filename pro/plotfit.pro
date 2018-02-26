pro plotfit,x,y,fitstr,ymod,ysig=ysig,conlev=conlev,consig=consig,$
	funcs=funcs,intens=intens,scol=scol,ecol=ecol,ccol=ccol,mcol=mcol,$
	ncol=ncol,imcol=imcol,clin=clin,mlin=mlin,nlin=nlin,$
	missing=missing,fwhm=fwhm,betap=betap,angle=angle,$
	submodel=submodel, _extra=e
;+
;procedure	plotfit
;	takes the output of FITLINES and makes a nice plot of it
;	with the data and the model(s) and everything.
;
;syntax
;	plotfit,x,y,fitstr,ymod,ysig=ysig,conlev=conlev,consig=consig,$
;	funcs=funcs,intens=intens,scol=scol,ecol=ecol,ccol=ccol,mcol=mcol,$
;	ncol=ncol,imcol=imcol,clin=clin,mlin=mlin,nlin=nlin,$
;	missing=missing,fwhm=fwhm,betap=betap,angle=angle,submodel=submodel,$
;	/nuthin,/noday,/notime,/nouser,/nopack,stacol=stacol,stasiz=stasiz,$
;	stathk=stathk,$
;	PLOT_KEYWORDS
;
;parameters
;	x	[INPUT; required] absissa, e.g., wavelength or energy
;	y	[INPUT; required] count spectrum Y(X)
;		* size of Y must match that of X
;	fitstr	[INPUT] output of FITLINES
;	ymod	[OUTPUT] model function
;
;keywords
;	ysig	[INPUT] errors on Y
;		* default is 0.
;		* plotted if supplied, otherwise not.
;	funcs	[INPUT] name of user defined function that takes as input
;		X and A (the model parameters), and returns Y(X;A) and
;		the partial derivatives dY/dA
;		* default is set to X3MODEL
;	intens	[INPUT] if set, Y(X) is assumed to be intensity, i.e.,
;		(counts|ergs|.../bin_unit/...)
;		* default is to take Y(X) to be (counts|ergs|.../bin/...)
;		* passed straight to FITLINES_EVENT
;		* WARNING: this has >not< been road-tested!
;	conlev	[INPUT] the continuum that was taken out of the spectrum
;		* CONLEV must match the size of X and Y else ignored
;		* if ignored, interpolated from FITSTR.CONLEV
;	consig	[INPUT] error on CONLEV
;		* default for CONSIG is sqrt(abs(CONLEV)+0.75)+1.
;		* currently not implemented.  i.e., it's not plotted even
;		  if given.
;	scol	[INPUT; default=6] color of spectrum plot
;	ecol	[INPUT; default=SCOL] color of error-bars
;	ccol	[INPUT; default=1] color of continuum
;	mcol	[INPUT; default=2] color of best-fit model+continuum
;	ncol	[INPUT; default=3] color of best-fit model alone
;	imcol	[INPUT; default=4] color to plot model components
;		* individual model components are not plotted unless
;		  IMCOL is defined and non-zero
;	clin	[INPUT; default=1] line style of continuum
;	mlin	[INPUT; default=0] line style of model+continuum
;	nlin	[INPUT; default=2] line style of model alone
;
;	missing	[INPUT] passed w/o comment to FUNCS
;	fwhm	[INPUT] passed w/o comment to FUNCS
;	betap	[INPUT] passed w/o comment to FUNCS
;	angle	[INPUT] passed w/o comment to FUNCS
;	submodel[OUTPUT] the model function for each component,
;		returned as a (NCOMP,NX) sized array
;		* calculated iff IMCOL is set
;
;	_extra	[INPUT ONLY] use this to pass defined keywords to
;		PLOT and STAMPLE
;
;restrictions
;	requires subroutines
;		FUNCS (X3MODEL [MK_3MODEL, [MK_GAUSS, MK_LORENTZ]])
;	PLOT cannot handle any keyword that is either strange or undefined
;
;known bugs
;	cannot handle spectra with varying bin sizes
;	INTENS has not been tested
;
;history
;	vinay kashyap (MarMM)
;	added call to STAMPLE (VK; FebMMI)
;	bug correction; changed default colors; added keyword IMCOL
;	  (VK; May'03)
;	added keyword SUBMODEL (VK; NovMMVI)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & ny=n_elements(y)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X undefined' else $
  if ny eq 0 then ok='Y undefined' else $
   if nx ne ny then ok='Y[X] and X are incompatible' else $
    if nx lt 2 then ok='Y[X] might as well be undefined'
if ok ne 'ok' then begin
  print,'Usage: plotfit,x,y,fitstr,ymod,ysig=ysig,conlev=conlev,consig=consig,$'
  print,'       funcs=funcs,intens=intens,scol=scol,ecol=ecol,ccol=ccol,mcol=mcol,$'
  print,'       ncol=ncol,imcol=imcol,clin=clin,mlin=mlin,nlin=nlin,$'
  print,'       missing=missing,fwhm=fwhm,betap=betap,angle=angle,$'
  print,'       /nuthin,/noday,/notime,/nouser,/nopack,$'
  print,'       stacol=stacol,stasiz=stasiz,stathk=stathk, PLOT KEYWORDS'
  if np ne 0 then message,ok,/info
  return
endif

;	check input
mk='ok' & mf=n_elements(fitstr) & nf=n_tags(fitstr)
if mf eq 0 then mk='output from FITLINES undefined' else $
 if nf eq 0 then mk='Fit parameters not in a structure' else begin
   namef=tag_names(fitstr)
   ipos=(where(namef eq 'POS'))[0]
   iwdt=(where(namef eq 'WDT'))[0]
   iflx=(where(namef eq 'FLX'))[0]
   ityp=(where(namef eq 'TYPE'))[0]
   icon=(where(namef eq 'CONLEV'))[0]
   iconx=(where(namef eq 'CONLEVX'))[0]
   if ipos lt 0 then mk='Missing positions in fit structure'
   if iwdt lt 0 then mk='Missing widths in fit structure'
   if iflx lt 0 then mk='Missing fluxes in fit structure'
 endelse

;	keywords
; error-bars on data
sigy=fltarr(ny) & nys=n_elements(ysig) & if nys eq ny then sigy=ysig
; /bin or /unit
dWVL=abs(median(x(1:*)-x)) & if keyword_set(intens) then dWVL=1.
if keyword_set(intens) then norm=0 else norm=1
; continuum
clev=fltarr(ny)				;default is 0
if mk eq 'ok' then begin		;fit structure exists..
  if iconx ge 0 and icon ge 0 then begin	;in correct format
    fclev=fitstr.(icon) & fclevx=fitstr.(iconx)
    clev=interpol(fclev,fclevx,x)		;interpolate
  endif
endif
;but use CONLEV if given
ncon=n_elements(conlev) & if ncon eq nx then clev=conlev
if ncon eq 1 then clev[*]=conlev[0]

;	set up the defaults for the plot
xmin=min(x,max=xmax) & xr=[xmin,xmax]
ymin=min(y,max=ymax) & yr=[ymin,ymax]
if not keyword_set(scol) then scol=6
if not keyword_set(ecol) then ecol=scol
if not keyword_set(ccol) then ccol=1
if not keyword_set(mcol) then mcol=2
if not keyword_set(ncol) then ncol=3
if not keyword_set(clin) then clin=1
if not keyword_set(mlin) then mlin=0
if not keyword_set(nlin) then nlin=2
plot,xr,yr,/nodata,xr=xr,yr=yr, _extra=e

;	stamp the plot
stample, _extra=e

;	overplot the spectrum
oplot,x,y,psym=10,color=long(scol[0])

;	overplot the error-bars
for i=0L,ny-1L do oplot,x[i]*[1,1],y[i]+sigy[i]*[1,-1],color=long(ecol[0])

;	overplot the continuum
oplot,x,clev,color=long(ccol[0]),linestyle=fix(clin[0])

;	now make the model curves
if mk eq 'ok' then begin
  if not keyword_set(funcs) then funcs='x3model'	;model function name
  pos=fitstr.(ipos) & wdt=fitstr.(iwdt) & flx=fitstr.(iflx)
  ncomp=n_elements(pos)
  if ityp lt 0 then type=strarr(ncomp)+'gauss' else type=fitstr.(ityp)
  a=fltarr(3L*ncomp) & ii=lindgen(ncomp)
  a[3L*ii+0L]=pos[*] & a[3L*ii+1L]=wdt[*] & a[3L*ii+2L]=flx[*]*dWVL
  call_procedure,funcs,x,a,ymod,norm=norm,type=type,$
	missing=missing,fwhm=fwhm,betap=betap,angle=angle
  oplot,x,ymod+clev,color=long(mcol[0]),linestyle=fix(mlin[0])
  oplot,x,ymod,color=long(ncol[0]),linestyle=fix(nlin[0])
  if keyword_set(imcol) then begin
    submodel=fltarr(ncomp,nx)
    for i=0L,ncomp-1L do begin
      a=fltarr(3L) & a[0]=pos[i] & a[1]=wdt[i] & a[2]=flx[i]*dWVL
      call_procedure,funcs,x,a,yimod,norm=norm,type=type[i],$
	missing=missing,fwhm=fwhm,betap=betap,angle=angle
      oplot,x,yimod,color=long((imcol[[i]])[0]),linestyle=fix(mlin[0])
      submodel[i,*]=yimod
    endfor
  endif
endif

return
end
