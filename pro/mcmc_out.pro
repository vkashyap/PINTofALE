;+
;MCMC_OUT
;	program to read in the output of MCMC_DEM for further
;	manipulation.
;
;vinay kashyap (Mar97)
;	numerous bug fixes (VK/LL; Dec'02)
;-

;	restore variables
if not keyword_set(savfil) then savfil='mcmc.sav'
c1='' ;& read,prompt='restore from file ['+savfil+']: ',c1
if c1 ne '' and c1 ne ' ' then savfil=c1
restore,savfil
if n_tags(e) gt 0 then begin
  enam=tag_names(e) & os=where(enam eq 'STEPS',mos)
  if mos gt 0 then steps=e.steps
endif

;	reset parameters
	ebound=0.68		;fraction of p(M|D) to include in error bar
;	hwhm=1			;compute upper+lower half-widths@half-max
	smooth=0		;smooth histograms or not?
if not keyword_set(parthresh) then parthresh=0.25

;	outputs
mapar=double(allpar)			;MAP estimates of parameters
;tmp=varsmooth(allpar(0:nt-1),lscal,_extra=e)
if not keyword_set(loopy) then $
  tmp=mk_dem('varsmooth',logT=logT,indem=allpar(0:nt-1),pardem=lscal,_extra=e) else $
  tmp=mk_dem('loop',indem=allpar(0:nt-1),pardem=logT,_extra=e)
mapar(0:nt-1)=tmp
upar=mapar				;upper confidence level
lpar=upar				;lower confidence level
fupar=upar				;upper HalfWidth@HalfMax
flpar=lpar				;lower HalfWidth@HalfMax
for i=0,npar-1 do begin			;{shtep thru parameters
  j=opar(i)
  if j lt nt then x=xdem else x=xab
  if j lt nt then xx=0.5*(xdem(1:*)+xdem) else xx=0.5*(xab(1:*)+xab)
  ;
  oh=where(storidx eq j,moh)
  if moh gt long(parthresh*nsim) then begin	;(any reason to bother?
    whee,spoke=spoke,/moveit
    g=storpar(oh) & nx=n_elements(xx) & gf=fltarr(nx)
    for k=0,nx-1 do begin		;probability distribution for parameter
      oo=where(g ge x(k) and g lt x(k+1),moo) & gf(k)=moo
    endfor
    f=gf
    if keyword_set(smooth) then begin
      ;whee,spoke=spoke,/moveit
      ;smth=findscale(gf) & o0=where(gf lt 1e-4*max(gf),mo0)
      ;;if mo0 gt 0 then f(o0)=varsmooth(gf(o0),smth(o0),steps=steps) else f=varsmooth(gf,smth,steps=steps)
      ;if not keyword_set(loopy) then begin
      ;  if mo0 gt 0 then $
      ;	  f(o0)=mk_dem('varsmooth',logT=logT,indem=gf(o0),pardem=smth(o0),steps=steps,_extra=e) else $
      ;	  f=mk_dem('varsmooth',logT=logT,indem=gf,pardem=smth,steps=steps,_extra=e)
      ;endif else begin
      ;  if mo0 gt 0 then $
      ;   f(o0)=mk_dem('loop',indem=gf(o0),pardem=logT,_extra=e) else $
      ;	  f=mk_dem('loop',indem=gf,pardem=logT,_extra=e)
      ;endelse
    endif
    f=float(f)/total(f)
    cf=f & for k=1,n_elements(f)-1 do cf(k)=cf(k-1)+f(k)
    cf=cf/max(cf)				;cumulative distribution
    ;
    ;tmp=max(f,ip) & mapar(j)=xx(ip)	;MAP estimate
    ;mapar(j)=total(g)/moh & tmp=min(abs(xx-mapar(j)),ip)
    if j lt nt then mapar(j)=alog10(bestdem(j)) & tmp=min(abs(xx-mapar(j)),ip)
    ;
    ;	error bounds
    ebu=(cf(ip)+ebound/2.)<1 & ebl=(cf(ip)-ebound/2.)>0
    ipu=min(where(cf ge ebu)) & ipl=max(where(cf le ebl))>0
    if ipu eq -1 then ipu=nx-1
    upar(j)=xx(ipu) & lpar(j)=xx(ipl)
    ;
    ;	HWHMs
    oo=where(f ge f(ip)/2.)
    ifl=min(oo,max=ifu)
    fupar(j)=xx(ifu) & flpar(j)=xx(ifl)
  endif					;MOO>0.25*NSIM)
endfor					;I=0,NPAR-1}

;	recast outputs
dem=mapar(0:nt-1) & abund=mapar(nt:*)
demerr=dblarr(nt,2) & demerr(*,0)=lpar(0:nt-1) & demerr(*,1)=upar(0:nt-1)
aberr=fltarr(nab,2) & aberr(*,0)=lpar(nt:*) & aberr(*,1)=upar(nt:*)
if keyword_set(hwhm) then begin
  demerr(*,0)=flpar(0:nt-1) & demerr(*,1)=fupar(0:nt-1)
  aberr(*,0)=flpar(nt:*) & aberr(*,1)=fupar(nt:*)
endif
dem=10.D^(dem) & abund=10.^(abund) & demerr=10.D^(demerr) & aberr=10.^(aberr)

;	plot
if !d.name eq 'X' then window,1
tmp=pred_flx(line,logT,ww,dem,Z=zz,fobs=fx,fsigma=sig,abund=abund)
;
if !d.name eq 'X' then window,0
oo=where(demerr(*,0) ne demerr(*,1),moo)
plot,logt(oo),dem(oo),yr=[min(demerr),max(demerr)],/yl,/ys,psym=10,$
	xr=[min(logt),max(logt)],title=savfil,$
	xtitle='log!d10!n(T [K])',ytitle='DEM [cm!u-5!n]';,col=150
	;xr=[min(logt(oo)),max(logt(oo))],title=savfil,$
for i=0,moo-1 do oplot,logt(oo(i))*[1,1],demerr(oo(i),*);,col=150

end
