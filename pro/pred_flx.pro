function pred_flx,line,logT,wvl,DEM,Z=Z,NH=NH,fobs=fobs,fsigma=fsigma,$
	ulim=ulim,verbose=verbose, _extra=e
;+
;function	pred_flx
;	return line fluxes for a given DEM and compare with observed fluxes.
;
;syntax
;	f=pred_flx(line,logT,wvl,DEM,Z=Z,NH=NH,fobs=fobs,fsigma=fsigma,$
;	ulim=ulim,verbose=verbose,effar=effar,wvlar=wvlar,fH2=fH2,He1=He1,$
;	HeII=HeII,/fano)
;
;parameters
;	line	[INPUT; required] line cooling emissivities (T,Wvl) in units
;		of 1e-23 erg cm^3/s, including ion balance, but not abundances.
;	logT	[INPUT; required] array of log10(Temperature [K]) at which
;		emissivities are given.
;	wvl	[INPUT; required] wavelength of matching spectral line(s) [Ang]
;	DEM	[INPUT] Differential Emission Measure at each T [cm^-5/logK]
;		* passed w/o change to LINEFLX
;
;keywords
;	Z	[INPUT] atomic number of element that generates each WVL
;		* passed w/o change to LINEFLX
;	NH	[INPUT] H column density [cm^-2]
;	fobs	[INPUT] observed fluxes; if supplied, displays comparison
;		plots of predicted v/s observed
;	fsigma	[INPUT] errors on observed fluxes; if not supplied, assumed
;		to be 1+sqrt(abs(FOBS)+0.75)
;	ulim	[INPUT] flag to indicate whether a given FOBS is an upper
;		limit (1: UL, 0: not)
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT] allows setting defined keywords to
;		LINEFLX [EFFAR, WVLAR, TEMP, ABUND, KEV, NOPH]
;		ISMTAU [FH2,HE1,HEII,FANO]
;
;restrictions
;	* requires subroutines
;	  GETABUND
;	  LINEFLX [WHEE]
;	  ISMTAU
;
;history
;	vinay kashyap (May97)
;	changed keyword name SIGMA to FSIGMA (VK; MMaug)
;	changed Fobs/Fpred from /ylog to plotting alog10() (VK; MMVsep)
;	changed Observed v/s predicted to plot psym=4
;	bug correction: matching up upper limit points (VK; MMVIapr)
;	added keyword VERBOSE (VK; MMVIIjun)
;-

;	usage
nl=n_elements(line) & nw=n_elements(wvl) & nt=n_elements(logT)
nd=n_elements(dem)
if nl eq 0 or nw eq 0 or nt eq 0 or nd eq 0 then begin
  print,'Usage: f=pred_flx(line,logT,wvl,DEM,Z=Z,NH=NH,fobs=fobs,fsigma=fsigma,$
  print,'       ulim=ulim,verbose=verbose)'
  print,'  returns predicted fluxes for given DEM'
  print,'also accepts defined keywords:'
  print,'  EFFAR,WVLAR (LINEFLX) and FH2,HE1,HEII,FANO (ISMTAU)'
  return,-1L
endif

;	error checks
szl=size(line) & nlt=szl(1) & nlw=szl(2)
if nlt ne nt then begin
  message,'temperature grid mismatch',/info & return,0*wvl
endif
if nlw ne nw then begin
  message,'wavelength mismatch',/info & return,0*wvl
endif

;	keywords, etc.
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
nz=n_elements(z) & zz=intarr(nw)+1
if nz lt nw then zz(0:nz-1)=z(*) else zz=z(0:nw-1)
ww=[abs(wvl)]

;	now to get counts...
fx=dblarr(nw)
if nt eq 1 then begin
  for i=0,nw-1 do fx(i)=lineflx(line(i),logT,ww(i),zz(i),DEM=DEM,_extra=e)
endif else fx=lineflx(line,logT,ww,zz,DEM=DEM, _extra=e)

;	get optical depths
tau=0*ww
if keyword_set(NH) then tau=ismtau(ww,NH=NH,fH2=0,/Fano,_extra=e) < 69.
fx=fx*exp(-tau)

;	compare fluxes
if n_elements(fobs) eq nw then begin		;{plot

  if n_elements(fsigma) ne nw then fsigma=1.+sqrt(abs(fobs)+0.75)
  if n_elements(ulim) ne nw then ulim=intarr(nw)

  pmult=!p.multi & !p.multi=[0,2,3]
  oo=sort(ww) & w=ww(oo) & flx=fx(oo) & ff=fobs(oo) & sig=fsigma(oo)
  emis=line(*,oo) & ul=ulim(oo)

  ou=where(ul gt 0,mou)

  ;plot #1 -- "spectrum"
  dw=abs(min(w(1:*)-w)) < (1./float(nw))
  if vv gt 0 then begin
    plot,w,ff,/yl,psym=7,title='"spectrum"',charsize=2
    for i=0,nw-1 do oplot,w(i)+dw*[-1,0,1],flx(i)*[0,1,0],psym=10
    for i=0,nw-1 do oplot,w(i)*[1,1],ff(i)+sig(i)*[-1,1]
    if mou gt 0 then oplot,w(ou),ff(ou),psym=4
  endif

  ;plot #2 -- "delta_chi (observed-predicted)/err(obs)"
  chi=(ff-flx)/sig & chi2=total(chi^2)
  if vv gt 0 then plot,w,chi,psym=10,title='(obs-pred)/err(obs) ('+strtrim(chi2,2)+')',$
	charsize=2

  ;plot #3 -- ratios
  rat=ff/flx & ratu=(ff+sig)/flx & ratl=(ff-sig)/flx
  ;plot,w,rat,psym=4,/yl,title='F!dobs!n/F!dpred!n',charsize=2,/ys
  ;for i=0,nw-1 do oplot,w(i)*[1,1],[ratu(i),ratl(i)]
  ;if mou gt 0 then oplot,w(ou),rat(ou),psym=6
  ;if mou gt 0 then oplot,w(ou),rat(ou),psym=7,col=100
  if vv gt 0 then begin
    plot,w,alog10(rat),psym=4,title='log!d10!n(F!dobs!n/F!dpred!n)',charsize=2,/ys
    for i=0,nw-1 do oplot,w(i)*[1,1],alog10([ratu(i),ratl(i)])
    if mou gt 0 then oplot,w(ou),alog10(rat(ou)),psym=6
    if mou gt 0 then oplot,w(ou),alog10(rat(ou)),psym=7,col=100
  endif

  ;plot #4 -- obs v/s thr
  if vv gt 0 then begin
    plot,ff,flx,psym=4,/xl,/yl,xtitle='F!dobs!n',ytitle='F!dpred!n',$
	title='observed v/s predicted',charsize=2
    oo=sort(ff)
    oplot,ff(oo),ff(oo)
    oplot,ff(oo),2.*ff(oo),col=100 & oplot,ff(oo),ff(oo)/2.,col=100
    if mou gt 0 then oplot,ff(ou),flx(ou),psym=1
  endif

  ;plot #5 -- observed & predicted v/s peak temperatures
  tmax=fltarr(nw)
  for i=0,nw-1 do begin
    lmx=max(emis(*,i),imx) & tmax(i)=logt(imx)
  endfor
  ot=sort(tmax)
  if vv gt 0 then begin
    plot,ff(ot),/yl,psym=7,title='F!dobs!n & F!dpred!n arranged for T!dmax!n',$
	charsize=2
    for i=0,nw-1 do oplot,[i,i],ff(ot(i))+sig(ot(i))*[-1,1],col=150
    for i=0,nw-1 do oplot,[i,i],flx(ot(i))*[1e-10,1]
  endif
  if mou gt 0 then begin
    tmaxu=tmax(ot) & chiu=chi(ot) & ffu=ff(ot) & ulu=ul(ot)
    ouu=where(ulu gt 0)
    if vv gt 0 then for i=0,mou-1 do oplot,[ouu(i)],[ffu(ouu(i))],psym=4
    ;for i=0,mou-1 do oplot,[ou(i)],[ff(ot(ou(i)))],psym=4
    ;tmaxu=tmax(ou) & chiu=chi(ou) & ffu=ff(ou) & out=sort(tmaxu)
  endif

  ;plot #6 -- "delta_chi" v/s peak temperatures
  if vv gt 0 then begin
    plot,tmax(ot),chi(ot),psym=7,title='deltaF/error v/s T!dmax!n',charsize=2
    if mou gt 0 then oplot,tmaxu(ouu),chiu(ouu),psym=6,col=100
  endif

  !p.multi=pmult
endif						;done plotting}

return,fx
end
