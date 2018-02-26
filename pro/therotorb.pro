function therotorb,wvl,lstr,x,logT=logT,DEM=DEM,vrot=vrot,vorb=vorb,$
	allem=allem,nx=nx,dx=dx,fwhm=fwhm,verbose=verbose,$
	therm_flux=therm_flux,therm_prof=therm_prof,thnorm=thnorm,$
	nonorm=nonorm, _extra=e
;+
;function	therotorb
;	compute the profile by which a delta-function would be broadened
;	and shifted due to thermal, rotational, and orbital motions
;
;syntax
;	p=therotorb(wvl,lstr,x,logT=logT,DEM=DEM,vrot=vrot,vorb=vorb,$
;	/allem,nx=nx,dx=dx,fwhm=fwhm,verbose=verbose,therm_flux=therm_flux,$
;	therm_prof=therm_prof,/thnorm,/nonorm)
;
;parameters
;	wvl	[INPUT; required] wavelength of line
;		* must be a scalar -- if vector, only first element is used
;	lstr	[INPUT; required] line emissivity structure for given line
;		* if there is more than one line in the structure, then
;		  -- if ALLEM is set, all the emissivities are summed
;		  -- if ALLEM is not set, then the closest line to WVL is
;		     looked for.  this may not be what is needed, e.g., for
;		     higher-order grating lines, in which case the structure
;		     should be pre-filtered in using CAT_LN().
;	x	[OUTPUT] wavelengths at which profile is computed
;		* uses NX and DX to determine grid
;
;keywords
;	logT	[INPUT] log10(T[k]) at which DEM is defined
;	DEM	[INPUT; default=constant] differential emission measure
;		[cm^-5/logK]
;		* size must match LOGT
;	vrot	[INPUT; default=0] rotational velocity in [km/s]
;	vorb	[INPUT; default=0] orbital velocity in [km/s]
;		* -VE values get RED shifted
;	allem	[INPUT] if set, sums up all the emissivities in LSTR.LINE_INT
;		into a single array
;	nx	[INPUT; default=1001L] number of bins in X
;	dx	[INPUT; default=WVL/500] range of X
;	fwhm	[OUTPUT] FWHM of resulting profile
;	verbose	[INPUT] controls chatter
;	therm_flux	[OUTPUT] flux at each temperature, the product
;			of emissivity and DEM
;			* NOTE that abundance is not included unless it
;			  is already folded in via the emissivity
;	therm_prof	[OUTPUT] thermally broadened line profile
;			* normalized to 1 if THNORM is set
;	thnorm	[INPUT] if set, ensures that the integral of THERM_PROF=1
;	nonorm	[INPUT] if set, the output broadened profile is not
;		renormalized in addition to THNORM
;		* the default, for historical reasons, is to set the
;		  sum of binned profile to 1
;	_extra	[JUNK] here only to avoid crashing the program
;
;example usage
;	lstr=rd_list('O7|22.1+-0.01|$CHIANTI',sep='|',n_e=1e10,/incieq)
;	p=therotorb(22.1,lstr,x,verbose=10,vrot=30,vorb=-50)
;	l=mk_lorentz(x,22.1,0.005,10.,betap=2.5)
;	pl=voorsmooth(l,1,x,usekern=p,eps=1e-4)
;	plot,x,pl,thick=2 & oplot,x,l & oplot,x,p*max(pl)/max(p),line=2
;
;subroutines
;	INICON
;	MK_DEM
;
;history
;	vinay kashyap (Oct01)
;	now does not crash if VORB,VROT=0 (VK; Apr06)
;	added keywords THERM_PROF,THNORM,NONORM (VK; Jul07)
;-

;	usage
ok='ok'
np=n_params() & nw=n_elements(wvl) & nl=n_tags(lstr)
if np lt 2 then ok='Insufficient parameters' else $
 if nw eq 0 then ok='WVL not defined' else $
  if nl eq 0 then ok='LSTR not a line emissivity structure'
if ok ne 'ok' then begin
  print,'Usage: p=therotorb(wvl,lstr,x,logT=logT,DEM=DEM,vrot=vrot,vorb=vorb,$'
  print,'       /allem,nx=nx,dx=dx,fwhm=fwhm,therm_flux=therm_flux,verbose=verbose,$'
  print,'       therm_prof=therm_prof,/thnorm,/nonorm)'
  print,'  compute thermal+rotational+orbital effects on line'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	check inputs
wave=abs(wvl[0])
tagnam=tag_names(lstr)
if tagnam[0] ne 'LINE_INT' or tagnam[1] ne 'LOGT' then begin
  message,'line emissivity structure is not in right format; exiting',/info
  return,-1L
endif
emis=lstr.LINE_INT & ylogT=lstr.LOGT & Z=lstr.Z & ww=abs(lstr.WVL)
mw=n_elements(ww)
if mw gt 1 then begin
  if keyword_set(allem) then begin	;(combine all lines
    emis=0.D*ylogT & Z=0.
    for i=0L,mw-1L do emis=emis+(lstr.LINE_INT)[*,i]
    for i=0L,mw-1L do Z=Z+total((lstr.LINE_INT)[*,i])*(lstr.Z)[i] & Z=long(Z/total(emis))
  endif else begin			;ALLEM=1)(pick nearest line
    dw=min(abs(wave-ww),i) & emis=reform((lstr.LINE_INT)[*,i]) & Z=(lstr.Z)[i]
  endelse				;ALLEM=0)
endif

;	keywords
;	thermal
nD=n_elements(DEM) & nT=n_elements(logT)
if nT eq 0 then begin
  message,'LOGT not defined; assuming default',/info
  xlogT=ylogT & nT=n_elements(xlogT)
endif else xlogT=logT
if nD eq 0 then begin
  message,'DEM not defined; assuming constant',/info
  xDEM=dblarr(nT)+1d10 & nD=n_elements(xDEM)
endif else xDEM=DEM
if nD ne nT then begin
  message,'DEM and LOGT are incompatible; assuming defaults',/info
  xlogT=ylogT & nT=n_elements(xlogT)
  xDEM=dblarr(nT)+1d10 & nD=n_elements(xDEM)
endif
yDEM=mk_dem('interpolate',logT=ylogT,indem=xDEM,pardem=xlogT)
;	rotational
nrot=n_elements(vrot) & rotv=0.
if nrot gt 0 then rotv=vrot[0]
;	orbital
norb=n_elements(vorb) & orbv=0.
if norb gt 0 then orbv=vorb[0]
;	grid
if not keyword_set(nx) then nx=1001L
if not keyword_set(dx) then dx=wave/500.
xmin=wave-dx/2. & xmax=wave+dx/2. & deltax=(xmax-xmin)/nx
x=(findgen(nx)-nx/2)*deltax+wave
;	verbose
vv=0 & if keyword_set(verbose) then vv=abs(verbose[0]) > 1

;	initialize
if n_tags(amu) eq 0 then inicon,atom=atom,amu=amu,fundae=fundae

;	initialize
if n_tags(amu) eq 0 then inicon,atom=atom,amu=amu,fundae=fundae
mu=(amu.(Z[0]-1))[0]

;	compute thermal broadening
therm_prof=fltarr(nx)
therm_flux=fltarr(n_elements(ylogT))
for it=0L,n_elements(ylogT)-1L do begin
  hw_therm=4.301e-7*wave*sqrt(10.^(ylogT[it])/mu)	;Gray 1992, eqn 11.5
  sig_therm=2.*hw_therm/2.355 & tmp=0.*x+100.
  if sig_therm gt 0 then tmp=((x-wave)/sig_therm)^2/2. < 60.
  flx=yDEM[it]*emis[it]
  therm_flux[it]=flx
  if sig_therm gt 0 then therm_prof=therm_prof+flx*exp(-tmp)/sqrt(!pi)/sig_therm
endfor
;therm_prof=therm_prof/total(therm_prof)
;p=therm_prof
if keyword_set(thnorm) then therm_prof=therm_prof/total(therm_prof*deltax)
p=therm_prof
if not keyword_set(nonorm) then p=p/total(p)
if vv gt 5 then plot,x,p

;	compute rotational broadening
if keyword_set(rotv) then begin
  hw_rot=abs((rotv*1e5/fundae.C)*wave)
  tmp=(1.-((x-wave)/hw_rot)^2) > 0
  rot_prof=dx*(2./!pi/hw_rot)*sqrt(tmp)
  ok=where(rot_prof gt 0,mok)
  if mok eq nx then begin
    c1='Rotational profile kernel is same size as wavelength grid;'
    c2='increase DX and NX so that the grid spans a wider range.'
    c3='For now, arbitrarily deleting first element and continuing.'
    message,c1,/informational
    message,c2,/informational
    message,c3,/informational
    ok=ok[1:*]
  endif
  if mok gt 0 then begin
    if keyword_set(nonorm) then begin
      ;	normalize the rotation profile so that it integrates to 1
      p=convol(p,rot_prof[ok]/total(rot_prof[ok]*deltax),/edge_wrap)
    endif else begin
      ;	normalize the rotation profile so that it sums to 1
      p=convol(p,rot_prof[ok]/total(rot_prof[ok]),/edge_wrap)
    endelse
  endif
  if vv gt 5 then oplot,x,p,line=1
endif

;	compute orbital shift
if keyword_set(orbv) then begin
  wshift=(orbv*1e5/fundae.C)*wave
  tmp=min(abs(x-wave),ishft0)
  tmp=min(abs(x-(wave+wshift)),ishft)
  p=shift(p,ishft0-ishft)
  if vv gt 5 then oplot,x,p,line=2
endif

;	find the FWHM from the final profile
pmax=max(p,imx)
if pmax eq 0 then begin
  message,'could not find profile peak; setting FWHM to 0',$
	/informational
  fwhm=0. & return,p
endif
x0=interpol(x[0:imx],p[0:imx],0.5*pmax)
x1=interpol(x[imx:*],p[imx:*],0.5*pmax)
fwhm=(x1-x0)
if not keyword_set(nonorm) then p=p/total(p)

return,p
end
