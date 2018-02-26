;+
;MAKE_SPEC.PRO
;	take the line and continuum emissivities and generate a
;	theoretical spectrum over the specified wavelength grid
;
;restrictions
;	output wavelength grid must be defined in WGRID
;	line database must be in structure LINSTR
;	continuum database must be in structure CONSTR
;
;subroutines
;	LINEFLX
;	    GETABUND
;	        SETABUND
;	    WHEE
;	ISMTAU
;	HASTOGRAM
;	REBINW
;
;vinay kashyap (Dec98)
;-

;	check for required inputs
if n_tags(linstr) eq 0 then stop,'Need line database structure, LINSTR'
	wlin=abs(linstr.WVL) & nlin=n_elements(wlin)
if n_tags(constr) eq 0 then stop,'Need continuum database structure, CONSTR'
	wcon=abs(constr.midWVL) & ncon=n_elements(wcon)
nwgrid=n_elements(wgrid)
if nwgrid lt 2 then begin
  print,'Using same grid as in continuum database'
  wgrid=constr.WVL
endif
dgrid=abs(wgrid(1:*)-wgrid)

;	make sure WGRID is in ascending order (because REBINW wants it so)
oo=where(dgrid lt 0,moo)
if moo gt 0 then begin
  print,'' & print,'Reversing WGRID!'
  wgrid=reverse(wgrid)
  dgrid=abs(wgrid(1:*)-wgrid)
endif

;	get predicted line fluxes	[ph/s/cm^2]
if n_elements(flx) ne nlin then $
  flx=lineflx(linstr.LINE_INT,linstr.LOGT,linstr.WVL,linstr.Z,$
  DEM=DEM,abund=abund, effar=areff,wvlar=arwvl) else $
  print,'bypassing new calculation of FLX'

;	get predicted continuum fluxes	[ph/s/cm^2/Ang]
if n_elements(fcx) ne ncon then $
  fcx=lineflx(constr.CONT_INT,constr.LOGT,constr.midWVL,$
  DEM=DEM,abund=abund, effar=areff,wvlar=arwvl) else $
  print,'bypassing new calculation of FCX'

;	get column density
etaul=fltarr(nlin)+1. & etauc=fltarr(ncon)+1.
if keyword_set(NH) then begin
  if NH lt 30 then NH=10.^(NH)
  if not keyword_set(fH2) then fH2=0.26
  if not keyword_set(He1) then He1=0.1*NH
  if not keyword_set(HeII) then HeII=0.1*He1
  etaul=exp(-ismtau(wlin,NH=NH,fH2=fH2,He1=He1,HeII=HeII,Fano=Fano))
  etauc=exp(-ismtau(wcon,NH=NH,fH2=fH2,He1=He1,HeII=HeII,Fano=Fano))
endif
fl=flx*etaul & fc=fcx*etauc

;	make line spectrum
spl=hastogram(wlin,wgrid,wts=fl)	;[ph/s/cm^2/bin]

;	make continuum spectrum by rebinning
spc=rebinw(fc,constr.WVL,wgrid)*dgrid	;[ph/s/cm^2/Ang]*[Ang/bin]

;	final spectrum	[ph/s/cm^2/bin]
spec=spl+spc

;	say what's which
help,linstr,wlin,etaul,flx,fl,/str
help,constr,wcon,etauc,fcx,fc,/str
help,wgrid,spl,spc,spec

;	plot
if not keyword_set(xtitle) then xtitle='[Ang]'
if not keyword_set(ytitle) then ytitle='[ph s!u-1!n cm!u-2!n bin!u-1!n]'
if not keyword_set(title) then title='Simulated Spectrum'
if not keyword_set(xlog) then xlog=0
if not keyword_set(ylog) then ylog=0
plot,wgrid,spec,/xs,/ys,xlog=xlog,ylog=ylog,$
	xtitle=xtitle,ytitle=ytitle,title=title

end
