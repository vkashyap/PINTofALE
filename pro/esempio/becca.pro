;+
;BECCA.PRO
;
;example program for calling LINESPEC and CONT_CIE
;							vinay kashyap
;-

if not keyword_set(fstr) then fstr=1	;to avoid calling RD_LINE many times
if not keyword_set(w0) then w0=10.	;wavelength range [AA]
if not keyword_set(w1) then w1=800.	;wavelength range [AA]
if not keyword_set(dbdir) then dbdir='$SPEX'	;line DB directory
if n_elements(abund) ne 30 then abund=getabund('anders & grevesse')
if not keyword_set(n_e) then n_e=1e9			;e-density
if not keyword_set(nbin) then nbin=1001			;number of bins
if not keyword_set(tlog) then tlog=7.0			;temperature
if not keyword_set(EM) then EM=1e20			;emission measure
;atom='Fe'	;stick to Iron
print,'remember to set FSTR=0 if you change any of'
help,w0,w1,dbdir

;	define DEM matching this EM
nT=81L & dlogT=0.05 & logT=findgen(nT)*dlogT+4. & tmp=min(abs(tlog-logT),it)
DEM=dblarr(nT) & DEM(it)=EM/alog(10.)/dlogT

;	get line spectrum [ergs/cm^2/s/A]
lsp=linespec(atom,wmin=w0,wmax=w1,nbin=nbin,ww=wwl,dbdir=dbdir,n_e=n_e,$
	fstr=fstr,DEM=DEM,/NOPH)
;	convert to [1e-23 ergs/cm^2/s/A]
lsp=1e23*lsp & dwl=wwl(1:*)-wwl
nrgw=6.626e-27*3e10*1e8/wwl

;	alternate line flux
tmp=min(abs(tlog-fstr.logT),iT)
lem=fstr.LINE_INT & fxl=0. & nrgl=6.626e-27*3e10*1e8/abs(fstr.wvl)
for iZ=1,30 do begin
  oZ=where(fstr.Z eq iZ,moZ)
  if moZ gt 0 then fxl=fxl+abund(iZ-1)*total(EM*lem(iT,oZ)/nrgl(oZ)/1e23)
endfor

;	get continuum spectrum [1e-23 ergs cm^3/s/A]
csp=cont_cie(atom,pres,tlog,wwc,wmn=w0,wmx=w1,nbin=nbin,n_e=n_e,abund=abund)
;	convert to [1e-23 ergs/cm^2/s/A]
csp=csp*EM & dwc=wwc(1:*)-wwc
nrgc=6.626e-27*3e10*1e8/(0.5*(wwc(1:*)+wwc))

;	alternate continuum flux
fxc=total(csp*dwc/nrgc/1e23)

;	total flux in line v/s continuum
print,'' & print,'--------------------------------------'
print,'Line v/s Cont',total(lsp*dwl),total(csp*dwc),' [1e-23 ergs/s/cm^2]'
print,total(lsp*dwl/nrgw/1e23),total(csp*dwc/nrgc/1e23),' [ph/s/cm^2]'
print,total(lsp*dwl/nrgw/1e23)/$
	(total(lsp*dwl/nrgw/1e23)+total(csp*dwc/nrgc/1e23)),' line/total {ph}'
;print,fxl,fxc,fxl/fxc

;	plot
loadct,3
wvl=0.5*(wwl(1:*)+wwl)
plot,wvl,csp,psym=10,xtitle='[A]',ytitle='[1e-23 ergs/s/cm^2/keV]',/xl,/yl
oplot,wvl,lsp,psym=10,col=150

print,'remember to set FSTR=0 if you change any of'
help,w0,w1,dbdir

end
