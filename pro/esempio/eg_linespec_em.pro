;+
;EG_LINESPEC_EM.PRO
;
;example program for calling LINESPEC_EM
;							vinay kashyap
;-

if not keyword_set(fstr) then fstr=1	;to avoid calling RD_LINE many times
w0=50. & w1=100.	;wavelength range in Ang
n_e=1e9		;e-density
dbdir='$SPEX'	;line DB directory
nbin=2501	;number of bins
;atom='Fe'	;stick to Iron
tlog=6.		;temperature
EM=1e20		;emission measure

;	get same spectrum [ph/s/A] in same range again
sp1=linespec_em(atom,wmin=w0,wmax=w1,nbin=nbin,ww=ww,fstr=fstr,tlog=tlog,EM=EM)

;	define DEM matching this EM
logT=fstr.logT & tmp=min(abs(tlog-logT),it)
nT=n_elements(logT) & dlogT=median(logT(1:*)-logT)
DEM=dblarr(nT) & DEM(it)=EM/dlogT/alog(10.)	;dlogT for obv. reasons, and
						;alog(10) to convert to dlnT

;	get same spectrum [ph/s/A] in same range again
sp2=linespec(atom,wmin=w0,wmax=w1,nbin=nbin,ww=ww,fstr=fstr,DEM=DEM)

;	plot
loadct,3
wvl=0.5*(ww(1:*)+ww)
plot,wvl,sp1,psym=10,title=atom,xtitle='Ang',ytitle='ph/s/cm^2/keV'
oplot,wvl,sp2,psym=1,col=150

end
