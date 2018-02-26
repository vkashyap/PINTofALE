;+
;EG_LINESPEC.PRO
;
;example program for calling LINESPEC
;							vinay kashyap
;-

fstr=1		;to avoid calling RD_LINE many times
w0=1. & w1=2.	;wavelength ranges in keV
n_e=1e9		;e-density
dbdir='$SPEX'	;line DB directory
nbin=1000	;number of bins
atom='Ni'	;stick to Nickel

;	get spectrum in [1,2] keV range [ph/s/keV]
spe=linespec(atom,wmin=w0,wmax=w1,nbin=nbin,ww=we,/kev,fstr=fstr,$
	n_e=n_e,dbdir=dbdir)

w0=12.3985/w0 & w1=12.3985/w1	;convert to A
;	get same spectrum [ph/s/A] in same range again
spw=linespec(atom,wmin=w0,wmax=w1,nbin=nbin,ww=ww,fstr=fstr,$
	n_e=n_e,dbdir=dbdir)

dwe=we(1:*)-we & dww=ww(1:*)-ww			;delta_keV and delta_Ang
nrg=0.5*(we+we(1:*)) & wvl=0.5*(ww+ww(1:*))	;mid-bin keV and Ang

;	plot
loadct,3
plot,nrg,spe,psym=10,title=atom,xtitle='keV',ytitle='ph/s/cm^2/keV'
;	convert [ph/s/cm^2/A] to [ph/s/cm^2/keV], Ang to keV, and overplot
oplot,12.3985/wvl,spw*dww/dwe,psym=1,col=150

end
