;+
;SOLARB_SPEC.PRO
;
;program for generating a "spectrum" for Solar-B
;							vinay kashyap
;-

;	initialize spectrum range &c.
w0=9.			;begin at 9 A
w1=19.			;end at 19 A
dw=2.			;bin size
nbin=long((w1-w0)/dw)	;number of bins

;	keywords
if not keyword_set(logP) then logP=16.		;n*T [cm^-3 K]
dbdir='/data/fubar/SCAR/emissivity/xinti'
chifil=1
logT=findgen(41)*0.05+5				;temperatures
;if not keyword_set(dem) then dem=fltarr(41)+2e14
if not keyword_set(abund) then abund=getabund('anders & grevesse')
loadct,3
fidgit,findgen(nbin)*dw+w0,randomu(seed,nbin)+1,dem=dem,logt=logt,$
	logP=logP,dbdir=dbdir,chifil=chifil,abund=abund
wdelete & wdelete

;	generate spectrum
spec=linespec(atom,wmin=w0,wmax=w1,nbin=nbin,ww=ww,$
	logP=logP,dbdir=dbdir,chifil=chifil,dem=dem,abund=abund)

;	plot
plot,findgen(nbin)*dw+w0+dw/2,spec,psym=10,xr=[w0,w1],/xs,$
	xtitle='[Ang]',ytitle='Relative Intensity',title='Solar B "spectrum"'

end
