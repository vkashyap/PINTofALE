;+
;MK_CONT_CIE.PRO
;	IDL program to call WRT_CONT_CIE repeatedly for various densities
;	to generate a database of continuum emission
;
;vinay kashyap (Dec97, based on MK_LCOOL_CHIANTI.PRO)
;-

abund=getabund('anders & grevesse')
eden=findgen(6)+8 & pres=findgen(8)+13
cieZ=['He','C','N','O','Ne','Na','Mg','Al','Si','S','Ar','Ca','Fe','Ni']
wrange=[1.,3000.]
wrange=[0.5,3000.]	;VK 12Dec06
nbin=-2000L		;VK 12Dec06
outdir='/data/fubar/SCAR/emissivity/contnew'
root='cie'
ciedir='/data/fubar/SCAR/CIE'

;	first, at constant electron density
for i=0,n_elements(eden)-1 do begin
  n_e=10.^(eden(i))
  wrt_cont_cie,cieZ,pr,logT,wvl,reH,n_e=n_e,wrange=wrange,$
	outdir=outdir,root=root,nbin=nbin,abund=abund,ciedir=ciedir
endfor

;	next, at constant pressure
for i=0,n_elements(pres)-1 do begin
  pr=10.^(pres(i))
  wrt_cont_cie,cieZ,pr,logT,wvl,reH,n_e=0,wrange=wrange,$
	outdir=outdir,root=root,nbin=nbin,abund=abund,ciedir=ciedir
endfor

end
