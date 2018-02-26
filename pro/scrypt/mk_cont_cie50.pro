;+
;MK_CONT_CIE50.PRO
;	IDL program to call WRT_CONT_CIE repeatedly for various densities
;	to generate a database of continuum emission
;
;vinay kashyap (Dec97, based on MK_LCOOL_CHIANTI.PRO)
;-

abund=getabund('anders & grevesse')
eden=findgen(6)+8 & pres=findgen(8)+13
cieZ=['He','C','N','O','Ne','Na','Mg','Al','Si','S','Ar','Ca','Fe','Ni']
wrange=[12.3985/50.,500.]
nbin=-600L
outdir='/data/fubar/SCAR/emissivity/cont'
root='cie50'
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
