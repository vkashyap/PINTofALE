;+
;ADDLINCON.PRO
;	example script to demonstrate adding lines to continuum spectra
;
;subroutines
;	GETABUND [SETABUND]
;	RD_CONT [SYMB2ZION [LAT2ARAB]
;	LINEFLX [GETABUND, WHEE]
;	LINESPEC [RD_LINE [KILROY, SYMB2ZION [LAT2ARAB],
;	  FOLD_IONEQ [WHEE, GETLOCMAX, RD_IONEQ [READ_IONEQ (CHIANTI)]]],
;	  LINEFLX [WHEE]
;	  HASTOGRAM [KILROY]]
;	REBINW [KILROY]
;
;vinay k (Jan 98)
;-

;	define keywords
if not keyword_set(ldbdir) then ldbdir='$SPEX'
if not keyword_set(cdbdir) then cdbdir='$CONT'
if not keyword_set(chifil) then chifil=$
	'ioneq/arnaud_rothenflug_lmf.ioneq'
if not keyword_set(n_e) then n_e=1e9
if not keyword_set(abund) then abund=getabund('anders & grevesse')
if not keyword_set(wmin) then wmin=1.	;Ang
if not keyword_set(wmax) then wmax=100.	;Ang
if not keyword_set(nbin) then nbin=1000L
if not keyword_set(flstr) then flstr=1
if not keyword_set(fcstr) then fcstr=1
if not keyword_set(noph) then noph=0	;[erg/...]
if not keyword_set(keV) then keV=0	;[.../Ang]

;	any changes?
if not keyword_set(dbdirl) then dbdirl=ldbdir
if not keyword_set(dbdirc) then dbdirc=cdbdir
if not keyword_set(filchi) then filchi=chifil
if not keyword_set(eden) then eden=n_e
if not keyword_set(oldwmin) then oldwmin=wmin
if not keyword_set(oldwmax) then oldwmax=wmax
if not keyword_set(oldnbin) then oldnbin=nbin

;	in case we have changed the ion balance
if filchi ne chifil then flstr=1

;	in case we have changed the database directory
if ldbdir ne dbdirl then flstr=1
if cdbdir ne dbdirc then fcstr=1

;	in case we have changed the density
if eden ne n_e then begin
  flstr=1 & fcstr=1
endif

;	in case we have changed the wavelength range
if keyword_set(ww) and keyword_set(oldww) then begin
  nww=n_elements(ww) & mww=n_elements(oldww)
  if nww ne mww then begin
    if nww gt 1 then begin
      wmin=min(ww,max=wmax) & flstr=1 & fcstr=1 & ww=0
    endif
  endif else begin
    dw=total(abs(ww-oldww))
    if dw gt 1e-5 then begin
      wmin=min(ww,max=wmax) & flstr=1 & fcstr=1 & ww=0
    endif
  endelse
endif
if keyword_set(oldwmin) then begin
  if oldwmin ne wmin or oldwmax ne wmax or oldnbin ne nbin then begin
    flstr=1 & fcstr=1 & ww=0
  endif
endif

if keyword_set(keV) then wrange=12.3985/[wmax,wmin] else wrange=[wmin,wmax]

help,ldbdir,cdbdir,chifil,n_e,abund,wmin,wmax,nbin,flstr,fcstr,noph,kev,dem

;	get continuum emissivities [1e-23 erg cm^3/s/Ang]
conem=rd_cont('cie',n_e=n_e,wrange=wrange,dbdir=cdbdir,abund=abund,$
	wvl=wvl,logT=logT,fcstr=fcstr)

;	get continuum spectrum [(ph|erg)/cm^2/s/(Ang|keV)]
if not keyword_set(dem) then begin
  message,'DEM not defined.  Assuming constant DEM',/info
  dem=dblarr(n_elements(logT))+1e19
endif
if n_elements(dem) ne n_elements(logT) then message,$
	'push DEM into same grid as logT first! (e.g., using MK_DEM)'
wave=0.5*(wvl(1:*)+wvl)
conspec=lineflx(conem,logT,wave,dem=dem,noph=noph)	;[(ph|erg)/cm^2/s/A]
if keyword_set(keV) then begin
  dw=wvl(1:*)-wvl & conspec=conspec*dw			;[../A] -> [../bin]
  wvl=12.3985/wvl & de=abs(wvl(1:*)-wvl) & conspec=conspec/de	; -> [../keV]
endif

;	make line spectrum [(ph|erg)/cm^2/s/(Ang|keV)]
linsp=linespec(atom,wmin=wmin,wmax=wmax,nbin=nbin,wlog=wlog,ww=ww,kev=kev,$
	fstr=flstr, n_e=n_e,dbdir=ldbdir,chifil=chifil,dem=dem,abund=abund,$
	noph=noph)

;	rebin continuum spectrum
consp=rebinw(conspec,wvl,ww)

;	add the two
spec=linsp+consp

;	plot
loadct,3
if keyword_set(keV) then xt='[keV]' else xt='[Ang]'
if keyword_set(keV) and keyword_set(noph) then yt='[ph/s/cm^2/keV]'
if not keyword_set(keV) and keyword_set(noph) then yt='[ph/s/cm^2/Ang]'
if keyword_set(keV) and not keyword_set(noph) then yt='[erg/s/cm^2/keV]'
if not keyword_set(keV) and not keyword_set(noph) then yt='[erg/s/cm^2/Ang]'
plot,ww,spec,/yl,yr=[1e-3,2]*max(spec),psym=10,xtitle=xt,ytitle=yt
oplot,ww,consp,psym=10,col=150
oplot,ww,linsp,psym=10,col=100

;	what if we want to change the wavelength range?
oldwmin=wmin & oldwmax=wmax & oldnbin=nbin & oldww=ww

;	what if we want to change the density?
eden=n_e

;	what if we want to change the database directory?
dbdirl=ldbdir & dbdirc=cdbdir

;	what if we want to change the ion balance?
filchi=chifil

end
