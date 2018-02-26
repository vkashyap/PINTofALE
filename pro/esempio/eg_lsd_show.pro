;+
;EG_LSD_SHOW.PRO
;
;	example program to find density sensitive lines in given wavelength
;	range with LSD.PRO and display it with SHOW_LINE.PRO
;
;vinay kashyap
;-

;	initialize for LSD
if n_elements(wrange) lt 2 then wrange=[80.,200.]
if not keyword_set(edens) then edens=[1e8,1e13]
if not keyword_set(elem) then elem=0
if not keyword_set(lindir) then lindir='$CHIANTI'
if not keyword_set(ceiling) then ceiling=0.5
if not keyword_set(floor) then floor=0.
if not keyword_set(chifil) then chifil='ioneq/arnaud_raymond.ioneq'
if n_elements(abund) lt 30 then abund=getabund('anders & grevesse')

;	initialize for SHOW_LINE
if not keyword_set(lambda) then begin
  print,'need wavelenghts for input spectrum in LAMBDA'
endif
if n_elements(spec) ne n_elements(lambda) then begin
  print,'need input spectrum SPEC(LAMBDA)'
endif
;
if n_tags(linstr) gt 0 then begin
  fx=lineflx(linstr.LINE_INT,linstr.LOGT,linstr.WVL,linstr.Z,$
	DEM=DEM,abund=abund,effar=effar,wvlar=wvlar,noph=noph,kev=kev)
  lWVL=linstr.WVL & lZ=linstr.Z & lION=linstr.ION
endif else begin
  lWVL=wrange
endelse
;
if not keyword_set(order) then order='1'

;	call LSD
ld=lsd(wrange,fluxes,wvls,elem=elem,edens=edens,dbdir=dbdir,ceiling=ceiling,$
	floor=floor,ratmax=ratmax,flxmax=flxmax,chifil=chifil,$
	DEM=DEM,abund=abund,effar=effar,wvlar=wvlar,noph=noph,kev=kev)

;	call SHOW_LINE
if ld(0) ne '' then begin
  show_line,lWVL,fx,Z=lZ,ion=lION,lambda=lambda,spec=spec,order=order,$
	wmark=wvls,fmark=flxmax,lmark=' '+ld,oxr=oxr,oyr=oyr,$
	sep=sep,squish=squish,dynrng=dynrng,markp=markp,markc=markc,$
	marko=marko,xtitle=xtitle,ytitle=ytitle
endif

end
