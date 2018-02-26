;+
;PIPELINE_01.PRO
;
;1. reads in ASCA spectrum and associated RSP
;2. reads in initial DEM
;3. reads in line emissivities for specified wavelength range
;4. enters a keyboard controlled loop
;	A. reads in continuum emissivities
;	B. constructs line spectrum
;	C. constructs continuum spectrum
;	D. convolves with the ASCA response
;	E. adds the two
;	F. plots the computed and observed spectra
;	G. checks for user input to do one of the following:
;		a. change metallicity (type in)
;		b. change abundances (use SETABUND)
;		c. edit DEM (use DEMACS)
;		d. choose channels (use PICKRANGE)
;		e. changes plotting styles
;		f. quit
;
;vinay kashyap (Jul97)
;-

;	user controlled variables
datadir='/data/drake6/vinay/22017000/gof'	;ASCA data directory
asca='s0b1'					;ASCA data file root
savdem='/data/drake7/vinay/SCAR/emissivity/tmp/dem.sav'	;save file w. DEM(T)
ldbdir='$SCAR'				;line DB directory
cdbdir='/data/drake7/vinay/SCAR/emissivity/emisscie'	;continuum emiss.
chifil='ioneq/arnaud_rothenflug.ioneq'			;ion balance
;NOTE: continua are computed with Arnaud & Rothenflug, not Arnaud & Raymond
abund=getabund('anders & grevesse')		;abundances
logP=16.		;log(pressure = nT [cm^-3 K])
eden=1e10		;electron density [cm^-3] (set to 0 to use logP)
croot='angr'		;prefix for continuum emissivities files
metal=1.		;initial metallicity
setplot='el'		;plot counts v/s energy on log scale
loadct,3		;color table (data in white, model in reddish)
window,0		;for data
window,1		;for model
window,2		;for DEM

;111111111111111111111111111111111111111

;	ASCA spectrum
phafil=datadir+'/'+asca+'.pha'
;	read spectrum
spec=mrdfits(phafil,1,hspec)
chan=spec.CHANNEL+1 & counts=spec.COUNTS & cterr=1.+sqrt(COUNTS+0.75)
exptim=sxpar(hspec,'EXPOSURE')			;exposure time
rspfil=strtrim(sxpar(hspec,'RESPFILE'),2)	;.rsp or .rmf
arffil=strtrim(sxpar(hspec,'ANCRFILE'),2)	;.arf
;	read RSP or RMF+ARF
if strlowcase(arffil) ne 'none' then begin
  rsp=marfrmf(datadir+'/'+arffil,datadir+'/'+rspfil,rstr)
endif else rsp=rdresp(datadir+'/'+rspfil,rstr)
;unpack RSTR
nnrg=rstr.NNRG & elo=rstr.ELO & ehi=rstr.EHI
nchan=rstr.NCHAN & emn=rstr.EMN & emx=rstr.EMX
;
wset,0 & plot,chan,counts/exptim,xtitle='PHA',ytitle='Counts s!u-1!n'
wset,1 & plot,chan,fold_resp(rsp,lindgen(nchan)),xtitle='PHA',ytitle='response'

;222222222222222222222222222222222222222

;	DEM
restore,savdem & dem0=DEM & logt0=logt
wset,2 & plot,logt0,dem0,xtitle='log!d10!n(T [K])',ytitle='DEM [cm!u-5!n]'

;333333333333333333333333333333333333333

;	line emissivities
;
wrange=[10.,0.2] & wrange=12.3985/wrange
;
;read emissivities
if n_tags(fstr) eq 0 then begin		;if FSTR is read in, don't redo
  fstr=1
  flin=rd_line(atom,logP=logP,n_e=eden,wrange=wrange,dbdir=ldbdir,fstr=fstr)
  if n_elements(flin) le 1 then stop,'type .con to continue anyway'
  ;
  ;fold in ion balance
  Z=fstr.Z & ion=fstr.ION & jon=fstr.JON logT=fstr.LOGT
  flin=fold_ioneq(flin,Z,jon,logT=logT,chifil=chifil)
  fstr.(0) = flin
endif
;
;get bin indices of wavelengths
ee=fstr.wvl & ee=12.3985/ee & iel=binersp(ee,bstr=rstr)
;
;get best guess at DEM
wset,2 & dem=demacs(logT,dem0=dem0,logt0=logt0,group=group,igroup=igroup)
dem0=dem & logt0=logt

;444444444444444444444444444444444444444

;	loop until asked to quit
ok='' & first='y' & ochan=lindgen(nchan)
while ok ne 'q' do begin			;{

  ;A: read in continuum emissivities
  fcstr=1
  fcon=rd_cont(croot,metal=metal,wrange=wrange,dbdir=cdbdir,fcstr=fcstr)
  if n_elements(fcon) le 1 then begin
    print,'try again' & goto,reask
  endif
  if first eq 'y' then begin
    ww=fcstr.midwvl & nrg=6.626176e-27*2.9979e10*1e8/ww
    ee=12.3985/ww & iec=binersp(ee,bstr=rstr)
    reH=fcstr.reH
    wvl=fcstr.wvl & dw=wvl(1:*)-wvl
    nt=n_elements(logT) & dlogT=total(logT(1:*)-logT)/(nt-1)
    first='n'
  endif

  ;B: fold in DEM and abundances with line emissivities
  abnd=abund & abnd(2:*)=abund(2:*)*metal
  wvl=fstr.wvl & nwvl=n_elements(wvl)
  flx=lineflx(flin,logT,wvl,Z,DEM=DEM,abund=abnd)
  wset,1 & mflx=max(flx)*1e-7
  plot,[0],xtitle='[keV]',ytitle='[ph/s/cm^2]',xr=[0.2,10],$
	yr=[mflx,max(flx)],/yl,/nodata
  for i=0,nwvl-1 do oplot,[1,1]*12.3985/abs(wvl(i)),[mflx,flx(i)],col=150

  ;C: fold in DEM with continuum spectrum
  fcx=0.*ww
  for it=0,nt-1 do begin
    fcx=fcx+DEM(it)*fcon(it,*)*dlogT*dw/nrg/reH(it)^2
    ;[1e-23 ph/s/cm^2]=[cm^-5]*[1e-23 ergs cm^3/s/A]*[A]/[ergs]/[(n_e/n_H)^2]
  endfor
  fcx=fcx/1e23		;[ph/s/cm^2]
  oplot,12.3985/ww,fcx,psym=10,col=200

  ;D: fold through ASCA response
  phl=fold_resp(rsp,iel,flx)
  phc=fold_resp(rsp,iec,fcx)

  ;E: add the two and compute chi^2
  pha=phl+phc
  chi2=total((pha(ochan)-counts(ochan)/cterr(ochan))^2)

  replot: ;F: plot
  if setplot eq 'e' or setplot eq 'el' then begin
    ee=rstr.emn & xt='[keV]'
    if setplot eq 'el' then xlog=1 else xlog=0
  endif
  if setplot eq 'c' or setplot eq 'cl' then begin
    ee=chan & xt='PHA'
    if setplot eq 'cl' then xlog=1 else xlog=0
  endif
  wset,0
  plot,ee,counts,psym=10,/yl,yr=[1,max(counts)],xtitle='Energy',$
	ytitle='Counts',title=phafil,subtitle='Chi2:'+strtrim(chi2,2),$
	xlog=xlog
  for i=0,nchan-1 do oplot,ee(i)*[1,1],counts(i)+cterr(i)*[1,-1]
  oplot,ee,pha,psym=10,col=150

  reask: ;G:	loop control
  print,'q to quit, m to change metallicity, d to edit DEM, r to renorm'
  print,'c to choose channels, a to change abundances, p for plot settings'
  ok=get_kbrd(1)
  case ok of
    'q':	;quit
    'm': read,prompt='new metallicity? ',metal
    'a': begin
      ab=setabund(init=abnd)
      metal=ab(25)/getabund('anders & grevesse',elem='Fe')
      print,'choosing new metallicity of',metal
    end
    'c': begin
      ochan=pickrange(rstr.emn,counts,dynrng=6)
      if ochan(0) eq -1 then ochan=lindgen(nchan)
    end
    'd': begin
      wset,2 & dem=demacs(logT,dem0=dem0,logT0=logT0,group=group,igroup=igroup)
      dem0=dem & logt0=logT
    end
    'r': begin
      rat=total(pha(ochan))/total(counts(ochan))
      if rat gt 0 then dem=dem/rat
      wset,2
      plot,logt,dem,xtitle='log!d10!n(T [K])',ytitle='DEM [cm!u-5!n]',/yl
    end
    'p': begin
      print,'e for counts v/s energy, c for counts v/s PHA'
      print,'el and cl for log scaled X-axes'
      c1='' & read,prompt='plotting style? ['+setplot+'] ',c1
      c2=strlowcase(strmid(c1,0,1))
      if c2 ne '' and c2 ne ' ' then setplot=c1
      goto, replot
    end
    else: goto, reask
  endcase

endwhile					;}

end
