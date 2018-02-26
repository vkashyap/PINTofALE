;+
;MK_LINES.PRO
;	IDL script to generate a line list from a test DEM
;	that may be used to exercise MCMC_DEM.PRO
;
;vinay kashyap
;-

;	initialize
root='/data/drake7/vinay/SCAR/emissivity/Sol97/sol4'	;save file prefix
logP=16							;pressure [cm^-3 K]
dbdir='/data/drake7/vinay/SCAR/emissivity/emiss/'	;line database
chifil='ioneq/arr.ioneq'				;ion balance file
abund=getabund('anders & grevesse')			;abundances
fstr=1							;store in structure
wmin=300. & wmax=400.					;wavelength range [Ang]
norm=1.							;[cm^2 s]
if n_elements(read_db) eq 0 then read_db=1						;read database?
maxline=50						;lines to keep
maxulim=00						;non-detections to keep
if not keyword_set(dem0) then begin			;initial DEM
  nt0=11 & dt0=0.2 & tmin0=5.
  logt0=findgen(nt0)*dt0+tmin0
  dem0=dblarr(nt0)+11.2
endif

if read_db eq 1 then begin
  ;	read in all lines in this range
  emis=rd_line(atom,logP=logP,wrange=[wmin,wmax],dbdir=dbdir,fstr=fstr,$
	/desig,/econf)

  ;	unpack FSTR
  logT=fstr.logt & nt=n_elements(logt)
  wvl=fstr.wvl & nw=n_elements(wvl)
  Z=fstr.Z
  ion=fstr.ion
  jon=fstr.jon
  desig=fstr.desig
  econf=fstr.config

  ;	how many to keep?
  print,'there are '+strtrim(nw,2)+' lines in specified wavelength range'
  ndet=(maxline<nw) & nul=(maxulim<(nw-ndet))>0
  print,"There'll be"+string(ndet)+' detections and'+string(nul)+' upper limits'
  nlines=ndet+nul

  ;	fold in ion-balance
  emis=fold_ioneq(emis,Z,jon,logT=logT,chifil=chifil)

  ;	find temperature at which each line has maximum contribution
  tmax=fltarr(nw) & emax=dblarr(nw)
  for iw=0,nw-1 do begin
    emax(iw)=max(reform(emis(*,iw)),imx) & tmax(iw)=logT(imx)
  endfor
  ;and some leetle stuff for SERTS...
  oz=where(Z eq 26 and ion gt 17,moz)
  if moz gt 0 then begin
    ul_line=emis(*,oz) & ul_wvl=wvl(oz) & ul_Z=Z(oz) & ul_ion=ion(oz)
  endif

  ;	don't need to go through this anymore
  read_db=0
endif

;	define DEM
dem=demacs(logt,dem0=dem0,logt0=logt0,group=group,igroup=igroup)
dem0=interpol(dem,logT,logt0)
dem=10.D^(dem)

;	get line fluxes
fx=lineflx(emis,logT,wvl,Z,DEM=dem,abund=abund)		;[ph/s/cm^2]
fx=fx*norm						;[ph]

;;	sort according to temperature maxima
;oo=sort(tmax) & flx=fx(oo) & line=emis(*,oo)
;ww=wvl(oo) & zz=Z(oo) & iion=ion(oo) & dlabl=desig(*,oo) & elabl=econf(*,oo)

;	sort the fluxes
oo=reverse(sort(fx)) & flx=fx(oo) & line=emis(*,oo)
ww=wvl(oo) & zz=Z(oo) & iion=ion(oo) & dlabl=desig(*,oo) & elabl=econf(*,oo)

;	throw away the weak lines
ii=lindgen(nlines)
flx=flx(ii) & ww=ww(ii) & zz=zz(ii) & iion=iion(ii) & line=line(*,ii)
ulim=lonarr(nlines) & if nul gt 0 then ulim(ndet:*)=1
sigma=1+sqrt(abs(flx)+0.75)
if nul gt 0 then sigma(ndet:*)=flx(ndet:*)/3.		;approx. 3-sigma ULs

;	now save stuff
c1='' & read,prompt='save to files with prefix ['+root+']? ',c1
if c1 eq '' or c1 eq ' ' or strlowcase(strmid(c1,0,1)) eq 'y' then $
	c1=root else root=c1

demfil=root+'_dem.sav' & blind=root+'.sav'
print,'saving all variables to file: '+demfil
save,file=demfil				;save everything
print,'saving fluxes to file: '+blind
save,file=blind,flx,ww,zz,iion,ulim,sigma,line,abund,logt

;	plot spectrum
plot,abs(ww),flx,psym=1,xtitle='[Ang]',ytitle='[ph]',title=demfil
for i=0,ndet-1 do oplot,abs(ww(i))*[1,1],[0,flx(i)]
for i=0,nul-1 do oplot,abs(ww(ndet+i))*[1,1],[0,flx(ndet+i)],line=1

;	verify
tmp=pred_flx(line,logt,ww,dem,Z=zz,fobs=flx,sigma=sigma)

end
