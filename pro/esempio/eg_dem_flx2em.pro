;+
;SERTS_DEM_FLX2EM.PRO
;
;program to call DEM_FLX2EM for SERTS data
;				-vinay kashyap
;-

;	first restore the appropriate SERTS data and IDs
restore,'act93.sav'		;AR 7563, 99.4s exp, 1993 Aug 17
;restore,'qut93.sav'		;202.3s exp, 1993 Aug 17

;	unpack
x=wvl & y=flx				;observed spectrum {[A],[ph/s/cm^2]}
NH=1e10					;column density [cm^-2]
exptim=1.				;exposure time [s]
xwvl=allid.wvl & nw=n_elements(xwvl)	;matched wavelengths
;SIG, LOGP, DBDIR, ABUND, CHIFIL	;also restored
dbdir='$SCAR'
abund=getabund('grevesse et al.')	;use different version
chifil='ioneq/arnaud_raymond.ioneq'	;use new version

;	initialize
if n_elements(onlypot) eq 0 then onlypot=1	;plot only Pottasch points?
rom=['I','II','III','IV','V','VI','VII','VIII','IX','X']
rom=[rom,'X'+rom,'XX'+rom]			;roman numerals from 1-30
atom=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',$
	'Ni','Cu','Zn']				;elements from 1-30
defEM=1e28					;default EM
level=0.1					;FullWidth@(LEVEL*MAX)
loadct,3
allEM=0*y & allZ=allEM & allI=allZ
eden=10.^(findgen(7)+8)

;	set up plot scales
plot,[0],xr=[7,14],yr=[5e26,1e29],/ylog,xtitle='log(N!de!n)',$
	ytitle='EM [cm!u-5!n]',/xs,/ys

for i=0,nw-1 do begin			;{for each match
  fx=y(where(xwvl(i) eq x))/exptim	;observed flux [ct/s/cm^2]
  sigma=sig(where(xwvl(i) eq x))/exptim	;error on FX
  snr=fx/sigma				;S/N
  wvl=(allid.(i+1).wvl)			;wavelengths
  Z=(allid.(i+1).z)			;atomic numbers
  ion=(allid.(i+1).ion)			;ionic states
  labl=(allid.(i+1).labl)		;labels
  nid=n_elements(wvl)			;number of IDs
  ;
  emden=dem_flx2em(fx,wvl,Z,ion,eden=eden,$
	abund=abund,NH=NH,defEM=defEM,dbdir=dbdir,chifil=chifil,level=level)
  ;
  print,i,wvl,fx(0),sigma(0),atom(z(0)-1)+rom(ion(0)-1)
  oplot,alog10(eden),emden,col=200
  xyouts,alog10(eden(0)),emden(0),atom(z(0)-1)+rom(ion(0)-1),col=150,align=-0.5
  ;
  ;c1='q to stop, any key to continue' & print,c1 & c1=get_kbrd(1)
  ;if c1 eq 'q' then stop
endfor					;i=0,nw-1}

end
