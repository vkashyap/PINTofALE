;+
;EG_FLX2EM.PRO
;
;example program to call FLX2EM
;				-vinay kashyap
;-

;	first call LINEID and get the IDs
restore,'eg_lineid.sav'		;cf. EG_LINEID.PRO

;	unpack
x=lamda & y=spec			;the observed spectrum
NH=1e18					;column density [cm^-2]
exptim=100000.				;exposure time [s]
;also restored: LOGP, DBDIR, SWEA, SW
xwvl=lid.wvl & nw=n_elements(xwvl)	;matched wavelengths
abund=getabund('anders & grevesse')	;abundances

;	initialize
rom=['I','II','III','IV','V','VI','VII','VIII','IX','X']
rom=[rom,'X'+rom,'XX'+rom]			;roman numerals from 1-30
atom=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',$
	'Ni','Cu','Zn']				;elements from 1-30

;	set up plot scales
plot,[0],xr=[5,8],yr=[1e10,1e20],/ylog,xtitle='logT',ytitle='EM [cm!u-5!n]',$
	title='logP='+strtrim(logP,2)+' [cm!u-3!n K]',/xs,/ys

;	get EMs
for i=0,nw-1 do begin			;{for each match
  fx=y(where(xwvl(i) eq x))/exptim	;observed flux [ct/s] (to within
					;factors of sqrt(2 pi fwhm))
  wvl=(lid.(i+1).wvl)			;wavelengths
  Z=(lid.(i+1).z)			;atomic numbers
  ion=(lid.(i+1).ion)			;ionic states
  nid=n_elements(wvl)			;number of IDs
  skiprd=0 & lfx=0			;don't use output of previous call

  ;	call FLX2EM
  em=flx2em(fx,wvl,Z,ion,lfx,logT=logT,abund=abund,NH=NH,skiprd=skiprd,$
	logP=logP,dbdir=dbdir,effar=swea,wvlar=sw)

  ;	plot
  for j=0,nid-1 do begin
    tmp=em(*,j) & oo=where(tmp gt 0,moo)
    if moo gt 0 then begin
      y0=min(tmp(oo),j0) & x0=logT(oo(j0))
      oplot,logT,tmp & xyouts,x0,0.9*y0,atom(z(j)-1)+rom(ion(j)-1)
    endif
  endfor

  c1='q to stop, any key to continue' & print,c1 & c1=get_kbrd(1)
  if c1 eq 'q' then stop

endfor					;i=0,nw-1}

end
