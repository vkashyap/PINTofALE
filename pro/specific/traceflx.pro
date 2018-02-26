;+
;TRACEFLX.PRO
;	interactively specify DEM and compute passband flux from
;	spectral lines
;
;vinay kashyap (Apr 97)
;	updated (99Jun)
;-

;	initialize
restore,'trace.save',/verb			;read in TRACE data
nt=8 & dt=0.2 & tlog=findgen(nt)*dt+5.5		;temperature grid
tmin=min(tlog) & tmax=max(tlog)
logP=16						;cm^-3 K
n_e=1e9						;[cm^-3]
ftrace=dblarr(3)				;TRACE fluxes
abund=getabund('anders & grevesse')		;abundances
dbdir='$SPEX'
dbdir='$CHIANTI'				;line database
chifil='ioneq/arnaud_raymond_lmf.ioneq'		;use new CHIANTI ion balances

;	for 173 A passband
if not keyword_set(fe173) then begin
  w0=min(w173,max=w1)			;wavelength range
  if not keyword_set(fe173) then begin
    fe173=rd_line(elem,logP=logP,n_e=n_e,dbdir=dbdir,wrange=[w0,w1],wvl=wvl173,$
	logT=logT,Z=Z173,ion=ion173,jon=jon173)	;read in all the lines
    fe173=fold_ioneq(fe173,Z173,jon173,logt=logt,chifil=chifil)
  endif
  oo=where(logt ge tmin and logt le tmax)
  logt=logt(oo) & fe173=fe173(oo,*)
  mt=n_elements(logt) & mw173=n_elements(wvl173)
endif

;	for 195 A passband
if not keyword_set(fe195) then begin
  w0=min(w195,max=w1)			;wavelength range
  if not keyword_set(fe195) then begin
    fe195=rd_line(elem,logP=logP,n_e=n_e,dbdir=dbdir,wrange=[w0,w1],wvl=wvl195,$
	logT=logT,Z=Z195,ion=ion195,jon=jon195)	;read in all the lines
    fe195=fold_ioneq(fe195,Z195,jon195,logt=logt,chifil=chifil)
  endif
  oo=where(logt ge tmin and logt le tmax)
  logt=logt(oo) & fe195=fe195(oo,*)
  mt=n_elements(logt) & mw195=n_elements(wvl195)
endif

;	for 284 A passband
if not keyword_set(fe284) then begin
  w0=min(w284,max=w1)			;wavelength range
  if not keyword_set(fe284) then begin
    fe284=rd_line(elem,logP=logP,n_e=n_e,dbdir=dbdir,wrange=[w0,w1],wvl=wvl284,$
	logT=logT,Z=Z284,ion=ion284,jon=jon284)	;read in all the lines
    fe284=fold_ioneq(fe284,Z284,jon284,logt=logt,chifil=chifil)
  endif
  oo=where(logt ge tmin and logt le tmax)
  logt=logt(oo) & fe284=fe284(oo,*)
  mt=n_elements(logt) & mw284=n_elements(wvl284)
endif

done=0 & dem0=0.D*tlog+23.5
while done eq 0 do begin	;{get DEM, get fluxes
  dem0=demacs(tlog,dem0=dem0,group=group,igroup=igroup)
  dem=interpol(dem0,tlog,logt) & dem=10.D^(dem)
  f173=lineflx(fe173,logt,wvl173,Z173,DEM=DEM,abund=abund,$
	effar=a173,wvlar=w173)
  f195=lineflx(fe195,logt,wvl195,Z195,DEM=DEM,abund=abund,$
	effar=a195,wvlar=w195)
  f284=lineflx(fe284,logt,wvl284,Z284,DEM=DEM,abund=abund,$
	effar=a284,wvlar=w284)
  ftrace=[total(f173),total(f195),total(f284)] & ftrace=ftrace/max(ftrace)
  ;	print out the relative fluxes in the 173, 195, 284 A bands
  print,'' & print,string(ftrace) & print,''
  c1='x to stop, q to quit, any key to continue' & print,c1
  c1=get_kbrd(1)
  if strlowcase(c1) eq 'q' then done=1
  if strlowcase(c1) eq 'x' then stop
endwhile			;DONE}

end
