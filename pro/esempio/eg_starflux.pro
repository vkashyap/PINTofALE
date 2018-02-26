;+
;EG_STARFLUX.PRO
;	example program to call STARFLUX
;
;vinay kashyap (1999May)
;-

message,'this program computes PSPC count rates',/info

;	initialize
if not keyword_set(n_e) then n_e=1e9	;[cm^-3]	e-density
if not keyword_set(ldir) then ldir='$SPEX'
if not keyword_set(cdir) then cdir='$CONT'
if not keyword_set(conroot) then conroot='cie'
if not keyword_set(abund) then abund=getabund('grevesse et al')
if not keyword_set(chifil) then chifil='ioneq/arnaud_raymond_lmf.ioneq'
if not keyword_set(emin) then emin=0.2	;[keV]	passband
if not keyword_set(emax) then emax=2.0	;[keV]	passband
if not keyword_set(VEM) then VEM=1d50	;[cm^-3]	emission measure
if not keyword_set(dist) then dist=10.	;[pc]	distance to star
if not keyword_set(NH) then NH=5e19	;[cm^-2]	H-column
if not keyword_set(fH2) then fH2=0.26	;frac H2
if not keyword_set(He1) then He1=NH/10.	;[cm^-2]	He I-column
if not keyword_set(HeII) then HeII=NH/100.	;[cm^-2]	He II-column
wmin=12.3985/emax & wmax=12.3985/emin	;WRANGE

;	read in effective areas
rd_pimms_file,get_pimms_file('ROSAT','PSPC',special='open'),pspc_a,pspc_e
pspc_w=12.3985/pspc_e

;	read in emissivities
if n_tags(lstr) eq 0 then $
  lstr=rd_list('ALL:'+strtrim(wmin,2)+'-'+strtrim(wmax,2)+':'+ldir,sep=':',$
	/incieq,chifil=chifil,n_e=n_e)
cstr=1 & tmp=rd_cont(conroot,n_e=n_e,wrange=[wmin,wmax],dbdir=cdir,$
	abund=abund,fcstr=cstr)

;	compute count rates
if n_elements(tlog) eq 0 then tlog=findgen(31)*0.1+5.
message,'calling STARFLUX',/info
flx=starflux(lstr,cstr,tlog=tlog,DEM=VEM,abund=abund,wrange=[wmin,wmax],$
	NH=NH,fH2=fH2,He1=He1,HeII=HeII,effar=pspc_a,wvlar=pspc_w,$
	dist=dist,lflux=lflux,cflux=cflux)
;[ct/s] == [ergs cm^3/s]*[ph/ergs]*[EM:cm^-3]*[DIST:cm^-2]*[AREA:cm^2]

;	plot
plot,tlog,flx,xtitle='log!d10!n(T [K])',ytitle='[ct/s]',psym=-1,$
	title=ldir+' ; '+cdir
oplot,tlog,lflux/4./!dpi/dist^2/(3.1d18)^2,line=1
oplot,tlog,cflux/4./!dpi/dist^2/(3.1d18)^2,line=2

end
