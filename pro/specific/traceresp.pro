;+
;TRACERESP.PRO
;	T v/s Line Flux for given EM for given density
;
;vinay kashyap (Jun 99)
;-

;	initialize
if not keyword_set(tracesavfil) then tracesavfil='trace.save'
restore,tracesavfil,/verb			;read in TRACE data
	;{include filter transmission
t173=interpol(trans,lambda,w173) & a173=a173*t173
t195=interpol(trans,lambda,w195) & a195=a195*t195
t284=interpol(trans,lambda,w284) & a284=a284*t284
	;done}
nt=41 & dt=0.05 & tlog=findgen(nt)*dt+5		;temperature grid
if not keyword_set(n_e) then n_e=1e9		;[cm^-3]
if not keyword_set(EM) then EM=1e19		;[cm^-5]
f173=dblarr(nt) & f195=f173 & f284=f173		;fluxes [DN]
abund=getabund('anders & grevesse')		;abundances
dbdir='$CHIANTI'				;line database
chifil='ioneq/arnaud_raymond_lmf.ioneq'		;use new CHIANTI ion balances

pmulti=!p.multi & !p.multi=[0,1,3]
if not keyword_set(ylog) then ylog=0

  ;	for 173 A passband
  w0=min(w173,max=w1)			;wavelength range
  if not keyword_set(fe173) then begin
    fe173=rd_line(elem,n_e=n_e,dbdir=dbdir,wrange=[w0,w1],wvl=wvl173,$
	logT=logT,Z=Z173,ion=ion173,jon=jon173)	;read in all the lines
    fe173=fold_ioneq(fe173,Z173,jon173,logt=logt,chifil=chifil)
  endif
  mt=n_elements(logt) & mw=n_elements(wvl173)
  ;
  for i=0,nt-1 do begin			;{step through temperatures
    temp=tlog(i) & delt=logt-temp & ww=wvl173 & zz=Z173
    dtmn=min(abs(delt),imn)		;best match...
    ft=reform(fe173(imn,*))		;...and corresponding emissivities
    if dtmn ne 0 then begin		;(interpolate
      imnm=(imn-1)>0
      imnp=(imn+1)<(mt-1)
      w=[delt(imnm),delt(imn),delt(imnp)] & w=1./abs(w) & w=w/total(w)
      ft=w(0)*reform(fe173(imnm,*))+w(1)*ft+w(2)*reform(fe173(imnp,*))
    endif					;min(delt).NE.0)
    kilroy; was here.
    for j=0,mw-1 do begin
      f=lineflx(ft(j),temp,ww(j),zz(j),DEM=EM,abund=abund,$
	effar=a173,wvlar=w173) & f173(i)=f173(i)+f
    endfor
  endfor					;I=0,NT-1}
  ;
  plot,tlog,f173,title='TRACE [173A]; n!d!e!n='+strtrim(n_e,2),$
	xtitle='log!d10!n(T [K])',ytitle='[DN]',ylog=ylog,$
	subtitle='EM='+string(EM,'(g7.2)')+' [cm!u-5!n]'

  ;	for 195 A passband
  w0=min(w195,max=w1)
  if not keyword_set(fe195) then begin
    fe195=rd_line(elem,logP=logP,n_e=n_e,dbdir=dbdir,wrange=[w0,w1],wvl=wvl195,$
	logT=logT,Z=Z195,ion=ion195,jon=jon195)	;read in all the lines
    fe195=fold_ioneq(fe195,Z195,jon195,logt=logt,chifil=chifil)
  endif
  mt=n_elements(logt) & mw=n_elements(wvl195)
  ;
  for i=0,nt-1 do begin			;{step through temperatures
    temp=tlog(i) & delt=logt-temp & ww=wvl195 & zz=Z195
    dtmn=min(abs(delt),imn)		;best match...
    ft=reform(fe195(imn,*))		;...and corresponding emissivities
    if dtmn ne 0 then begin		;(interpolate
      imnm=(imn-1)>0
      imnp=(imn+1)<(mt-1)
      w=[delt(imnm),delt(imn),delt(imnp)] & w=1./abs(w) & w=w/total(w)
      ft=w(0)*reform(fe195(imnm,*))+w(1)*ft+w(2)*reform(fe195(imnp,*))
    endif					;min(delt).NE.0)
    kilroy; was here.
    for j=0,mw-1 do begin
      f=lineflx(ft(j),temp,ww(j),zz(j),DEM=EM,abund=abund,$
	effar=a195,wvlar=w195) & f195(i)=f195(i)+f
    endfor
  endfor					;I=0,NT-1}
  ;
  plot,tlog,f195,title='TRACE [195A]; n!d!ne='+strtrim(n_e,2),$
	xtitle='log!d10!n(T [K])',ytitle='[DN]',ylog=ylog,$
	subtitle='EM='+string(EM,'(g7.2)')+' [cm!u-5!n]'

  ;	for 284 A passband
  w0=min(w284,max=w1)
  if not keyword_set(fe284) then begin
    fe284=rd_line(elem,logP=logP,n_e=n_e,dbdir=dbdir,wrange=[w0,w1],wvl=wvl284,$
	logT=logT,Z=Z284,ion=ion284,jon=jon284)	;read in all the lines
    fe284=fold_ioneq(fe284,Z284,jon284,logt=logt,chifil=chifil)
  endif
  mt=n_elements(logt) & mw=n_elements(wvl284)
  ;
  for i=0,nt-1 do begin			;{step through temperatures
    temp=tlog(i) & delt=logt-temp & ww=wvl284 & zz=Z284
    dtmn=min(abs(delt),imn)		;best match...
    ft=reform(fe284(imn,*))		;...and corresponding emissivities
    if dtmn ne 0 then begin		;(interpolate
      imnm=(imn-1)>0
      imnp=(imn+1)<(mt-1)
      w=[delt(imnm),delt(imn),delt(imnp)] & w=1./abs(w) & w=w/total(w)
      ft=w(0)*reform(fe284(imnm,*))+w(1)*ft+w(2)*reform(fe284(imnp,*))
    endif					;min(delt).NE.0)
    kilroy; was here.
    for j=0,mw-1 do begin
      f=lineflx(ft(j),temp,ww(j),zz(j),DEM=EM,abund=abund,$
	effar=a284,wvlar=w284) & f284(i)=f284(i)+f
    endfor
  endfor					;I=0,NT-1}
  ;
  plot,tlog,f284,title='TRACE [284A]; n!d!ne='+strtrim(n_e,2),$
	xtitle='log!d10!n(T [K])',ytitle='[DN]',ylog=ylog,$
	subtitle='EM='+string(EM,'(g7.2)')+' [cm!u-5!n]'

!p.multi=pmulti

end
