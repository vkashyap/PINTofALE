;+
;SOLARB_SPEC.PRO
;
;program for generating Intensity v/s T responses for Solar-B
;"passbands"
;							vinay kashyap
;-

;	initialize
nt=41 & dt=0.05 & tlog=findgen(nt)*dt+5.5	;temperature grid
logP=16						;cm^-3 K
EM=1e19						;cm^-5
wmin=10 & dw=2. & nw=6 & w=findgen(nw)*dw+wmin+dw/2	;mid-bin wavelengths
w(nw-1)=w(nw-1)+1.
flx=dblarr(nw,nt)				;fluxes in 2A bins from 8-20A
abund=getabund('coronal')			;abundances
dbdir='/data/fubar/SCAR/emissivity/xinti'	;line database
chifil='ioneq/arr.ioneq'			;CHIANTI, but new Fe data
if not keyword_set(rdin) then rdin=intarr(nw)	;read in line database or no?

for iw=0,nw-1 do begin			;{for each wavelength "bin"
;for iw=nw-1,0,-1 do begin			;{for each wavelength "bin"
  if rdin(iw) eq 0 then begin		;(read in line info
    w_0=w(iw)-dw/4. & w_1=w(iw)+dw/4.
    tmp=rd_line(elem,logP=logP,dbdir=dbdir,wrange=[w_0,w_1],wvl=wtmp,$
	logT=logT,Z=Ztmp,ion=itmp,jon=jtmp)	;read in all the lines
    tmp=fold_ioneq(tmp,Ztmp,jtmp,logt=logt,chifil=chifil)
    rdin(iw)=1
    if iw eq 0 then begin & f0=tmp & w0=wtmp & z0=ztmp & i0=jtmp & endif
    if iw eq 1 then begin & f1=tmp & w1=wtmp & z1=ztmp & i1=jtmp & endif
    if iw eq 2 then begin & f2=tmp & w2=wtmp & z2=ztmp & i2=jtmp & endif
    if iw eq 3 then begin & f3=tmp & w3=wtmp & z3=ztmp & i3=jtmp & endif
    if iw eq 4 then begin & f4=tmp & w4=wtmp & z4=ztmp & i4=jtmp & endif
    if iw eq 5 then begin & f5=tmp & w5=wtmp & z5=ztmp & i5=jtmp & endif
  endif else begin			;)(RDIN(IW)=1
    if iw eq 0 then begin & tmp=f0 & wtmp=w0 & ztmp=z0 & jtmp=i0 & endif
    if iw eq 1 then begin & tmp=f1 & wtmp=w1 & ztmp=z1 & jtmp=i1 & endif
    if iw eq 2 then begin & tmp=f2 & wtmp=w2 & ztmp=z2 & jtmp=i2 & endif
    if iw eq 3 then begin & tmp=f3 & wtmp=w3 & ztmp=z3 & jtmp=i3 & endif
    if iw eq 4 then begin & tmp=f4 & wtmp=w4 & ztmp=z4 & jtmp=i4 & endif
    if iw eq 5 then begin & tmp=f5 & wtmp=w5 & ztmp=z5 & jtmp=i5 & endif
  endelse				;RDIN(IW)=1)
  mt=n_elements(logt) & mw=n_elements(wtmp)

  wvlar=findgen(11)*0.1*dw/2.+w(iw)-dw/2./2.
  effar=(wvlar-w(iw))^2/2./(dw/5.)^2 & effar=exp(-effar)

  for i=0,nt-1 do begin			;{step through temperatures
    temp=tlog(i) & delt=logt-temp & ww=wtmp & zz=ztmp
    dtmn=min(abs(delt),imn)		;best match...
    ft=reform(tmp(imn,*))		;...and corresponding emissivities
    if dtmn ne 0 then begin		;(interpolate
      imnm=(imn-1)>0
      imnp=(imn+1)<(mt-1)
      wt=[delt(imnm),delt(imn),delt(imnp)] & wt=1./abs(wt) & wt=wt/total(wt)
      ft=wt(0)*reform(tmp(imnm,*))+wt(1)*ft+wt(2)*reform(tmp(imnp,*))
    endif				;min(delt).NE.0)
    kilroy; was here.
    for j=0,mw-1 do begin
      ftmp=lineflx(ft(j),temp,ww(j),zz(j),DEM=EM,abund=abund,$
	effar=effar,wvlar=wvlar)
      flx(iw,i)=flx(iw,i)+ftmp
    endfor
  endfor				;I=0,NT-1}
endfor					;IW=0,NW-1}

;	plot
pmulti=!p.multi & !p.multi=[0,2,3]
for iw=0,nw-1 do plot,tlog,flx(iw,*),/yl,yr=[1e-4,1]*max(flx(iw,*)),$
	title='SolarB count rates ['+strtrim(w(iw),2)+' A]',$
	subtitle='EM='+string(EM,'(g7.2)')+' [cm!u-5!n]',ytitle='[ph/s/cm^2]',$
	xtitle='log!d10!n(T [K])'
!p.multi=pmulti

end
