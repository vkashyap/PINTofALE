;+
;RADLOSS.PRO
;	Radiative loss function
;
;	use SPEX line emissivities from 1-50 AA, and CHIANTI from 50-3000 AA,
;	and add CIE continuum emissivities to it.
;
;vinay k. (Jul98)
;-

;	initialize SCAR keywords
if not keyword_set(ldir1) then ldir1='$SPEX'
if not keyword_set(ldir2) then ldir2='$CHIANTI'
if not keyword_set(cdir) then cdir='$CONT'
chifil='ioneq/mazzotta_etal.ioneq'
;
n_e=1e10	;[cm^-3]
wrl1=[1.,50.] & wrl2=[50.,3000.] & wrc=[1.,500.]
;
;abund=getabund('allen')
abund=getabund('grevesse et al.')
;abund=getabund('feldman et al.')
;abund=getabund('chromospheric')
fcs=1
;
lnlst=[	'0	'+strtrim(wrl2(0),2)+'-'+strtrim(wrl2(1),2)+'	'+ldir2,$
	'0	'+strtrim(wrl1(0),2)+'-'+strtrim(wrl1(1),2)+'	'+ldir1]
;(in apparent reverse order because of a quirk in CAT_LN, which has a
;DO-loop over the elements in the second part)
;lnlst=[	'0	'+strtrim(wrl1(0),2)+'-'+strtrim(wrl2(1),2)+'	'+ldir2]

;	lines
fls=rd_list(lnlst,/incieq,n_e=n_e,chifil=chifil)
nlin=n_elements(fls.WVL)
logT=fls.logT & nT=n_elements(logT)

;	continuum
f=rd_cont('cie',n_e=n_e,wrange=wrc,dbdir=cdir,abund=abund,fcstr=fcs)
ncon=n_elements(fcs.midWVL) & WW=fcs.WVL & dW=WW(1:*)-WW

;	P(T)
powl=dblarr(nT) & powc=powl & zz=fls.Z & zab=abund(zz-1)
for it=0,nT-1 do begin
  powl(it)=total( ( ((fls.LINE_INT))(it,*) ) * zab )
  powc(it)=total( ( ((fcs.CONT_INT))(it,*) ) * dW )
endfor
;
pow=powl+powc

;	compare with R-S?
powrs=dblarr(nT)
for it=0,nT-1 do begin
  thermalp,10.^(logT(it)),p,metals=1. & powrs(it)=p*1e23
endfor

plot,logt,powl+powc,/yl,ytitle='[1e-23 ergs cm^3/s]',xtitle='log(T)',$
	title='Radiative Loss Function',yr=[1e-4,2e2],/ys
oplot,logT,powl,psym=10,line=1 & oplot,logT,powc,psym=10,line=2
oplot,logt,powrs,line=3
oplot,[6,6.5],0.01*[1,1] & xyouts,6.6,0.01,'SCAR'
oplot,[6,6.5],0.006*[1,1],line=1 & xyouts,6.6,0.006,'lines'
oplot,[6,6.5],0.004*[1,1],line=2 & xyouts,6.6,0.004,'continuum'
oplot,[6,6.5],0.002*[1,1],line=3 & xyouts,6.6,0.002,'Raymond-Smith'

;save,file='radloss.sav',logT,pow,powl,powc,powrs,abund,n_e,wrl1,wrl2,wrc,$
;	chifil,ldir1,ldir2,cdir,lnlst

end
