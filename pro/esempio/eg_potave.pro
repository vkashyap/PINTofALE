;+
;EG_POTAVE.PRO
;	example program to call POTAVE
;
;vinay kashyap (Jun97)
;-

loadct,3
;	initialize keywords
dbdir='$SPEX'						;line DB directory
logP=16.						;P=n*T [cm^-3 K]
n_e=0.							;[cm^-3]
abund=getabund('grevesse et al.')			;abundances
chifil='ioneq/arnaud_raymond.ioneq'			;ion balance
ewt=0							;error wtd. averages?

;	read in SERTS data
serts='act93'
;serts='qut93'
filnam='serts_'+serts+'.tab'
rd_sertstab,filnam,elem,wvl,flx,sigma,fwhm,fwhmsig
nw=n_elements(wvl)

;	initialize arrays
atom=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',$
	'Ni','Cu','Zn']			;elements from 1-30
rom=['I','II','III','IV','V','VI','VII','VIII','IX','X']
rom=[rom,'X'+rom,'XX'+rom,'XXX'+rom]	;roman numerals from 1-40
;	these are "density and temperature independent" lines (I suppose
;	Brosius means that the appropriate line ratios are d & T indep.)
if serts eq 'act93' then xdxt_w=$
	[274.168, 284.135, 291.961, 292.773, 296.108, 311.783, 312.171,$
	312.569, 313.744, 314.358, 315.029, 316.223, 317.018, 319.852,$
	320.568, 320.809, 321.479, 327.045, 334.191, 335.418, 341.136,$
	341.987, 345.148, 345.753, 346.867, 348.199, 349.895, 352.131,$
	352.699, 358.694, 359.852, 360.782, 364.494, 369.205, 417.296]
if serts eq 'qut93' then xdxt_w=$
	[274.157, 284.142, 292.852, 296.110, 311.574, 311.783, 312.171,$
	313.732, 314.353, 315.022, 316.219, 317.026, 319.843, 320.802,$
	321.464, 327.041, 334.188, 335.417, 339.017, 341.141, 341.982,$
	344.974, 345.748, 346.871, 348.196, 349.885, 352.127, 352.694,$
	358.662, 359.851, 360.785, 364.490, 365.570, 369.158, 417.294]
intersect,wvl,xdxt_w,y1y2,y1n2,n1y2

;	these are Jeremy's quality grades for each line (1=good)
grade=[4,1,2,2,2,1,1,1.5,1,2,2,1.5,2,2,1,1,1,1,2,1.5,1.5,1,3,1.5,2,3,$
       1,1,1,1,1,1,1,1,1,1,1,1,3,1,1.5,1,2.5,1,1,1,1,3,2,3,2,3,1,1,3,$
       3,1,1.5,2.5,1,2]

;	FLX is in units of [ergs/s/cm^2/sr] -- convert to [ph/s/cm^2]
nrg=6.626176e-27*2.9979e10*1e8/abs(wvl)
flx=4.*!pi*flx/nrg & sigma=4.*!pi*sigma/nrg

;	restore ID info
restore,'id_'+serts+'.sav'		;restores ALLID
;WARNING: not all of WVL may have a match!
x=(allid.wvl) & nx=n_elements(x)	;these wavelengths have matches
y=fltarr(nx) & sy=y			;corresponding fluxes & errors

;	pick subset of wavelengths
ow=lindgen(nw)
;ow=where(grade le 1.5)
;ow=y1y2
wvl=wvl(ow) & flx=flx(ow) & sigma=sigma(ow)

;	extract ID info
mwvl=[0.] & midx=[-1L] & Z=[0] & ion=Z
for ix=0,nx-1 do begin
  oo=where(x(ix) eq wvl,moo)
  ;if moo eq 0 then message,'ID for non-existent wavelength?!'
  if moo gt 0 then begin
    y(ix)=flx(oo(0)) & sy(ix)=sigma(oo(0))
    mx=allid.(ix+1).wvl & mz=allid.(ix+1).Z & mi=allid.(ix+1).ion
    m_id=lonarr(n_elements(mx))+ix
    mwvl=[mwvl,mx] & midx=[midx,m_id] & Z=[Z,mz] & ion=[ion,mi]
  endif
endfor
mwvl=mwvl(1:*) & midx=midx(1:*) & Z=Z(1:*) & ion=ion(1:*)
nwvl=n_elements(mwvl)

;	figure out how many unique ions
atm=100*Z+ion & uatm=atm(uniq(atm,sort(atm))) & nuatm=n_elements(uatm)
em=dblarr(2,nuatm) & tmx=dblarr(nuatm)

;	set up plot
plot,[5,7],[1e26,6e28],/nodata,/xs,/ys,/yl,xtitle='log!d10!n(T [K])',$
	ytitle='EM [cm!u-5!n]',title='Ion Averaged EMs'

;************************************************************
;	call POTAVE
;************************************************************
for i=0,nuatm-1 do begin			;{for each unique ion
  zz=long(uatm(i)/100) & ii=long(uatm(i)-100*zz)
  elem=atom(zz-1)+' '+rom(ii-1)
  print,'working on '+elem
  ;
  oo=where(Z eq zz and ion eq ii,moo)		;matches for all such ELEM
  if moo eq 0 then message,'bug!'
  m_wvl=mwvl(oo) & m_idx=midx(oo)		;select 'em
  ;
  oz=m_idx(uniq(m_idx,sort(m_idx)))		;corresponding lines
  f=y(oz) & w=x(oz) & s=sy(oz)			;select 'em
  ;
  mindex=0*m_idx				;match from IDs to lines
  for iz=0,n_elements(oz)-1 do mindex(where(m_idx eq oz(iz)))=iz
  ;
  ;************************************************************
  tmp=potave(elem,f,w,m_wvl,mindex,sigma=s,ewt=ewt,tmax=tmax,$
	abund=abund,logP=logP,n_e=n_e,dbdir=dbdir,chifil=chifil)
  ;************************************************************
  print,tmp,tmax
  em(*,i)=tmp & tmx(i)=total(tmax)/n_elements(tmax)
  oplot,[tmax],[tmp(0)],psym=1
  oplot,tmax*[1,1],tmp(0)+tmp(1)*[1,-1]
  xyouts,tmax,tmp(0),elem,align=-0.5,col=150
endfor						;I=0,NUATM-1}

end
