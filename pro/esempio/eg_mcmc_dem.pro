;+
;EG_MCMC_DEM
;	create an arbitrary set of line fluxes from an arbitrary DEM
;	and see if we can get back the DEM
;-

;	read in test lines and IDs from (say) SERTS
serts='act93' & filnam=serts+'.sav'
restore,filnam
x=allid.wvl & nx=n_elements(x) & allz=intarr(nx)
for i=0,nx-1 do allz(i)=(allid.(i+1).z)(0)

;	make example spectrum and tweak DEM
avdem=total(dem)/n_elements(dem) & dem(*)=avdem/1e14
loadct,3
chifil='ioneq/arr.ioneq'	;ion balance including new Fe data
fidgit,lamda,spec/1e14,elem,dem=dem,logt=logt,sigma=sig,ww=ww,$
  logP=logP,dbdir=dbdir,abund=abund,chifil=chifil
wset,0

;	get new line fluxes
fl=lineflx(line,logT,x,allz,DEM=DEM,abund=abnd)
ulim=long(0*fl)				;upper limits
sigma=sig/1e14 > (1.+sqrt(fl+0.75))
oo=where(fl le sigma,moo) & if moo gt 0 then ulim(oo)=1
print,'there are ',moo,' upper limits'
if moo gt 0 then flx(oo)=sig(oo)

;	get a subset of the wavelengths to consider...
;e.g., ow=lindgen(nx), ow=where(allz eq 26), ow=where(x lt 302 and x gt 305),
;ow=where(allz ne 2), ow=where(allz ne 2 and ulim eq 0), etc.
ow=where(allz ne 2 and ulim eq 0)
mow=n_elements(ow)
wvl=x(ow) & flx=fl(ow) ;& sigma=sig(ow)/1e14 > (1.+sqrt(flx+0.75))
zz=allz(ow) & emis=line(*,ow)

;	get a subset of the temperatures to consider...
;e.g., ot=lindgen(nt), ot=where(logt ge 5 and logt le 7), ot=lindgen(41)+20,
;ot=lindgen(21)*2+20
ot=where(logt ge 5 and logt le 7)
mot=n_elements(ot) & tlog=logt(ot) & diffem=dem(ot) & emis=emis(ot,*)

;	keywords
nbatch=10               ;number of simulations/batch
nsim=100		;number of batches
ngroup=0		;group adjoining DEM bins
nosrch=0		;don't carry out initial grid search (1) or do (0)?
savfil='test_mcmc.sav'	;save results in this file
sphere=0		;change parameters all at once (1) or one by one (0)?
bound=0.68		;fraction of simulations that should lie between upper
			;(or lower) bound and MAP estimate
hwhm=0			;return upper and lower HWHMs instead of bounds (1)?
rchi=0			;reduced chi^2 (1) or plain ol' chi^2 (0)?

;	DEM ranges
mindem=min(diffem,max=maxdem)
demrng=[0*diffem+mindem/1e3,0*diffem+maxdem*1e3] & demrng=reform(demrng,mot,2)

;	ABUND ranges
minab=min(abund,max=maxab) & nab=n_elements(abund)
abrng=[abund,abund] & abrng=reform(abrng,nab,2)

;	call MCMC_DEM
diffem(*)=total(diffem)/mot
;dem_2k=mcmc_dem1(wvl,flx,emis,z=zz,logT=tlog,diffem=diffem,abund=abund,$
;	sigma=sigma,demrng=demrng,abrng=abrng,nsim=nsim,nbatch=nbatch,$
;	nosrch=nosrch,savfil=savfil,bound=bound,hwhm=hwhm,sphere=sphere,$
;	demerr=demerr,aberr=aberr,rchi=rchi)
dem_2k=mcmc_dem(wvl,flx,emis,z=zz,logT=tlog,diffem=diffem,abund=abund,$
	sigma=sigma,demrng=demrng,abrng=abrng,nsim=nsim,nbatch=nbatch,$
	ngroup=ngroup,nosrch=nosrch,savfil=savfil,bound=bound,hwhm=hwhm,$
	demerr=demerr,aberr=aberr,rchi=rchi,ulim=ulim)

;	plot results
!p.multi=0
oo=where(demerr(*,0) ne demerr(*,1))
plot,tlog(oo),dem_2k(oo),yr=[min(demerr),max(demerr)],/yl,$
	xr=[min(tlog),max(tlog)],col=150
for i=0,mot-1 do oplot,tlog(i)*[1,1],demerr(i,*),col=150
oplot,logt,dem,line=1

end
