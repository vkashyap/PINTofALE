;+
;POA_DEMO
;	set up data and instrument arrays for use in PoA demos
;-

;	first, Chandra grating data of HR1099

;	read in spectra
restore,'/fubar/kashyap/CXC/HETG/acisf62538_000N001_tg2.sav'
oo=sort(bin_lo09)
spc_m1p_v711Tau=counts09[oo]
spe_m1p_v711Tau=1.+sqrt(abs(spc_m1p_v711Tau)+0.75)
ww_m1p_v711Tau=[bin_lo09[oo],max(bin_lo09)]
lam_m1p_v711Tau=0.5*(ww_m1p_v711Tau[1:*]+ww_m1p_v711Tau)

;	read in ARFs
restore,'/fubar/kashyap/CXC/HETG/HR1099/Demo/garf.sav'
effar_m1p=arf_megp1
wvlar_m1p=wvl_megp1
exptim_v711Tau=exptim

;	read in DEM
restore,'/data/drake7/vinay/SCAR/emissivity/hr1099dem.sav',/verbose
DEM_v711Tau=mk_dem('interpolate',logT=logT_v711Tau,$
	indem=DEM(oo),pardem=logT(oo))

;	plot data
plot,lam_m1p_v711Tau,spc_m1p_v711Tau,/nodata,xr=[5.,25.],/xl,/xs,$
	yr=[0.1,1200.],/ys,xtitle='Wavelength ['+string(byte(197))+']',$
	ytitle='[counts/bin]',title='HR1099: HETG/MEG +1',color=1
oplot,lam_m1p_v711Tau,spc_m1p_v711Tau,psym=10,color=2
for i=0L,n_elements(spe_m1p_v711Tau)-1L do oplot,lam_m1p_v711Tau[i]*[1,1],$
	spc_m1p_v711Tau[i]+spe_m1p_v711Tau[i]*[-1,1],color=3
stample,stacol=4
setcolor,'DarkGreen',1
setcolor,'blue',2
setcolor,'red',3
setcolor,'yellow',4

;	save to file
save,file='/fubar/SCAR/ardb/demo_v711Tau.save',spc_m1p_v711Tau,$
	spe_m1p_v711Tau,ww_m1p_v711Tau,lam_m1p_v711Tau,$
	effar_m1p,wvlar_m1p,exptim_v711Tau,logT_v711Tau,DEM_v711Tau

end
