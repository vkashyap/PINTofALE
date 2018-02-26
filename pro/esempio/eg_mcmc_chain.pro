;+
;EG_MCMC_DEM
;
;example program to call MCMC_CHAIN
;
;vinay kashyap (Aug2008)
;-

;	set up the color table
peasecolr & loadct,3 & peasecolr

;	filenames of PHA spectrum, ARF, and RMF
phafil='/data/snafu/acisshetg/1014/tg_reprocess/acisf01014N002_MEG_1_NONE.pha'
arfil='/data/snafu/acisshetg/1014/tg_reprocess/acisf01014N002MEG_1_garf.fits'
rmfil='/data/snafu/acisshetg/1014/tg_reprocess/MEG1.rmf'

;	read in the data (and background) and set up the x-axis in Energy [keV]
if n_tags(t) eq 0 then t=mrdfits(phafil,1,h)
x=12.3985/(0.5*(t.BIN_LO+t.BIN_HI)) & y=t.COUNTS & bgy=(t.BACKGROUND_UP+t.BACKGROUND_DOWN)  
bkgscal=sxpar(h,'BACKSCUP')+sxpar(h,'BACKSCDN')
if n_tags(arstr) eq 0 then arf=rdarf(arfil,arstr)
if n_tags(rmf) eq 0 then rmf=rd_ogip_rmf(rmfil,effar=effar)

;	set up the energy filter
xfilter=where(x ge 0.4 and x le 7.5)

;	define the background model
btype=['power']				;model type = broken power law
bparval=[1e2,-3,-3,10]			;NORM,GAMMA1,GAMMA2,BREAK
bparsig=[-0.05,-0.05,0.1,0.1]		;guess at uncertainties
					;(-ve == fractional error, +ve == sigma)
bfreeze=[2,3]				;freeze GAMMA2 and BREAK
bties=['a0=(a0)>(1e-5)','a1=a1<0']	;parameter constraints

type=['absorb','power']			;ISMabs * broken-power-law
parval=[1e20,3e3,-2.5,-2.5,2]		;NH, NORM,GAMMA1,GAMMA2,BREAK
parsig=[-0.1,-0.05,0.1,0.1,0.1]		;guess at uncertainties
ties=['a0=a0>1e19','a1=a1>0.001','a2=a2<0','a3=a3<0','a4=((a4>1)<15)']
					;parameter constraints

;	call MCMC_CHAIN
parstr=mcmc_chain(x,y,parval,parsig,verbose=11,$
	type=type,/ikev,arfil=arfil,rmfil=rmfil,$
	ties=ties,xfilter=xfilter,/onepar,/cash,$
	bgy=bgy,bkgscal=bkgscal,bparval=bparval,bparsig=bparsig,$
	btype=btype,bfreeze=bfreeze,bties=bties)

;	all outputs are in
help,parstr,/str

end
