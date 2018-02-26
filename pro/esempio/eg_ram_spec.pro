;+
;eg_ram_spec.pro
;	make a spectrum in the 2-200 AA region (for RAM, which is
;	supposed to have a calorimeter as a detector) with a resolution
;	of 2 eV
;
;usage
;	ldbdir='$CHIANTI' & verbose=5 & abref='grevesse et al'
;	edens=1e9 & NH=1e10 & betap=2.5
;	outfil='/tmp/ram_spec.save'
;	.run eg_ram_spec
;
;vinay kashyap (28may02; Ed Deluca wanted it)
;-

;	initialize variables
if not keyword_set(verbose) then verbose=5
if not keyword_set(abref) then abref='grevesse et al'
if not keyword_set(edens) then edens=1e9
if not keyword_set(NH) then NH=1e10
if not keyword_set(ldbdir) then ldbdir='$CHIANTI'
if not keyword_set(cdbdir) then cdbdir='$CONT'
if not keyword_set(ceroot) then ceroot='cie'
if not keyword_set(ioneqf) then ioneqf='ioneq/mazzotta_etal.ioneq'
if n_elements(abund) ne 30 then abund=getabund(abref)
if n_tags(fundae) eq 0 then inicon,fundae=fundae
if not keyword_set(outfil) then outfil='/tmp/ram_spec.save'

;	define the DEM
brosdem=[22.0,20.7,20.4,20.5,21.1,21.5,21.1,21.3,19.7]
broslgT=[5.0, 5.4, 5.6, 5.7, 6.0, 6.2, 6.4, 6.7, 7.0]
logT=findgen(81)*0.05+4
DEM=mk_dem('spline',logT=logT,indem=brosdem,pardem=broslgT)
oo=where(logT lt 5 or logT gt 7,moo) & if moo gt 0 then DEM[oo]=0.
oo=where(logT ge 5 and logT le 7,moo) & if moo gt 0 then DEM[oo]=10.D ^ (DEM[oo])

;	define the wavelength grid
wmin=2. & wmax=200. & emin=fundae.keVAng/wmax & emax=fundae.keVAng/wmin
dE=2e-3/20.	;2 eV resolution for the calorimeter, oversample by 20
if not keyword_set(betap) then betap=2.5	;smoothing function is a beta-profile
nbinE=long((emax-emin)/dE)+1L & egrid=findgen(nbinE)*dE+emin
wgrid=fundae.keVAng/egrid & os=sort(wgrid) & wgrid=wgrid[os]

;	get the line contribution function
lconf=rd_line(atom,n_e=edens[0],wrange=[wmin,wmax],dbdir=ldbdir[0],$
	verbose=verbose[0],wvl=lwvl,logT=llogT,Z=lZ,ion=lion,jon=ljon,$
	fstr=lstr,/desig,/econf)
lconf=fold_ioneq(lconf,lZ,ljon,chidir=chidir,logT=llogT,eqfile=ioneqf[0],$
	verbose=verbose[0])

;	get the continuum contribution function
cconf=rd_cont(ceroot[0],n_e=edens[0],wrange=[wmin,wmax],dbdir=cdbdir[0],$
	abund=abund,verbose=verbose[0],wvl=cww,logT=clogT,fcstr=cstr)
cwvl=0.5*(cww[1:*]+cww) & cdw=cww[1:*]-cww

;	derive intensities
linint=lineflx(lconf,logT,lwvl,lZ,DEM=DEM,abund=abund)	;[ph/s/cm^2]
conint=lineflx(cconf,logT,cwvl,DEM=DEM)*cdw		;[ph/s/cm^2/Ang]*[Ang]

;	ISM absorption
ltau=ismtau(abs(lwvl),NH=NH,fH2=fH2,He1=He1,HeII=HeII,$
	Fano=Fano,wam=wam,/bam,abund=abund)
ctau=ismtau(abs(cwvl),NH=NH,fH2=fH2,He1=He1,HeII=HeII,$
	Fano=Fano,wam=wam,/bam,abund=abund)
ltrans=exp(-ltau) & ctrans=exp(-ctau)

;	fluxes
linflx=linint*ltrans
conflx=conint*ctrans

;	rebin to form spectrum
linspc0=hastogram(abs(lwvl),wgrid,wts=linflx)	;[ph/s/cm^2/bin]
conspc0=rebinw(conflx,cwvl,wgrid,/perbin)	;[ph/s/cm^2/bin]

;	smooth with LRF
hwhm=10. & smscl=hwhm/sqrt(2.^(1./betap)-1.)
linspc=voorsmooth(linspc0,smscl,type='beta',betap=betap)
conspc=voorsmooth(conspc0,smscl,type='beta',betap=betap)
totspc=linspc+conspc

;	plot
egrid=fundae.keVAng/wgrid
plot,egrid,totspc/dE,/yl,/xl,/xs,thick=3,$
	xtitle='Energy [keV]',ytitle='flux [ph/s/cm^2/A]'
oplot,egrid,linspc/dE,color=1 & oplot,egrid,conspc/dE,color=2
stample & setkolor,'red',1 & setkolor,'green',2
help=[	'plot,egrid,totspc/dE,/yl,/xl,/xs,thick=3,$',$
	"xtitle='Energy [keV]',ytitle='flux [ph/s/cm^2/A]",$
	'oplot,egrid,linspc/dE,color=1 & oplot,egrid,conspc/dE,color=2',$
	"stample & setkolor,'red',1 & setkolor,'green',2" ]

;	save variables
save,file=outfil,help,wgrid,totspc,linspc,conspc,dE,egrid,betap,hwhm,smscl,$
	linspc0,conspc0,verbose,abref,edens,NH,ldbdir,cdbdir,ceroot,ioneqf,abund

end
