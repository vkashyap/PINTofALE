;+ 
;
;script rhessi_synth 
;    compute synthetic rhessi spectrum 
; 
;01/04/05
;-

;--------------------------------------------------------------------------
;	control variables
;--------------------------------------------------------------------------
;	environment
 !LDBDIR = '$CHIANTI'		;'$CHIANTI','$SPEX','$APED','/path/to/dir'
 pimmsdir='/soft/pimms/data'    ; path to local PIMMS installation

;	source parameters
 !NH           = 3e20			; H column density [cm^-2]
 EDENS1        = 1e12   		; electron number density [cm^-3] 
 !ABUND        = getabund('feldman')   ; elemental abundances (SEE:getabund.pro)
 T_components  = [7.3,7.4,7.8]		; log(T[K]) components in the EM
 EM_components = [6e25,6e25,4e34]	; Emission Measure [cm^-3]
 obs_rate    = 1.0                    ; ROSAT/PSPC count rate [cts/s] -- SET TO ZERO IF NOT REQUIRED
 sr2oau2     = (!fundae.rsun^2)/(!fundae.au^2) ; 
 restore, '/home/llin/solar/mcmc/flare.v5/v16dem.save'

;	observation parameters
 EXPTIME       = 0.001			 ; nominal exposure time [ks]
  
;file_GRMF     = '/home/llin/solar/dFLARE/hsi_srm_20030426_030000_det4.fits' ; full path name to grating RSP
;file_GRMF     = '/home/llin/solar/dFLARE/hsi_srm_20030426_022500_det4.fits'
file_GRMF     = '/home/llin/solar/dFLARE/rhessi_OGIP/hsi_qsrm_20030426_024500.fits'
; restore, 'v.dem.save'

;       analysis parameters
 !WMIN         = 0.12398521             ; minimum wavelength for output spectrum [ang]
 !WMAX         = 4.1328402              ; maximum wavelength for output spectrum [ang] 
 nwbin         = 600L			; number of bins in theoretical spectrum [ang] 

;restore, 'v10.save',/v

;	make the dem 
;oo = where(logt lt 6.5) 
;bestdem(oo) = bestdem(oo)-bestdem(oo)+1e7*max(bestdem)
!DEM=mk_dem('interpolate',logT=!LOGT,pardem=logt,indem=bestdem)
;!DEM=mk_dem('delta', logT = !LOGT, pardem=T_components, indem=EM_components)
 ;A] Read in line cooling emissivities and calculate line intensities
	
               ;Read  line cooling emissivities of all possible
	       ;lines in the ROSAT/PSPC wavelength range from the atomic data base. 
	       !EDENS = EDENS1
	       lconf=rd_line(atom,n_e=!EDENS,wrange=[!wmin,!wmax],$
			     dbdir=!LDBDIR,verbose=!VERBOSE,wvl=LWVL,logT=LLOGT,Z=Z,$
			          ion=ION,jon=JON,fstr=lstr)

	       ;The output of rd_line.pro will only include level population,
               ;and not ion balances. We will use fold_ioneq.pro to fold ion balances.
               if strpos(strlowcase(!LDBDIR),'aped',0) lt 0 then lconf=$
               fold_ioneq(lconf,Z,JON,chidir=!CHIDIR,$
                                logT=LLOGT,eqfile=!IONEQF,verbose=!VERBOSE)      

               if strpos(strlowcase(!LDBDIR),'aped',0) ge 0 then  v_ABUND =$
               !ABUND/getabund('anders & grevesse') else v_ABUND=!ABUND

               ;And now calculate line intensities using lineflx.pro.

               linint=sr2oau2*lineflx(lconf,!LOGT,LWVL,Z,DEM=!DEM,abund=v_ABUND) ;[ph/s]

 ;B] Read in continuum emissivities and calculate continuum intensities
 	  
	       ;We can read in continuum emissivities using rd_cont.pro.
	       ;It is important to note that the output emissivities of rd_cont.pro
	       ;are in [1e-23 erg cm^3/s/Ang] and not [1e-23 erg cm^3/s] as with rd_line.pro
	       !EDENS=EDENS1
               cconf=rd_cont(!CEROOT,n_e=!EDENS,wrange=[!wmin, !wmax],$
                      dbdir=!CDBDIR,abund=!ABUND,verbose=!VERBOSE,$
                      wvl=CWW,logT=ClogT)

	       ;The continuum intensities per angstrom can be calculated again using
	       ;lineflx.pro. Note that CWW contains the wavelength bin boundaries for 
	       ;the emissivity array.
	       CWVL=0.5*(CWW[1:*]+CWW )
	       conint=sr2oau2*lineflx(cconf,!LOGT,CWVL,DEM=!DEM)   ;[ph/s/Ang]

	       ;Now to get just continuum intensity, we must multiply by an array
	       ;containing the bin widths. If we define this array simply
	       ;with:  CDW=CWW[1:*]-CWW, we will get an ugly 'saw-toothed' figure. 
	       ;(a side-effect of the way the data-base is constructed) To work 
	       ;around this, we can use CWVL, the mid-bin values, and mid2bound.pro, 
	       ;which gives intelligent bin-boundary values given mid-bin values:

	       CWB=mid2bound(CWVL) & CDW=CWB[1:*]-CWB
	       conint=conint*CDW       ;[ph/s/Ang]*[Ang]


      ;D] Bin spectra and fold in effective area 

	       ;Rebin to form theoretical line spectrum using hastrogram.pro
	       linspc = hastogram(abs(LWVL),wgrid,wts=linint)  ;[ph/s/cm^2/bin]

	       ;Rebin to form theoretical continuum spectrum using rebinw.pro
	       conspc = rebinw(conint,CWVL,wgrid,/perbin)      ;[ph/s/cm^2/bin]

	       ;Derive predicted flux spectrum.
	       WVLS=0.5*(WGRID[1:*]+WGRID)
;	       newEffAr=(interpol(effar,pspc_wvlar,WVLS) > 0) < (max(pspc_effar))
	       flxspc = (linspc + conspc)

	       ;Derive predicted counts spectrum.
	       flxspc=flxspc*EXPTIME*1e3*39.5              ;[ct/bin]
               EGRID=!fundae.KEVANG/WGRID

	 ;E] Convolve with RMF 

		conv_rmf,EGRID,flxspc,CHAN,CTSPC,file_gRMF,rmfstr=grmf,effar=effar
		
restore, '/home/llin/solar/dFLARE/rhessi_SAVE/hsi_qspec_20030426_024500.resik.save',/v 
hsi_rmf = rd_ogip_rmf('/home/llin/solar/dFLARE/rhessi_OGIP/hsi_qsrm_20030426_024500.fits',effar=effar) 
enrgrid = rhstr.egrid & egrid1 = (hsi_rmf.elo+hsi_rmf.ehi)/2
wvlar = !fundae.kevang/egrid1
o1 = where(egrid ge 3.0) 
enrgrid = rhstr.egrid(o1) 
spec  = rhstr.spec(o1) 
specerr = rhstr.specerr(o1) 

nspec = rebinw(spec,enrgrid,chan) 
plot, chan,ctspc*0.83, /yl, xr = [0,10] 
oplot, chan, nspec, linestyle = 1


        ;F] The final step is a simulation of counts based on the above spectrum
	
;	nbin=n_elements(ctspc)
;	sim_spec=intarr(nbin)
;	CTSIM=intarr(nbin) & CTSIM0=intarr(nbin)
;        for i=0L,nbin-1L do if ctspc[i] gt 0 then $
;	     CTSIM[i]=randomu(seed,poisson=double(ctspc[i]))
        
;end
