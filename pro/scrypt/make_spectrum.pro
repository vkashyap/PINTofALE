;+
;script 	make_spectrum
;	compute a spectrum using line and continuum emissivities
;
;warning
;	all program variables have names that begin as "v_" or "jnk_"
;	if you have any preexisiting variables with these names, they
;	will be overwritten.
;
;inputs
;	ATOMIC:
;	  line emissivities and wavelengths
;	  continuum emissivities and wavelength grid
;	  (subsumed in database directory names)
;	  ion fractions
;	STELLAR:
;	  coronal electron density
;	  coronal abundances
;	  Differential Emission Measure
;	  distance
;	  Galactic coordinates (or absorption column density)
;	OBSERVATIONAL:
;	  instrument effective area
;	  synthetic spectrum wavelength grid
;	  photon redistribution matrix
;	  exposure time
;
;outputs
;	v_linint	line intensities @ source
;		(@ wavelengths v_lwvl; from species v_Z,v_ion)
;	v_conint	continuum intensities @ source
;		(@ wavelength grid v_cww)
;	v_linflx	theoretical line fluxes @ earth
;		(using ISM transmission factors v_ltrans)
;	v_conflx	theoretical continuum fluxes @ earth
;		(using ISM transmission factors v_ctrans)
;	v_linspc	theoretical line spectrum @ earth
;	v_conspc	theoretical continuum spectrum @ earth
;	v_flxspc	predicted spectrum through telescope
;	v_wvls  	mid-bin values of predicted spectrum wavelength grid
;	v_ctspc 	predicted count spectrum for detector
;	v_chan  	energy channels for count spectrum
;
;how to use
;	1. initialize the PINTofALE environment using INITALE
;	2. the parameters that will be set with INITALE are:
;	   --	ATOMIC: !LDBDIR, !IONEQF, !CDBDIR, !CEROOT, !CHIDIR
;	   --	STELLAR: !EDENS, !ABUND, !NH, !FH2, !HE1, !HEII
;	   verify that they are initialized to the right values
;	3. define DEM as a function of LOGT (use DEMACS() or MK_DEM())
;	   --	place in variables v_LOGT [K] and v_DEM [cm^-3]
;	4. define distance in [pc]
;	   --	place in variable v_DIST
;	5. read in the instrument effective area
;	   --	place in variables v_EFFAR [cm^2] and v_WVLAR [Ang]
;	6. read in (or define) the synthetic spectrum wavelength grid
;	   --	place in variable v_WGRID [Ang]
;	   (must be _bin boundaries_)
;	7. define the exposure time [ks]
;	   --	place in variable v_EXPTIME
;	8. define the name of the OGIP-style RMF
;	   --	place in variable v_RMF
;	   (set to 'NONE' if, e.g., grating spectrum is being analyzed)
;	9. run this script
;
;history
;	vinay kashyap (OctMM)
;-

;	{	check that all variables are defined

;2:
v='LDBDIR' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_LDBDIR,/getval
v='IONEQF' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_IONEQF,/getval
v='CDBDIR' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_CDBDIR,/getval
v='CEROOT' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_CEROOT,/getval
v='CHIDIR' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_CHIDIR,/getval
v='EDENS' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_EDENS,/getval
v='ABUND' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_ABUND,/getval
v='NH' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_NH,/getval
v='FH2' & defsysv,'!'+v,exists=i
if i eq 0 then setsysval,'FH2',0.26,verbose=2
v='HE1' & defsysv,'!'+v,exists=i
if i eq 0 then setsysval,'HE1',0.,verbose=2
v='HEII' & defsysv,'!'+v,exists=i
if i eq 0 then setsysval,'HEII',0.,verbose=2
v='VERBOSE' & defsysv,'!'+v,exists=i
if i eq 0 then setsysval,'VERBOSE',5,verbose=2

;	expected units
message,'v_DEM must be in [cm^-3]',/info
message,'v_DIST must be in [pc]',/info
message,'v_EFFAR must be in [cm^2]',/info
message,'v_WVLAR must be in [Ang]',/info
message,'v_WGRID must be in [Ang] and must be sorted in ascending order',/info
message,'v_EXPTIME must be in [ksec]',/info
message,'v_RMF must be an OGIP-compatible response matrix file name',/info

;3:
nDEM=n_elements(v_DEM) & nlogT=n_elements(v_LOGT)
if nDEM eq 0 then message,'v_DEM: DEM is undefined'
if nlogT eq 0 then message,'v_LOGT: temperature grid is undefined'
if nDEM ne nlogT then message,'v_DEM(v_LOGT) and v_LOGT are incompatible'

;4:
if not keyword_set(v_DIST) then begin
  message,'v_DIST: Distance not defined.  Setting to 1 [pc]',/info
  v_DIST=1.
endif

;5:
nea=n_elements(v_EFFAR) & nwa=n_elements(v_WVLAR)
if nea eq 0 then message,'v_EFFAR: Effective area is undefined'
if nwa eq 0 then message,'v_WVLAR: wavelength scale is undefined'
if nea ne nwa then message,'v_EFFAR(v_WVLAR) and v_WVLAR are incompatible'

;6:
nww=n_elements(v_WGRID)
if nww eq 0 then message,'v_WGRID: synthetic spectrum grid is undefined'
if nww eq 1 then message,'v_WGRID: must be defined to be bin _boundaries_'

;7:
if not keyword_set(v_EXPTIME) then begin
  message,'v_EXPTIME: exposure time not defined.  Setting to 1 [ks]',/info
  v_EXPTIME=1.
endif

;8:
nr=n_elements(v_RMF) & szr=size(v_RMF) & nszr=n_elements(szr)
if nr eq 0 then message,'v_RMF: RMF undefined'
if szr[nszr-2] ne 7 then message,'v_RMF: must be filename'

;	OK, all set!	}

;	now define some useful variables
;	(these are used to figure out if databases must be reread)
v_WMIN=min(v_WGRID,max=v_WMAX)
if not keyword_set(jnk_LDBDIR) then jnk_LDBDIR=v_LDBDIR+' old'
if not keyword_set(jnk_CDBDIR) then jnk_CDBDIR=v_CDBDIR+' old'
if not keyword_set(jnk_CEROOT) then jnk_CEROOT=v_CEROOT+' old'
if not keyword_set(jnk_IONEQF) then jnk_IONEQF=v_IONEQF+' old'
if not keyword_set(jnk_EDENS) then jnk_EDENS=v_EDENS+1e10
if not keyword_set(jnk_ABUND) then jnk_ABUND=v_ABUND*0.9
if not keyword_set(jnk_WMIN) then jnk_WMIN=v_WMIN+100
if not keyword_set(jnk_WMAX) then jnk_WMAX=v_WMAX+100

;	and now follow the MAKE_SPECTRUM thread
;	(see $PoA/doc/make_spectrum.ps)

;	obtain line contribution function
;	[1e-23 ergs cm^3/s]
if (v_LDBDIR ne jnk_LDBDIR) or (v_IONEQF ne jnk_IONEQF) or $
   (v_EDENS ne jnk_EDENS) or (v_WMIN ne jnk_WMIN) or $
   (v_WMAX ne jnk_WMAX) then begin
  ;	read in line database
  v_lconf=rd_line(atom,n_e=v_EDENS[0],wrange=[v_WMIN,v_WMAX],$
	dbdir=v_LDBDIR[0],verbose=v_VERBOSE,$
	wvl=v_LWVL,logT=v_LLOGT,Z=v_Z,ion=v_ION,jon=v_JON,fstr=v_lstr)
  ;	include ion balance
  if strpos(strlowcase(v_IONEQF),'none',0) lt 0 then begin
    v_lconf=fold_ioneq(v_lconf,v_Z,v_JON,chidir=v_CHIDIR,$
	logT=v_LLOGT,eqfile=v_IONEQF,verbose=v_VERBOSE)
    v_lstr.LINE_INT = v_lconf
  endif
endif
if n_elements(v_LLOGT) ne n_elements(!LOGT) then message,$
	'!LOGT incompatible with input grid v_LLOGT; cannot continue.'

;	obtain continuum contribution function
;	[1e-23 ergs cm^3/s/Ang]
if (v_CDBDIR ne jnk_CDBDIR) or (v_CEROOT ne jnk_CEROOT) or $
   (v_EDENS ne jnk_EDENS) or (v_WMIN ne jnk_WMIN) or $
   (v_WMAX ne jnk_WMAX) or (not arrayeq(v_ABUND,jnk_ABUND)) then begin
  ;	read in continuum database
  v_cconf=rd_cont(v_CEROOT[0],n_e=v_EDENS[0],wrange=[v_WMIN,v_WMAX],$
	dbdir=v_CDBDIR[0],abund=v_ABUND,verbose=v_VERBOSE,$
	wvl=v_CWW,logT=v_ClogT)
endif
if n_elements(v_CLOGT) ne n_elements(!LOGT) then message,$
	'!LOGT incompatible with input grid v_CLOGT; cannot continue.'

;	derive line intensity
;	place DEM into proper grid
if nDEM ne n_elements(v_LLOGT) then begin
  message,'Regridding the DEM',/info
  v_DEM=mk_dem('interpolate',logT=!LOGT,indem=v_DEM,pardem=v_LOGT)
  ;oo=where(v_LLOGT lt 5 or v_LLOGT gt 7) & v_DEM[oo]=0.
endif
v_linint=lineflx(v_lconf,!LOGT,v_LWVL,v_Z,DEM=v_DEM,abund=v_ABUND)
	;[ph/s]

;	derive continuum intensity
;	place DEM into proper grid
if nDEM ne n_elements(v_CLOGT) then begin
  message,'Regridding the DEM',/info
  v_DEM=mk_dem('interpolate',logT=!LOGT,indem=v_DEM,pardem=v_LOGT)
  ;oo=where(v_CLOGT lt 5 or v_LLOGT gt 7) & v_DEM[oo]=0.
endif
v_CWVL=0.5*(v_CWW[1:*]+v_CWW) & v_CDW=v_CWW[1:*]-v_CWW
;	hack to avoid annoying saw-toothed CDW
v_CWB=mid2bound(v_CWVL) & v_CDW=v_CWB[1:*]-v_CWB
v_conint=lineflx(v_cconf,!LOGT,v_CWVL,DEM=v_DEM)*v_CDW
	;[ph/s/Ang]*[Ang]

;	distance correction factor
v_LogDCOR = alog10(4.*!DPI)+2.*alog10(v_DIST*3.1e18)

;	derive ISM absorptions
v_ltau=ismtau(v_LWVL,NH=v_NH,fH2=v_fH2,He1=v_He1,HeII=v_HeII,$
	Fano=Fano,wam=wam,/bam,abund=v_ABUND,verbose=v_VERBOSE)
v_ctau=ismtau(v_CWW,NH=v_NH,fH2=v_fH2,He1=v_He1,HeII=v_HeII,$
	Fano=Fano,wam=wam,/bam,abund=v_ABUND,verbose=v_VERBOSE)
v_ltrans=exp(-v_ltau) & v_ctrans=exp(-v_ctau)

;	derive theoretical line fluxes
;	[ph/s/cm^2]
v_linflx = v_linint * v_ltrans / 10.^(v_LogDCOR)

;	derive theoretical continuum fluxes
;	[ph/s/cm^2]
v_conflx = v_conint * v_ctrans / 10.^(v_LogDCOR)

;	rebin to form theoretical line spectrum
;	[ph/s/cm^2/bin]
v_linspc = hastogram(abs(v_LWVL),v_WGRID,wts=v_linflx)

;	rebin to form theoretical continuum spectrum
;	[ph/s/cm^2/bin]
v_conspc = rebinw(v_conflx,v_CWVL,v_WGRID,/perbin)

;	derive predicted flux spectrum
v_WVLS=0.5*(v_WGRID[1:*]+v_WGRID)
v_newEffAr=(interpol(v_EFFAR,v_WVLAR,v_WVLS) > 0) < (max(v_EFFAR))
;	[ct/s/bin] (if DEM is [cm^-5]: [ct/s/cm^2/bin])
v_flxspc = (v_linspc + v_conspc) * v_newEffAr

;	plot
plot,v_WVLS,v_flxspc,psym=10,xtitle='[Ang]',ytitle='[ct/s/bin]',$
	title='flux spectrum'
oplot,v_WVLS,v_linspc*v_newEffAr,psym=10,col=200
oplot,v_WVLS,v_conspc*v_newEffAr,psym=10,col=150
stample

;	derive predicted counts spectrum
v_flxspc=v_flxspc*v_EXPTIME*1e3		;[ct/bin]
if strlowcase(v_RMF) ne 'none' then begin
  v_NRG=12.3985/v_WGRID
  conv_rmf,v_NRG,v_flxspc,v_CHAN,v_CTSPC,v_RMF,rmfstr=v_rstr,shift1=shift1
  ;	plot
  plot,v_CHAN,v_CTSPC,psym=10,xtitle='[keV]',ytitle='[ct]',$
	title='count spectrum'
  stample
  ;	report
  message,'outputs are in:',/info
  message,'	count spectrum:	v_CTSPC[v_CHAN] [ct/bin]',/info
endif else begin
  ;	report
  message,'outputs are in:',/info
endelse
message,'	flux spectrum: v_FLXSPC[v_WVLS] [ct/bin]',/info
message,'	line spectrum: v_LINSPC[v_WVLS] [ph/s/cm^2/bin]',/info
message,'	continuum spectrum: v_CONSPC[v_WVLS] [ph/s/cm^2/bin]',/info
message,'	line fluxes: v_LINFLX[v_LWVL] [ph/s/cm^2]',/info
message,'	continuum fluxes: v_CONFLX[v_LWVL] [ph/s/cm^2]',/info

;	set the junk variables
jnk_LDBDIR=v_LDBDIR
jnk_CDBDIR=v_CDBDIR
jnk_CEROOT=v_CEROOT
jnk_IONEQF=v_IONEQF
jnk_EDENS=v_EDENS
jnk_WMIN=v_WMIN
jnk_WMAX=v_WMAX

end
