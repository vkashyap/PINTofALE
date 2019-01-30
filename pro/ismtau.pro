function ismtau,w,NH=NH,fH2=fH2,He1=He1,He2=He2,Fano=Fano,ikev=ikev,$
    vion=vion,wam=wam,bam=bam,mam=mam, noHeH=noHeH,$
    tauH1=tauH1,tauH2=tauH2,tauHe2=tauHe2,tauHe1=tauHe1,tauX=tauX,tauTB=tauTB,$
    tauFano=tauFano,tauFM=tauFM,$
    icrstr=icrstr,$
    EBV=EBV,R_V=R_V,LMC2=LMC2,AVGLMC=AVGLMC,ExtCurve=ExtCurve,$
    fmgamma=fmgamma,fmx0=fmx0,fmc1=fmc1,fmc2=fmc2,fmc3=fmc3,fmc4=fmc4,$
    HeII=HeII,tauHeII=tauHeII,$
    _extra=e
;+
;function	ismtau
;	returns the optical depth for a given photon energy at specified
;	column density.  i.e., the $\tau$ of $e^{-\tau}$.
;
;	ingredients:
;	* in the UV/Opt/IR regime, calls IDLastro's FM_UNRED.PRO
;	* in the EUV regime, cannibalizes ISMEUV.PRO for
;	  -- HI photoionization (Spitzer 19??, "Physical Processes in the
;	     Interstellar Medium", p105), between 0.03 keV (413.3A) and 911.75A
;	  -- HeII photoionization (ibid), between 0.03 keV (413.3A) and 911.75A
;	  -- HeI with auto-ionization resonances (Marr & West 1976;
;	     Oza 1986, Phys.Rev.A 33 824; Fernley et al. 1987, J.Phys.B 20,
;	     6457; Rumph, Bowyer, & Vennes 1994, AJ 107, 2108) below 503.97A
;	  -- NOTE: these are extended down to 43.657 AA when the X-ray
;	     calculations are not being done.
;	* in the X-ray regime (0.03 - 10 keV), either
;	  -- Morrison & McCammon (1983, ApJ 270, 119) polynomial fit between
;	     0.03 keV and 10 keV, or
;	  -- Balucinska-Church & McCammon (1992, ApJ 400, 699) cross-sections,
;	     which allows abundance variations, or
;	  -- Wilms, Allen, & McCray (2000, ApJ, in press) calculations
;	     that use cross-sections from Verner & Yakovlev (1995, A&AS,
;	     109, 125), which also allows abundance variations
;	  -- Verner, Ferland, Korista, & Yakovlev (1996, ApJ, 465, 487)
;	     ground state photoionisation cross-sections that allows for
;	     abundance variations (down to 911.75 Ang)
;	* in the hard X-ray regime and beyond (>10 keV),
;	  -- Tanaka & Bleeker (1977, ???) power-law extended past 10 keV
;	* Molecular H (Cruddace et al. 1977, ApJ 187, 497; Kashyap et
;	  al. 1992, ApJ 391, 684) between 12.3985A (1 keV) to 911.75A
;
;	to extract only the cross-sections, set NH=1
;	to obtain transmission factors, use exp(-ismtau(...))
;
;syntax
;	tau=ismtau(w,NH=NH,fH2=fH2,He1=He1,He2=He2,Fano=Fano,/ikev,$
;	/vion,/wam,/bam,/mam, abund=abund,noHeH=noHeH,verbose=verbose,$
;	tauH1=tauH1,tauH2=tauH2,tauHe2=tauHe2,tauHe2=tauHe2,tauHe1=tauHe1,$
;	tauX=tauX,tauTB=tauTB,tauFano=tauFano,tauFM=tauFM,$
;	R_V=R_V,LMC2=LMC2,AVGLMC=AVGLMC,ExtCurve=ExtCurve,$
;       fmgamma=fmgamma,fmx0=fmx0,fmc1=fmc1,fmc2=fmc2,fmc3=fmc3,fmc4=fmc4,$
;	icrstr=icrstr,ionfrac=ionfrac,vfkydir=vfkydir,nconsec=nconsec,$
;	DEM=DEM,logT=logT,chidir=chidir,eqfile=eqfile)
;
;parameters
;	w	[INPUT; required] photon energy values at which to compute
;		the optical depth.
;		* default units are [Ang], unless IKEV is set, when they are
;		  assumed to be [keV]
;
;keywords
;	NH	[INPUT; default=1] atomic H I column density [cm^-2]
;		* when keyword VION is set, IONABS() is called to deal
;		  with partially ionized media.  In this case, NH is
;		  assumed to be total H column, i.e., H I + protons
;		* note that when NH=1, the function essentially returns
;		  the absorption cross-section in units of [cm^2]
;		* to set the "column" for UV/optical/IR, use the keyword
;		  EBV, which defines the E(B-V) and gets passed down to
;		  the routine FM_UNRED
;		  - this means that the "feature" of getting the cross-section
;		    out for NH=1 that is available for high-energy photons is
;		    not available at wavelengths > 911.75 Ang
;	fH2	[INPUT; default=0.] N(H2)/N(HI)
;		* originally was set to a default of 0.26, the Galactic
;		  average -- comes from the asymptotic value of
;		  N(H2)/N(HI) computed with the 3 component model (cold
;		  HI, warm HI, H2) of Bloemen, J.B.G.M. (1987, ApJ 322, 694)
;		  -- but has now been changed to 0. explicitly set it to a
;		  float number to include H2 cross-sections
;		* if NoHeH is set to exclude H, fH2 is assumed to be 0
;	He1	[INPUT; default=0.1*NH] neutral He column density [cm^-2]
;	He2	[INPUT; default=0.1*He1] ionized He column density [cm^-2]
;		* if NoHeH is set to exclude He, both HeI and HeII cross-sections
;		  will be excluded regardless of what the keywords He1 and
;		  He2 are set to
;	HeII	[OBSOLETE; same as HE2] kept around for backwards compatibility
;	Fano	[INPUT] if set, the 4 strongest auto-ionizing resonances
;		of He I are included; the shape of these resonances are
;		given by a Fano profile (cf. Rumph, Bowyer, & Vennes 1994,
;		AJ 107, 2108).
;		* input wavelength grid between 190 A and 210 A should be
;		  ~0.01 A for this to be useful
;	ikev	[INPUT] if set, assumes that W are in units of keV
;	vion	[INPUT] if set, uses Verner et al. (1996) cross-sections
;		* calls IONABS()
;	wam	[INPUT] if set, uses Wilms, Allen, & McCray cross-sections
;		* NOT IMPLEMENTED YET -- if set, defaults to Balucinska-Church
;		  & McCammon.
;	bam	[INPUT] if set, uses Balucinska-Church and McCammon updates
;		to Morrison & McCammon cross-sections in 30 eV - 10 keV range
;		* calls BAMABS()
;	mam	[INPUT] if set, uses Morrison and McCammon photoionization
;		cross-sections in the 30 eV - 10 keV range
;		* This is the default
;		* VION takes precedence over
;		  WAM takes precedence over
;		  BAM takes precedence over MAM.
;		* if MAM is explicitly set to 0 and none of VION, BAM, or WAM are
;		  set, then the original ISMEUV code is used from 43 AA longwards
;	noHeH	[INPUT] passed on unchanged to BAMABS() and IONABS()
;		if set to any number other than 1 or 2, excludes H and He
;		cross-sections in BAMABS (and ensures that even if fH2 and
;		He1 are not zero, they are properly excluded within ISMTAU)
;		* if set to 1, excludes only H
;		* if set to 2, excludes only He
;		* NOTE: if set, will supersede fH2, He1 and HeII keyword
;		  values in that H, H2, He1, He2 are all appropriately
;		  excluded
;	tauH1	[OUTPUT] optical depth due to H
;	tauH2	[OUTPUT] optical depth due to molecular H2
;	tauHe2	[OUTPUT] optical depth due to HeII
;	tauHeII	[OBSOLETE] same as TAUHE2
;		* both present here for backwards compatibility
;	tauHe1	[OUTPUT] optical depth due to HeI
;		* NOTE: if BAM or VION is set, BAMABS() and IONABS() compute
;		  H, HeI, and HeII cross-sections, so the parts of tauH1,
;		  tauHe1, tauHe2 that overlap the relevant range get set to zero
;	tauX	[OUTPUT] optical depth in EUV and X-ray regime (0.03-10 keV)
;	tauTB	[OUTPUT] optical depth in hard X-ray regime (>10 keV)
;	tauFano	[OUTPUT] Fano correction to optical depth due to HeI
;		* NOTE: this is included after VION
;	tauFM	[OUTPUT] optical depth due to Fitzpatrick & Massa reddening
;	icrstr	[I/O] passed without comment to IONABS()
;	EBV	[INPUT] color excess, E(B-V), passed on to FM_UNRED
;		* default is 0
;		* note that -ve values are not forbidden in FM_UNRED,
;		  but are ignored here because we don't want any
;		  confusion with the other optical depths
;	R_V	[INPUT] passed without comment to FM_UNRED
;	LMC2	[INPUT] passed without comment to FM_UNRED
;	AVGLMC	[INPUT] passed without comment to FM_UNRED
;	ExtCurve[INPUT] passed without comment to FM_UNRED
;	fmgamma	[INPUT] passed without comment to FM_UNRED
;	fmx0	[INPUT] passed without comment to FM_UNRED
;	fmc1	[INPUT] passed without comment to FM_UNRED
;	fmc2	[INPUT] passed without comment to FM_UNRED
;	fmc3	[INPUT] passed without comment to FM_UNRED
;	fmc4	[INPUT] passed without comment to FM_UNRED
;		* note: all of these keywords are explicitly included
;		  because FM_UNRED doesn't have an _EXTRA keyword
;	_extra	[INPUT ONLY] pass defined variables to subroutines
;		BAMABS: ABUND,VERBOSE
;		IONABS: ABUND,IONFRAC,VFKYDIR,NCONSEC,DEM,LOGT,CHIDIR,EQFILE
;
;restrictions
;	requires BAMABS and IONABS (part of PINTofALE) and
;	FM_UNRED and associated subroutines (from IDLastro)
;
;history
;	written by Vinay Kashyap (Feb97), based on
;	  ISMEUV.PRO (W.Landsman [Oct94], via ISM.C@cea-ftp.cea.berkeley.edu
;	  (Pat Jelinsky, Todd Rumph, et al.)) and ISMABS.PRO (VK; Dec91)
;	changed default on HeII from 0 to 0.1*He1 (VK; Jan MM)
;	changed default on NH from 1e20 to 1.0 (VK; MarMM)
;	changed keyword KEV to IKEV; also for some reason H and He
;	  photoionization from ISMEUV had been cut off at the C edge,
;	  now corrected; added keywords WAM, BAM, and MAM, and call to BAMABS
;	  (VK; OctMM)
;	bug correction: BAM was not filtering on energy; allowed setting
;	  MAM=0 to recover original ISMEUV calc (VK; Jun03)
;	follow-up bug correction (VK; Jul03)
;	added keyword NOHEH to streamline behavior of handling H and He
;	  cross-sections with BAMABS (thanks Christian Schneider; VK Feb07)
;	added keywords TAUH1,TAUH2,TAUHE2,TAUHE1,TAUX,TAUTB, VION,ICRSTR,;
;	  changed all array index notation from old "()" to new "[]";
;	  **WARNING** changed behavior of FH2 -- new default is 0, not 0.26;
;	  added call to IONABS();
;	  (VK; May09)
;	added call to FM_UNRED for wvl > 911.75 Ang (VK; Aug11)
;	bug fixes (Jon Slavin; Sep13):
;	  /vion was double counting tauH,tauHe1,tauHe2
;	  (WARNING: changes behavior of tauH, tauHe1,tauHe2 on output)
;	added keyword TAUHE2 (VK; Sep13)
;	added keyword He2 as alternative for HeII; added keywords tauFano and
;	  tauFM; reincluded Fano correction after /vion; zeroed out He opacities
;	  after call to BAMAMS (VK; Jan19)
;-

;	usage
nw=n_elements(w)
if nw eq 0 then begin
  print,'Usage: tau=ismtau(W,NH=NH,fH2=fH2,He1=He1,He2=He2,/Fano,/ikev,$'
  print,'       /vion,/wam,/bam,/mam, noHeH=noHeH, abund=abund,verbose=verbose,$'
  print,'       tauH1=tauH1,tauH2=tauH2,tauHe2=tauHe2,tauHe1=tauHe1,$'
  print,'       tauX=tauX,tauTB=tauTB, icrstr=icrstr,ionfrac=ionfrac,$'
  print,'       EBV=EBV,R_V=R_V,LMC2=LMC2,AVGLMC=AVGLMC,ExtCurve=ExtCurve,$'
  print,'       fmgamma=fmgamma,fmx0=fmx0,fmc1=fmc1,fmc2=fmc2,fmc3=fmc3,fmc4=fmc4,$
  print,'       vfkydir=vfkydir,nconsec=nconsec,DEM=DEM,logT=logT,$'
  print,'       chidir=chidir,eqfile=eqfile)'
  print,'  returns optical depth at given photon energy for given column'
  return,0.
endif

nrg=[1.0*(abs(w))]	;avoid touching the input, convert scalar to vector,
			;avoid -ve values in wvl/keV, ensure at least R*4

;	check keywords
if n_elements(NH) eq 0 then NH=1.			;N(HI)
 ;if n_elements(fH2) eq 0 then fH2=0.26			;N(H2)/N(HI)
if not keyword_set(fH2) then fH2=0.			;N(H2)/N(HI)
if n_elements(He1) eq 0 then He1=0.1*NH			;N(HeI)
if n_elements(He2) eq 0 then begin			;N(HeII)
  if n_elements(HeII) ne 0 then He2=HeII[0] else He2=0.1*He1
endif
;	if input photon "energy" is in [Ang], convert to [keV]
if not keyword_set(ikeV) then begin
  wvl=nrg & oo=where(nrg gt 0,moo)
  if moo gt 0 then nrg[oo]=12.3985/nrg[oo] else return,0*w
endif else begin
  wvl=0*nrg & oo=where(nrg gt 0,moo)
  if moo gt 0 then wvl[oo]=12.3985/nrg[oo] else return,0*w
endelse
;
iMAM=1 & iBAM=0 & iWAM=0 & iVION=0
if n_elements(MAM) gt 0 then begin
  if MAM[0] eq 0 then begin
    iMAM=0 & iBAM=0 & iWAM=0 & iVION=0
  endif
endif
if keyword_set(BAM) then begin
  iMAM=0 & iBAM=1 & iWAM=0 & iVION=0
endif
if keyword_set(WAM) then begin
  message,'Wilms, Allen, McCray cross-sections not implemented yet',/info
  message,'using Balucinska-Church & McCammon data',/info
  iMAM=0 & iBAM=1 & iWAM=0 & iVION=0
endif
if keyword_set(VION) then begin
  iMAM=0 & iBAM=0 & iWAM=0 & iVION=1
endif
exclH=0 & exclHe=0
if keyword_set(noHeH) then begin
    exclH=1 & exclHe=1
    if noHeH[0] eq 1 then exclHe=0	;only exclude H
    if noHeH[0] eq 2 then exclH=0	;only exclude He
endif

;	initialize
tauH=0*wvl & tauH2=tauH & tauHeII=tauH & tauHeI=tauH & tauX=tauH & tauTB=tauH
tauFano=tauH & tauFM=tauH

;	HI photoionization (from ISMEUV)
ev30 = 12.3985/(30.e-3)
	;Ang corresponding to 30 eV, boundary of ISMEUV applicability
	;when WAM, BAM, or MAM is set
if iMAM eq 0 and iBAM eq 0 and iWAM eq 0 then ev30=43.657
r=wvl/911.75
;DO NOT USE THIS: oo=where(r ge 43.657/911.75 and r lt 1,moo)
oo=where(r ge ev30/911.75 and r lt 1,moo)
if moo gt 0 then begin
  z=0*r & z[oo]=sqrt(r[oo]/(1.-r[oo]))
  tauH[oo] = NH[0] * 3.44e-16 * (r[oo]^4) *$
    exp(-4.0*z[oo]*atan(1./z[oo]))/(1.0-exp(-2*!pi*z[oo]))
endif

;	HeII photoionization (from ISMEUV); just like above, but charge=2
r=4*wvl/911.75
;WRONG-> oo=where(r gt 43.657/911.75 and r lt 1,moo)
oo=where(r ge ev30/911.75 and r lt 1,moo)
if moo gt 0 then begin
  z=0*r & z[oo]=sqrt(r[oo]/(1.-r[oo]))
  tauHeII[oo]=He2[0]*3.44e-16*(r[oo]^4)*$
    exp(-4.0*z[oo]*atan(1./z[oo]))/((1.0-exp(-2*!pi*z[oo]))*4)
endif

;	HeI (from ISMEUV)
;	for wvl>46A (Marr & West 1976)
c1=[ -2.953607d+01, 7.083061d+00,  8.678646d-01, -1.221932d+00,$
      4.052997d-02, 1.317109d-01, -3.265795d-02,  2.500933d-03 ]
;	for wvl<46A (Marr & West 1976)
c2=[ -2.465188d+01, 4.354679d+00, -3.553024d+00, 5.573040d+00,$
     -5.872938d+00, 3.720797d+00, -1.226919d+00, 1.576657d-01 ]
;	parameters for auto-ionization resonances
;	from Oza (1986, Phys.Rev.A 33, 824)
q=[ 2.81d, 2.51d, 2.45d, 2.44d ]
;	from Fernley et al. (1987, J.Phys.B 20,6457)
nu=[ 1.610d, 2.795d, 3.817d, 4.824d ]
fano_gamma=[ 2.64061d-03, 6.20116d-04, 2.56061d-04, 1.320159d-04 ]
esubi=3.0d - 1.0d/nu^2 + 1.807317d
;
;INCORRECT-> oo=where(wvl gt 0 and wvl lt 503.97,moo)
oo=where(wvl gt ev30 and wvl lt 503.97,moo) & x=0*wvl & y=x & yf=y & eps=x
if moo gt 0 then begin
  x[oo]=alog10(wvl[oo])
  o1=where(wvl ge 46.0 and wvl lt 503.97,mo1)
  o2=where(wvl gt 0 and wvl lt 46.0,mo2)
  if mo1 gt 0 then y[o1]=poly(x[o1],c1)
  if mo2 gt 0 then y[o2]=poly(x[o2],c2)
  if keyword_set(Fano) then begin
    eps[oo]=911.2671/wvl[oo]
    for i=0,3 do begin
      x[oo]=2.0*((eps[oo]-esubi[i])/fano_gamma[i])
      yf[oo]=yf[oo]+alog10((x[oo]-q[i])^2/(1.+x[oo]*x[oo]))
    endfor
  endif
  tauHeI[oo]=He1[0] * 10.^(y[oo])
  tauFano[oo]=10.^(yf[oo])
endif

;	H2 (from ISMABS)
if not keyword_set(iwam) then begin
  ; Wilms, Allen & McCray have their own computation of H2 cross-section
  ; so don't use this if WAM is set
  oo=where(nrg gt 0 and nrg lt 1,moo)
  if moo gt 0 then tauH2[oo]=NH[0]*abs(fh2[0])*0.16e-22*nrg[oo]^(-3.54)
endif

;	X-Ray regime (from ISMABS)

if keyword_set(iMAM) then begin
  ;	Morrison & McCammon (1983, ApJ 270, 119) fit coefficients
  c0 = [ 17.3,      34.6,  78.1,  71.4,  95.5,  308.9, 120.6, 141.3, 202.7,$
       342.7,    352.2, 433.9, 629.0, 701.2 ]
  c1 = [ 608.1,    267.9,  18.8,  66.8, 145.8, -380.6, 169.3, 146.8, 104.7,$
       18.7,      18.7,  -2.4,  30.9,  25.2 ]
  c2 = [ -2150.0, -476.1,   4.3, -51.4, -61.1,  294.0, -47.7, -31.5, -17.0,$
       0.0,        0.0,  0.75,   0.0,   0.0 ]
  el = [ 0.030,    0.100, 0.284, 0.400, 0.532,  0.707, 0.867, 1.303, 1.840,$
       2.471,    3.210, 4.038, 7.111, 8.331 ]
  eh = [ 0.100,    0.284, 0.400, 0.532, 0.707,  0.867, 1.303, 1.840, 2.471,$
       3.210,    4.038, 7.111, 8.331, 10.00 ]
  nbin=n_elements(el) & sig0=1e-24
  for i=0,nbin-1 do begin
    oo=where(nrg gt el[i] and nrg le eh[i],moo)
    if moo gt 0 then tauX[oo]=NH[0]*sig0*$
      (c0[i]+c1[i]*nrg[oo]+c2[i]*nrg[oo]^2)/(nrg[oo])^3
  endfor
  oy=where(nrg ge 0.0136 and nrg lt 0.03,moy)
  if moy gt 0 then begin
    message,'Morrison & McCammon X-ray opacities between 13.6 to 300 eV set to 0',/informational
  endif
endif

if keyword_set(iBAM) then begin
  ;	Balucinska-Church & McCammon, 1992, ApJ 400, 699
  oo=where(nrg ge 0.03 and nrg le 10.,moo)
  if moo gt 0 then tauX[oo]=bamabs(nrg[oo],/ikeV,Fano=Fano,noHeH=noHeH, _extra=e)*NH[0]

  oy=where(nrg ge 0.0136 and nrg lt 0.03,moy)
  if moy gt 0 then begin
    message,'Balucinska-Church & McCammon X-ray opacities between 13.6 to 300 eV set to 0',/informational
  endif

  ;	BUG FIX? Jeremy Drake (2014-jul-30; same as the one pointed out by Jon Slavin 2013-sep-12)
  ;	these are already included in bamabs (when noHeH is not set [and when it is set,
  ;	makes no difference]), so no need to add them in again below
  if moo gt 0 then tauH[oo]=0.
  if moo gt 0 then tauHeI[oo]=0.
  if moo gt 0 then tauHeII[oo]=0.
endif

if keyword_set(iVION) then begin
  ;	Verner, Ferland, Korista, & Yakovlev, 1996, ApJ, 465, 487
  oo=where(nrg ge 0.0136 and nrg le 10.,moo)
  if moo gt 0 then tauX[oo]=ionabs(nrg[oo],/ikeV,icrstr=icrstr,noHeH=noHeH, _extra=e)*NH[0]

  ;	BUG FIX: Jon Slavin (2013-sep-12)
  ;	these are already included in ionabs (when noHeH is not set [and when it is set,
  ;	makes no difference]), so no need to add them in again below
  if moo gt 0 then tauH[oo]=0.
  if moo gt 0 then tauHeI[oo]=0.
  if moo gt 0 then tauHeII[oo]=0.
  if moo gt 0 then tauX[oo]=tauX[oo]+(1.-exclHe)*He1[0]*tauFano[oo]	;include Fano correction because VION does not include it
endif

;	Hard X-ray regime (from ISMABS)
;	Tanaka & Bleeker coefficients
c0=2e-22 & c1=-2.5 & oo=where(nrg gt 10,moo)
if moo gt 0 then tauTB[oo]=NH[0]*c0*(nrg[oo])^(c1)

;	UV/optical/IR (call FM_UNRED)
oo=where(nrg gt 0 and nrg le 12.3985/911.75,moo)
if moo gt 0 then begin
  ;if IDLastro is in the path (I don't know how to test for that generally so ignoring it for now)
    wave=12.3985/nrg[oo] & flux=fltarr(moo)+1.
    if not keyword_set(EBV) then zEBV=0. else zEBV=ebv[0] > 0
    fm_unred,wave,flux,zEBV,funred,$
  	R_V=R_V,LMC2=LMC2,AVGLMC=AVGLMC,ExtCurve=ExtCurve,$
	gamma=fmgamma,x0=fmx0,c1=fmc1,c2=fmc2,c3=fmc3,c4=fmc4
    tauFM[oo]=-alog(flux/funred)
  ;endif)
endif

;	add up ever'thin
tau = (1.-exclH)*tauH + (1.-exclH)*tauH2 + (1.-exclHe)*tauHeII + (1.-exclHe)*tauHeI*tauFano + tauX + tauTB + tauFM

;	strictly for keyword compatibility
;	because IDL doesn't understand that tauHeII and tauHeI are different keyword names
tauHe2=tauHeII
tauHe1=tauHeI

return,tau
end
