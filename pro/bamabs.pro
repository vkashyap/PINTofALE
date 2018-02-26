function bamabs,w,abund=abund,ikeV=ikeV,Fano=Fano,noHeH=noHeH,$
	verbose=verbose, _extra=e
;+
;function	bamabs
;	return photoelectric absorption cross-sections in energy range
;	30 eV - 10 keV using the polynomial fit coefficients determined
;	by Monika Balucinska-Church and Dan McCammon (Balucinska-Church,
;	M.\ and McCammon, D.\ 1992, ApJ 400, 699).  This is an update
;	to Morrison & McCammon 1983.
;
;syntax
;	sigabs=bamabs(w,abund=abund,/ikeV,/Fano,verbose=v)
;
;parameters
;	w	[INPUT; required] photon energy values at which to compute
;		the absorption
;		* default units are [Ang], unless IKEV is set, when they are
;		  assumed to be [keV]
;
;keywords
;	abund	[INPUT] abundances relative to H=1
;		* default is to use Anders & Grevesse 1982
;	ikeV	[INPUT] if set, assumes that W are in units of keV
;	Fano	[INPUT] if set, the 4 strongest auto-ionizing resonances
;		of He I are included; numbers come from Oza (1986, Phys
;		Rev A, 33, 824), and Fernley et al. (1987, J.Phys B 20, 6457).
;	noHeH	[INPUT] if set to any number other than 1 or 2, excludes
;		H and He from the cross-sections
;		* if set to 1, only excludes H
;		* if set to 2, only excludes He
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to avoid crashing the program
;
;subroutines
;	GETABUND
;	INICON
;
;restrictions
;	works only in energy range 30 eV - 10 keV
;	be warned that the edges are idealized and thus unrealistic
;	  for high resolution spectra because they do not include
;	  any resonance structure
;
;history
;	vinay kashyap (OctMM; based on TOTLXS.FOR and XSCTNS.FOR of
;	  Balucinska-Church & McCammon, obtained from ADC catalog VI/62A,
;	  ftp://adc.gsfc.nasa.gov/pub/adc/archives/catalogs/6/6062A/ )
;	bug corrections (Brian Kern, 5Dec2001): Ne cross-section was
;	  missing "^3"; Fano calc had EPSj[j] (VK; Dec2001)
;	force calculations only in valid energy range; added keyword
;	  NOHEH (VK; Jun03)
;	bug correction: "exclude=0" means "include" (VK; Jul03)
;	bug correction: lambda selection for He (Andy Ptak; VK Jun05)
;	edge lambda correction: O, N, C changed acc. to numbers in
;	  Atomic Data Booklet (2nd Edition), Thompson, A., et al., 2001,
;	  LBL/PUB-490 Rev. 2 (Jeremy Drake; VK Sep2006)
;-

;	usage
ok='ok' & np=n_params() & nw=n_elements(w)
if np eq 0 then ok='Insufficient parameters' else $
 if nw eq 0 then ok='W is undefined'
if ok ne 'ok' then begin
  print,'Usage: sigabs=bamabs(w,abund=abund,/ikeV,/Fano,noHeH=noHeH,verbose=v)'
  print,'  return photoelectric absorption cross-sections following'
  print,'  Balucinska-Church & McCammon 1992 update to Morrison & McCammon 1983'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	inputs
nrg=[float(abs(w))]
nab=n_elements(abund)
if nab lt 30 then abund=getabund('anders & grevesse')
if n_elements(ikeV) eq 0 then ikeV=0
if not keyword_set(ikeV) then begin	;convert from [Ang] to [keV]
  oo=where(nrg gt 0,moo)
  if moo gt 0 then nrg[oo]=12.3985/nrg[oo] else return,0*w
endif
v=0 & if keyword_set(verbose) then v=fix(abs(verbose[0])) > 0

exclH=0 & exclHe=0
if keyword_set(noHeH) then begin
  exclH=1 & exclHe=1
  if noHeH[0] eq 1 then exclHe=0	;only exclude H
  if noHeH[0] eq 2 then exclH=0		;only exclude He
endif

;	initialize
inicon,amu=amu,atom=atom,fundae=f
     ;H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Cl, Ar, Ca, Cr, Fe, Ni
useZ=[1,  2, 6, 7, 8, 10, 11, 12, 13, 14,16, 17, 18, 20, 24, 26, 28]
nZ=n_elements(useZ)
EE = nrg*1e3	;keV -> eV
;ok=where(EE gt 0 and EE le 1e4,mok) & elog=fltarr(nw)-90.
ok=where(EE ge 30 and EE le 1e4,mok) & elog=fltarr(nw)-90.
if mok gt 0 then elog[ok]=alog(EE[ok])

;	output
sig=0.*nrg	;photoelectric absorption cross-section in [cm^2]

for i=0L,nZ-1L do begin		;{step through encoded elements
  iZ=useZ[i] & aw=(amu.(iZ-1L))[0] & elem=atom[iZ-1L]
  case elem of

    'H': begin
      ;HYDRO (ANGIE BETKER, JEFF BLOCH, MBC 1991)
      hydro=fltarr(nw)
      if mok gt 0 then begin
        x=21.46941 + (3. - 2.060152)*Elog - 0.1492932*Elog*Elog +$
	  5.4634294e-3*(Elog^3)
        hydro[ok]=Exp(X[ok])/(EE[ok]^3)
      endif
      ;
      sig=sig + aw*hydro/f.NA*abund[iZ-1L]*(1-exclH)
    end

    'He': begin
      ;HELIUM (MBC 1997)	
      helium=fltarr(nw)
      C1=[-4.7416, 14.8200, -30.8678, 37.3584, -23.4585, 5.9133]
			;polynomial coefficients for Yan et al data
      IP=n_elements(C1)
      Q=[2.81, 2.51, 2.45, 2.44]
			;parameters Q for resonances (Fernley et al. 1987)
      EION=24.58	;ionization edge in eV
      NU=[1.610, 2.795, 3.817, 4.824]
			;parameters NU for resonances (Oza 1986)
      GAMMA=[2.64061E-3, 6.20116E-4, 2.56061E-4, 1.320159E-4]
			;parameters GAMMA for resonances (Oza 1986)
      ;
      lambda=fltarr(nw)+12398.54
      if mok gt 0 then begin	;(non-zero energies
        lambda[ok]=12398.54/EE[ok]
        X=EE/EION
	oo=where(lambda ge 1.239854 and lambda le 503.97,moo)
	;	was "or", which, as Andy Ptak points out, is always true
	y=fltarr(nw)+1. & sig_yan=fltarr(nw)
	if moo gt 0 then begin
	  for j=1L,ip do y[oo]=y[oo]+c1[j-1L]/(x[oo]^(j/2.))
	  sig_yan[oo]=733.E-24/(EE[oo]/1e3)^(3.5) * Y[oo]
	  if keyword_set(Fano) then begin	;(Fano profile correction
	    ;FANO (MBC 1993)
	    nf=n_elements(Q)
	    for j=0L,nf-1L do begin
	      eps = 911.2671/LAMBDA
	      epsj = 3.0 - 1./(NU[j]^2) + 1.807317
	      XX = 2.0*(EPS-EPSj)/GAMMA[j]
	      fano_factor=(XX-Q[j])^2/(1.+XX^2)
	      sig_yan = sig_yan * fano_factor
	    endfor
	  endif					;Fano)
	endif
	helium = sig_yan * f.NA/aw
      endif			;MOK>0)
      ;
      sig = sig + aw*helium/f.NA*abund[iz-1L]*(1-exclHe)
    end

    'C': begin
      ;CARBON
      carbon=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	;o1=where(EE gt 0 and EE lt 284.,mo1)
	;o2=where(EE ge 284. and EE le 1e4,mo2)
	o1=where(EE gt 0 and EE lt 284.2,mo1)	;JJD: cf. ADB rev2
	o2=where(EE ge 284.2 and EE le 1e4,mo2)	;JJD: cf. ADB rev2
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] =$
	  8.74161 + (7.13348*Elog[o1]) + (-1.14604*Elog[o1]^2) +$
	  (0.0677044*Elog[o1]^3)
	if mo2 gt 0 then X[o2] =$
	  3.81334 + (8.93626*Elog[o2]) + (-1.06905*Elog[o2]^2) +$
	  (0.0422195*Elog[o2]^3)
	carbon[ok]=exp(x[ok])/(ee[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*carbon/f.NA*abund[iz-1L]
    end

    'N': begin
      ;NITRO
      nitro=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	;o1=where(EE gt 0 and EE lt 401.,mo1)
	;o2=where(EE ge 401. and EE le 1e4,mo2)
	o1=where(EE gt 0 and EE lt 409.9,mo1)	;JJD: cf. ADB rev2
	o2=where(EE ge 409.9 and EE le 1e4,mo2)	;JJD: cf. ADB rev2
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] =$
	  9.24058 + (7.02985*Elog[o1]) + (-1.08849*Elog[o1]^2) +$
	  (0.0611007*Elog[o1]^3)
	if mo2 gt 0 then X[o2] =$
	  -13.0353 + (15.4851*Elog[o2]) + (-1.89502*Elog[o2]^2) +$
	  (0.0769412*Elog[o2]^3)
	;
	nitro[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*nitro/f.NA*abund[iz-1L]
    end

    'O': begin
      ;OXYGEN
      oxygen=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	;o1=where(EE gt 0 and EE lt 531.7,mo1)
	;o2=where(EE ge 531.7 and EE le 1e4,mo2)
	o1=where(EE gt 0 and EE lt 543.1,mo1)	;JJD: cf. ADB rev2
	o2=where(EE ge 543.1 and EE le 1e4,mo2)	;JJD: cf. ADB rev2
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] =$
	  2.57264 + (10.9321*Elog[o1]) + (-1.79383*Elog[o1]^2) +$
	  (0.102619*Elog[o1]^3)
	if mo2 gt 0 then X[o2] =$
	  16.53869 + (0.6428144 + 3.)*Elog[o2] - 0.3177744*Elog[o2]^2 +$
	  7.9471897e-3*(Elog[o2]^3)
	;
	oxygen[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*oxygen/f.NA*abund[iz-1L]
    end

    'Ne': begin
      ;NEON
      neon=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 867.,mo1)
	o2=where(EE ge 867. and EE le 1e4,mo2)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] =$
	  -3.04041 + (13.0071*Elog[o1]) + (-1.93205*Elog[o1]^2) +$
	  (0.0977639*Elog[o1]^3)
	if mo2 gt 0 then X[o2] =$
	  17.6007 + (3.29278*Elog[o2]) + (-0.263065*Elog[o2]^2) +$
	  (5.68290E-3*Elog[o2]^3)
	;
	neon[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*neon/f.NA*abund[iz-1L]
    end

    'Na': begin
      ;SODIUM
      sodium=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 1071.7,mo1)
	o2=where(EE ge 1071.7 and EE le 1e4,mo2)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] =$
	  -2737.598 + (2798.704 + 3.)*Elog[o1] - 1009.892*Elog[o1]^2 +$
	  87.16455*(Elog[o1]^3) + 43.20644*(Elog[o1]^4) -$
	  15.27259*(Elog[o1]^5) + 2.180531*(Elog[o1]^6) -$
	  0.1526546*(Elog[o1]^7) + 4.3137977E-3*(Elog[o1]^8)
	if mo2 gt 0 then X[o2] = $
	  1.534019 + (6.261744 + 3.)*Elog[o2] - 0.9914126*Elog[o2]^2 +$
	  3.5278253E-2*(Elog[o2]^3)
	;
	sodium[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*sodium/f.NA*abund[iz-1L]
    end

    'Mg': begin
      ;MAGNES
      magnes=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 49.45,mo1)
	o2=where(EE ge 49.45 and EE lt 1303.4,mo2)
	o3=where(EE ge 1303.4 and EE le 1e4,mo3)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] = 7.107172 + (0.7359418 + 3.)*Elog[o1]
	if mo2 gt 0 then X[o2] = $
	  -81.32915 + (62.2775 + 3.)*Elog[o2] - 15.00826*Elog[o2]^2 +$
	  1.558686*(Elog[o2]^3) - 6.1339621E-2*(Elog[o2]^4)
	if mo3 gt 0 then X[o3] = $
	  -9.161526 + (10.07448 + 3.)*Elog[o3] - 1.435878*Elog[o3]^2 +$
	  5.2728362E-2*(Elog[o3]^3)
	;
	magnes[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*magnes/f.NA*abund[iz-1L]
    end

    'Al': begin	;ALUM
      alum=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 72.78,mo1)
	o2=where(EE ge 72.78 and EE lt 1559.9,mo2)
	o3=where(EE ge 1559.9 and EE le 1e4,mo3)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] =$
	  26.90487 + (3. - 9.135221)*Elog[o1] + 1.175546*Elog[o1]^2
	if mo2 gt 0 then X[o2] = -38.1232 + 29.5161*Elog[o2] -$
	  4.45416*Elog[o2]^2 + 0.226204*Elog[o2]^3
	if mo3 gt 0 then X[o3] =$
	  14.6897 + 4.22743*Elog[o3] - 0.344185*Elog[o3]^2 +$
	  8.18542E-3*Elog[o3]^3
	;
	alum[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*alum/f.NA*abund[iz-1L]
    end

    'Si': begin
      ;SILICN
      silicn=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 100.6,mo1)
	o2=where(EE ge 100.6 and EE lt 1840.0,mo2)
	o3=where(EE ge 1840.0 and EE le 1e4,mo3)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] = $
	  -3.066295 + (7.006248 + 3.)*Elog[o1] - 0.9627411*Elog[o1]^2
	if mo2 gt 0 then X[o2] = $
	  -182.7217 + (125.061 + 3.)*Elog[o2] - 29.47269*Elog[o2]^2 +$
	  3.03284*(Elog[o2]^3) - 0.1173096*(Elog[o2]^4)
	if mo3 gt 0 then X[o3] = $
	  -33.39074 + (18.42992 + 3.)*Elog[o3] -$
	  2.385117*Elog[o3]^2 + 8.887583e-2*(Elog[o3]^3)
	;
	silicn[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*silicn/f.NA*abund[iz-1L]
    end

    'S': begin	;SULFUR
      sulfur=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 165.0,mo1)
	o2=where(EE ge 165.0 and EE lt 2470.5,mo2)
	o3=where(EE ge 2470.5 and EE le 1e4,mo3)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] = $
	  598.2911 + (3. - 678.2265)*Elog[o1] + 308.1133*Elog[o1]^2 -$
	  68.99324*(Elog[o1]^3) + 7.62458*(Elog[o1]^4) - 0.3335031*(Elog[o1]^5)
	if mo2 gt 0 then X[o2] = $
	  3994.831 + (3. - 3693.886)*Elog[o2] + 1417.287*Elog[o2]^2 -$
	  287.9909*(Elog[o2]^3) + 32.70061*(Elog[o2]^4) -$
	  1.968987*(Elog[o2]^5) + 4.9149349E-2*(Elog[o2]^6)
	if mo3 gt 0 then X[o3] = $
	  -22.49628 + (14.24599+ 3.)*Elog[o3] - 1.848444*Elog[o3]^2 +$
	  6.6506132E-2*(Elog[o3]^3)
	;
	sulfur[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*sulfur/f.NA*abund[iz-1L]
    end

    'Cl': begin
      ;CHLORN
      chlorn=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 202.0,mo1)
	o2=where(EE ge 202.0 and EE lt 2819.6,mo2)
	o3=where(EE ge 2819.6 and EE le 1e4,mo3)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] = $
	  6253.247 +(3. - 8225.248)*Elog[o1] + 4491.675*Elog[o1]^2 -$
	  1302.145*(Elog[o1]^3) + 211.4881*(Elog[o1]^4) -$
	  18.25547*(Elog[o1]^5) + 0.6545154*(Elog[o1]^6)
	if mo2 gt 0 then X[o2] = $
	  -233.0502 + (143.9776 + 3.)*Elog[o2] - 31.12463*Elog[o2]^2 +$
	  2.938618*(Elog[o2]^3) - 0.104096*(Elog[o2]^4)
	if mo3 gt 0 then X[o3] = $
	  -23.74675 + (14.50997 + 3.)*Elog[o3] - 1.857953*Elog[o3]^2 +$
	  6.6208832E-2*(Elog[o3]^3)
	;
	chlorn[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*chlorn/f.NA*abund[iz-1L]
    end

    'Ar': begin
      ;ARGON
      argon=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 245.0,mo1)
	o2=where(EE ge 245.0 and EE lt 3202.9,mo2)
	o3=where(EE ge 3202.9 and EE le 1e4,mo3)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] = $
	  -330.3509 + (267.7433 + 3.)*Elog[o1] - 78.90498*Elog[o1]^2 +$
	  10.35983*(Elog[o1]^3) - 0.5140201*(Elog[o1]^4)
	if mo2 gt 0 then X[o2] = $
	  -5.71870 + (8.85812*Elog[o2]) + (-0.307357*Elog[o2]^2) +$
	  (0.00169351*(Elog[o2]^3)) + (-0.0138134*(Elog[o2]^4)) +$
	  (0.00120451*(Elog[o2]^5))
	if mo3 gt 0 then X[o3] = $
	  19.1905 + (2.74276*Elog[o3]) + (-0.164603*Elog[o3]^2) +$
	  (0.00165895*Elog[o3]^3)
	;
	argon[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*argon/f.NA*abund[iz-1L]
    end

    'Ca': begin
      ;CALC
      calci=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 349.31,mo1)
	o2=where(EE ge 349.31 and EE lt 4038.1,mo2)
	o3=where(EE ge 4038.1 and EE le 1e4,mo3)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] = $
	  -873.972 + (865.5231 + 3.)*Elog[o1] - 339.678*Elog[o1]^2 +$
	  66.83369*(Elog[o1]^3) - 6.590398*(Elog[o1]^4) + 0.2601044*(Elog[o1]^5)
	if mo2 gt 0 then X[o2] = $
	  -3449.707 + (2433.409 + 3.)*Elog[o2] - 682.0668*Elog[o2]^2 +$
	  95.3563*(Elog[o2]^3) - 6.655018*(Elog[o2]^4) + 0.1854492*(Elog[o2]^5)
	if mo3 gt 0 then X[o3] = $
	  18.89376 + (3. - 0.2903538)*Elog[o3] - 0.1377201*Elog[o3]^2
	;
	calci[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*calci/f.NA*abund[iz-1L]
    end

    'Cr': begin
      ;CHROM
      chrom=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 598.0,mo1)
	o2=where(EE ge 598.0 and EE lt 691.0,mo2)
	o3=where(EE ge 691.0 and EE lt 5988.8,mo3)
	o4=where(EE ge 5988.8 and EE le 1e4,mo4)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] = $
	  -0.4919405 + (12.66939 + 3.)*Elog[o1] - 5.199775*Elog[o1]^2 +$
	  1.086566*(Elog[o1]^3) - 0.1196001*(Elog[o1]^4) +$
	  5.2152011E-3*(Elog[o1]^5)
	if mo2 gt 0 then X[o2] = 27.29282 +(3. - 2.703336)*Elog[o2]
	if mo3 gt 0 then X[o3] = -15.2525 + (13.23729 + 3.)*Elog[o3] -$
	  1.966778*Elog[o3]^2 + 8.062207E-2*(Elog[o3]^3)
	if mo4 gt 0 then X[o4] = 8.307041 + (2.008987 + 3.)*Elog[o4] -$
	  0.2580816*Elog[o4]^2
	;
	chrom[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*chrom/f.NA*abund[iz-1L]
    end

    'Fe': begin
      ;IRON
      iron=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 707.4,mo1)
	o2=where(EE ge 707.4 and EE lt 7111.2,mo2)
	o3=where(EE ge 7111.2 and EE le 1e4,mo3)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] = $
	  -15.07332 + (18.94335 + 3.)*Elog[o1] - 4.862457*Elog[o1]*Elog[o1] +$
	  0.5573765*(Elog[o1]^3) - 3.0065542E-2*(Elog[o1]^4) +$
	  4.9834867E-4*(Elog[o1]^5)
	if mo2 gt 0 then X[o2] = $
	  -253.0979 + (135.4238 + 3.)*Elog[o2] - 25.47119*Elog[o2]*Elog[o2] +$
	  2.08867*(Elog[o2]^3) - 6.4264648E-2*(Elog[o2]^4)
	if mo3 gt 0 then X[o3] = $
	  -1.037655 + (4.022304 + 3.)*Elog[o3] - 0.3638919*Elog[o3]*Elog[o3]
	;
	iron[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*iron/f.NA*abund[iz-1L]
    end

    'Ni': begin
      ;NICKEL
      nickel=fltarr(nw)
      if mok gt 0 then begin	;(non-zero energies
	o1=where(EE gt 0 and EE lt 853.6,mo1)
	o2=where(EE ge 853.6 and EE lt 8331.6,mo2)
	o3=where(EE ge 8331.6 and EE le 1e4,mo3)
	X=fltarr(nw)-200.
	if mo1 gt 0 then X[o1] = $
	  -7.919931 + (11.06475 + 3.)*Elog[o1] - 1.935318*Elog[o1]^2 +$
	  9.3929626e-2*(Elog[o1]^3)
	if mo2 gt 0 then X[o2] = $
	  3.71129 + (8.45098*Elog[o2]) + (-0.896656*Elog[o2]^2) +$
	  (0.0324889*Elog[o2]^3)
	if mo3 gt 0 then X[o3] = 28.4989 + (0.485797*Elog[o3])
	;
	nickel[ok] = EXP(X[ok])/(EE[ok]^3)
      endif			;MOK>0)
      ;
      sig = sig + aw*nickel/f.NA*abund[iz-1L]
    end

    else: message,elem+': cross-sections not available',/info
  endcase
  if v gt 0 then begin
    if nw eq 1 then print,elem+': ',sig[0],' @'+strtrim(nrg[0],2)
    if v gt 10 and mok gt 1 then plot,w[ok],sig[ok],/yl,title=elem,psym=3
  endif
endfor				;i=0,nZ-1}

return,sig
end
