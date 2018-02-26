pro wrt_cont_cie,elem,pres,logT,wvl,reH,n_e=n_e,wrange=wrange,$
	outdir=outdir,root=root,eqfile=eqfile,abund=abund, _extra=e
;+
;procedure	wrt_cont_cie
;	writes out continuum emissivities (NLOGT,NBIN) using CIE (subset
;	of SPEX) [1e-23 ergs cm^3/s/A]
;
;	WARNING: The continuum emissivities differ from the line
;	emissivities in that abundances and ion balance are *included*!
;
;usage
;	wrt_cont_cie,elem,pres,logT,wvl,reH,n_e=n_e,wrange=wrange,$
;	outdir=outdir,root=root,wmn=wmn,wmx=wmx,nbin=nbin,abund=abund,$
;	ciedir=ciedir,eqfile=eqfile,chidir=chidir
;
;parameters
;	elem	[INPUT; required] name of element e.g., "He", "C", etc.
;		* used as prefix for datafiles
;		* passed w/o comment to CONT_CIE
;	pres	[INPUT; default: 1e15 cm^-3 K] electron pressure at which
;		to compute intensities
;		* passed w/o comment to CONT_CIE
;	logT	[INPUT] array of log10(Temperature [K]) at which
;		to compute intensities.
;		* passed w/o comment to CONT_CIE
;	wvl	[OUTPUT] all the bin-beginning values and the final bin
;		ending value in the spectrum [Ang]
;		* passed w/o comment to CONT_CIE
;	reH	[OUTPUT] n_e/n_H for each LOGT
;		* passed w/o comment to CONT_CIE
;
;keywords
;	n_e	[INPUT] electron density [/cm^3]
;		* OVERRIDES values determined using PRES and LOGT
;		* if set, appends a "D" to OUTDIR
;	wrange	[INPUT] wavelength range over which to compute the spectrum
;		* default: [1.,180.] A
;		* min(WRANGE)=0.1 A	*** hardcoded ***
;	outdir	[INPUT; default: /data/fubar/SCAR/emissivity/emisscie]
;		directory to place outputs in
;	root	[INPUT; default: root] root prefix for output files
;	eqfile	[INPUT] pathname, relative to CHIDIR, of file containing
;		ionization equilibrium values
;		* default: ioneq/arnaud_rothenflug.ioneq
;		* passed w/o comment to CONT_CIE
;	abund	[INPUT] abundances (default is from Anders & Grevesse)
;		* passed w/o comment to CONT_CIE
;	_extra	[INPUT] pass defined keywords to subroutines
;		WMN, WMX, NBIN, CIEDIR, CHIDIR
;
;description
;	for each temperature, computes the emissivity using CIEDIR/cie
;	and places them all in the files
;		ROOT_wvl	- all the bin-beginning values and the final
;				  bin ending value [Ang]
;		ROOT_tmp	- the temperature values [log(T)]
;		ROOT_src	- description of parameters used
;		ROOT_##.#	- 2D array (tmp,wvl) of H+He emissivities
;				  I(tem,wvl)@(logP|log(n_e))
;		ROOTATOM_##.#	- As above, for H+He+ATOM
;	a note on the units:
;	in SPEX, the output is power [1e44 ph/s/keV], for an EM of 1e64 /m^3
;	in CIE, the output is "emissivity" [ph/m^3/s/keV] (later converted
;	  to [ph/m^3/s/Ang]).  the question is, where does the "/m^3" come
;	  from in CIE, and how do we get to [ergs cm^3/s/A]?
;	answer: there is an extra multiplication by nH^2 [(cm^-3)^2] and the
;	  EM is 1e50 /cm^3, and the constant in front is 3e-9, not 3e-1.
;	CIE/nH^2 :: [ph/m^3/s/keV]/[cm^-6]->[ph/m^3/s/keV]*[cm^6]
;		->[1e-12 ph m^3/s/keV]
;	correcting for the constant in front :: [1e-20 ph m^3/s/keV]
;	EM*CORR*CIE/nH^2 :: [1e-20 ph m^3/s/keV]*[1e44 m^-3]
;		->[1e24 ph/s/keV]
;	correct for SPEX's default EM :: [1e24 ph/s/keV]*[1e64/1e44]
;		->[1e44 ph/s/keV] -> voila!
;
;	so, take CIE, divide by 1e6 [(cm/m)^3], divide by nH^2 [(cm^-3)^2],
;	multiply by energy of photon at this wavelength,
;	multiply by [1e23] (to make numbers large),
;	to end up with [1e-23 ergs cm^3/s/(whatever)]
;
;restrictions
;	* requires CIE to have been compiled and accessible
;	* requires ln -s CIEDIR/*.dat $cwd/.
;	* works on UNIX.  only.
;	* requires subroutines
;	  CONT_CIE, NENH, RD_IONEQ, GETABUND, SYMB2ZION, LAT2ARAB
;		
;history
;	vinay kashyap (Dec97)
;	changed default OUTDIR to point to fubar (VK; Nov98)
;-

;initialize
allZ=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',$
	'Ni','Cu','Zn']				;elements from 1-30

;	usage
nelem=n_elements(elem)
if nelem lt 1 then begin
  print,'Usage: wrt_cont_cie,elem,pres,logT,wvl,reH,n_e=n_e,wrange=wrange,$'
  print,'       outdir=outdir,root=root,wmn=wmn,wmx=wmx,nbin=nbin,abund=abund,$'
  print,'       ciedir=ciedir,eqfile=eqfile,chidir=chidir'
  print,'  writes out continuum emissivities using CIE'
  return
endif

;	check output names
if not keyword_set(outdir) then outdir='/data/fubar/SCAR/emissivity/emisscie'
if not keyword_set(root) then root='root'
dbdir=outdir & if keyword_set(n_e) then dbdir=outdir+'D'
;	suffix for emissivity tables
if n_elements(pres) eq 0 then sufx=15. else sufx=alog10(pres)
if keyword_set(n_e) then sufx=alog10(n_e)
if sufx ge 0 and sufx lt 10 then sufx='0'+string(sufx,'(f3.1)')
if sufx lt 0 or sufx ge 10 then sufx=string(sufx,'(f4.1)')

;	if WRANGE is set, figure out wavelength range for CONT_CIE
if keyword_set(wrange) then begin
  nwr=n_elements(wrange)
  if nwr eq 2 then begin
    wmn=(wrange(0) < (wrange(1))) > (0.1)
    wmx=(wrange(1) > (wrange(0))) > (0.1)
  endif else begin
    wmn=1. & wmx=180.
  endelse
endif

;	first get the emissivity for H+He, to serve as a baseline
em0=cont_cie('He',pres,logT,wvl,reH,wmn=wmn,wmx=wmx,n_e=n_e,$
	eqfile=eqfile,abund=abund, _extra=e)
wfil=dbdir+'/'+root+'_wvl'	;wavelengths [Ang]
tfil=dbdir+'/'+root+'_tmp'	;log(temperatures [K])
sfil=dbdir+'/'+root+'_src'	;parameter description
efil=dbdir+'/'+root+'_'+sufx	;emissivities
nw=n_elements(wvl) & nt=n_elements(logT)
openw,u,wfil,/get_lun & writeu,u,nw,wvl & close,u & free_lun,u
openw,u,tfil,/get_lun & writeu,u,nT,logT & close,u & free_lun,u
openw,u,efil,/get_lun & writeu,u,nT,nw,em0 & close,u & free_lun,u
openw,u,sfil,/get_lun
  printf,u,'CIE'		;where this comes from
  printf,u,dbdir		;output directory
  if keyword_set(eqfile) then printf,u,eqfile else $
    printf,u,'ioneq/Arnaud_Rothenflug.ioneq'	;ion balance
  if keyword_set(abund) then printf,u,abund else $
    printf,u,getabund('anders & grevesse')	;adopted abundances
close,u & free_lun,u

;	next get the emissivities for all the other elements, and write to file
for iz=0,nelem-1 do begin			;{for each element
  symb2zion,elem(iz),zz,ii & atm=allZ(zz-1)
  message,'working on '+atm,/info
  emis=cont_cie(['He',atm],pres,logT,wvl,reH,wmn=wmn,wmx=wmx,n_e=n_e,$
	eqfile=eqfile,abund=abund, _extra=e)
  emis=emis-em0			;leaving just the contribution from ELEM(IZ)
  if total(emis) gt 1e-6 then begin	;(write iff something to write!
    efil=dbdir+'/'+root+atm+'_'+sufx
    openw,u,efil,/get_lun & writeu,u,nT,nw,emis & close,u & free_lun,u
  endif					;total(emis)>1e-6)
endfor						;IZ=0,NELEM-1}

return
end
