function cont_cie,elem,pres,logT,wvl,reH,wmn=wmn,wmx=wmx,nbin=nbin,$
	n_e=n_e,abund=abund,ciedir=ciedir,eqfile=eqfile,$
	cieinp=cieinp, ciespec=ciespec, _extra=e
;+
;FUNCTION	cont_cie
;	returns continuum emissivities (NLOGT,NBIN) using CIE (subset
;	of SPEX) [1e-23 ergs cm^3/s/A]
;
;	WARNING: The continuum emissivities returned here differ from the
;	PoA philosophy of line emissivities in that abundances and ion
;	balance are *included*!
;
;SYNTAX
;	emis=cont_cie(elem,pres,logT,wvl,reH,wmn=wmn,wmx=wmx,nbin=nbin,$
;	n_e=n_e,abund=abund,ciedir=ciedir,eqfile=eqfile,$
;	cieinp=cieinp,ciespec=ciespec,chidir=chidir)
;
;PARAMETERS
;	elem	[INPUT; required] name of element e.g., "He", "C", etc.
;		* may specify ionic state (e.g., 'Fe13'), but will be ignored
;		* only the elements that CIE can handle will be included
;		  in the final calculation.  As of Dec97, these were
;		  	C,N,O,Ne,Na,Mg,Al,Si,S,Ar,Ca,Fe,Ni
;		  in addition to H and He, which are ALWAYS included
;		* if invalid, ALL of the above will be included
;	pres	[INPUT; default: 1e15 cm^-3 K] electron pressure at which
;		to compute intensities
;	logT	[INPUT/OUTPUT] log(Temperature[K]) at which to compute the
;		intensities.  if not given, then LOGT=4.0+FINDGEN(81)*0.05
;		* best results if evenly spaced in log(T)
;	wvl	[OUTPUT] all the bin-beginning values and the final
;		bin ending value in the spectrum [Ang]
;	reH	[OUTPUT] n_e/n_H for each LOGT
;
;KEYWORDS
;	wmn	[INPUT; default=1] minimum in wavelength range [Ang]
;		* absolute minimum is 0.1 ***hardcoded***
;	wmx	[INPUT; default=180] maximum in wavelenth range [Ang]
;	nbin	[INPUT; default=1000] number of bins in the spectrum
;		* if -ve, then binning is logarithmic
;	n_e	[INPUT] electron density [/cm^3]
;		* OVERRIDES values determined using PRES and LOGT
;	abund	[INPUT] abundances (default is from Anders & Grevesse)
;	ciedir	[INPUT] directory containing the CIE executable
;		* default = /data/fubar/SCAR/CIE
;	eqfile	[INPUT] pathname, relative to CHIDIR, of file containing
;		ionization equilibrium values
;		* default: ioneq/arnaud_rothenflug.ioneq
;	cieinp	[INPUT; default='/tmp/cie_inp'] command file name
;	ciespec	[INPUT; default='tmp.spec'] file containing output spectrum
;
;	_extra	[INPUT ONLY] use to pass defined keywords to NENH
;		* CHIDIR
;
;description
;	for each temperature, computes the emissivity using CIEDIR/cie
;	and returns a density "insensitive" value
;
;a note on the units:
;	in SPEX, the output is power [1e44 ph/s/keV], for an EM of 1e64 /m^3
;	in CIE, the output is "emissivity" [ph/m^3/s/keV] (later converted
;	  to [ph/m^3/s/Ang]).  the question is, where does the "/m^3" come
;	  from in CIE, and how do we get to [ergs cm^3/s/A]?
;	answer: there is an extra multiplication by nH^2 [(cm^-3)^2] and the
;	  EM is 1e50 /cm^3, and the constant in front is 3e-9, not 3e-1.
;	CIE/(n_e*n_H) :: [ph/m^3/s/keV]/[cm^-6]->[ph/m^3/s/keV]*[cm^6]
;		->[1e-12 ph m^3/s/keV]
;	correcting for the constant in front :: [1e-20 ph m^3/s/keV]
;	EM*CORR*CIE/(n_e*n_H) :: [1e-20 ph m^3/s/keV]*[1e44 m^-3]
;		->[1e24 ph/s/keV]
;	correct for SPEX's default EM :: [1e24 ph/s/keV]*[1e64/1e44]
;		->[1e44 ph/s/keV] -> voila!
;
;	so, take CIE, divide by 1e6 [(cm/m)^3], divide by n_e*n_H [(cm^-3)^2],
;	multiply by energy of photon at this wavelength,
;	multiply by [1e23] (to make numbers large),
;	to end up with [1e-23 ergs cm^3/s/(whatever)]
;
;restrictions
;	* requires CIE to have been compiled and accessible
;	* requires ln -s CIEDIR/*.dat $cwd/.
;	* works on UNIX.  only.
;	* uses CHIANTI compilation of ion balance
;	* requires subroutines
;	  -- GETABUND [SETABUND]
;	  -- SYMB2ZION [LAT2ARAB]
;	  -- NENH [RD_IONEQ [READ_IONEQ (a CHIANTI routine)]
;	  -- KILROY
;
;history
;	vinay kashyap (Dec97)
;	minimum abundance now always > 0 (due to CIE bug); added keywords
;	  CIEINP and CIESPEC (VK; Jan98)
;-

;initialize
allZ=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',$
	'Ni','Cu','Zn']				;elements from 1-30
cieZ=[	'He','C','N','O','Ne','Na','Mg','Al','Si','S','Ar','Ca','Fe','Ni']
ncZ=n_elements(cieZ) & cZ=intarr(ncZ)
for i=0,ncZ-1 do begin
  symb2zion,cieZ(i),zz,ii & cZ(i)=zz
endfor

;usage
np=n_params(0)
if np eq 0 then begin
  print,'Usage: emis=cont_cie(elem,pressure,logT,wvl,ne_to_nH,wmn=wmn,wmx=wmx,$'
  print,'       nbin=nbin,n_e=n_e,abund=abund,ciedir=ciedir,eqfile=eqfile,$'
  print,'       cieinp=cieinp,ciespec=ciespec,chidir=chidir)'
  print,'  returns 2D array of continuum emissivities for given element(s) at'
  print,'  specified wavelengths and given temperatures'
  return,-1L
endif

;figure out which elements
nelem=n_elements(elem) & elm=['']
for i=0,nelem-1 do begin
  symb2zion,elem(i),zz,ii
  if zz ne 0 then elm=[elm,allZ(zz-1)]
endfor
if n_elements(elm) eq 1 then elm=cieZ else elm=elm(1:*)
elm=['H','He',elm]					;gotta keep H and He
;
;choose only subset of ELM that is relevant
nelem=n_elements(elm) & elmZ=[''] & Z=[0]
for i=0,nelem-1 do begin
  oo=where(elm(i) eq cieZ,moo)
  if moo gt 0 then begin
    elmZ=[elmZ,cieZ(oo(0))]
    symb2zion,elm(i),zz,ii & Z=[Z,zz]
  endif
endfor
elmZ=elmZ(1:*) & Z=Z(1:*)
oZ=uniq(Z,sort(Z)) & elmZ=elmZ(oz) & Z=Z(oz) & nelem=n_elements(elmZ)

;pressure
if n_elements(pres) eq 0 then pr=1e15 else pr=pres

;temperature array
nt=n_elements(logT)
if nt eq 1 then t=[float(logT)]
if nt lt 1 then begin			;undefined
  nt=81 & t=4.0+findgen(nt)*0.05 & logT=t
endif else t=float(logT)
TkeV=[t(*)] & nt=n_elements(TkeV)
TkeV=10.^(TkeV)*1.380662e-16/1.6021892e-9	;[keV]=[K]*[erg/K]/[erg/keV]

;where is CIE?
if not keyword_set(ciedir) then ciedir='/data/fubar/SCAR/CIE'

;which ion balance to use?
if not keyword_set(eqfile) then eqfile='ioneq/arnaud_rothenflug.ioneq'

;electron density
if keyword_set(n_e) then begin		;(density is specified
  pr = n_e(0)*10.^(t)		;override pressure
  eden = 0*t + n_e(0)
endif else eden = pr/10.^(t)		;n_e)

;abundances
defabu=getabund('anders & grevesse')
if n_elements(abund) lt 30 then abund=defabu
;set all unspecified element abundances to 0
abnd=abund & abnd(Z-1)=-abnd(Z-1) & abnd=-abnd
;	abnd=(abnd>0)
oa=where(abnd lt 0,moa) & if moa gt 0 then abnd(oa)=-(abnd(oa)/1e8) > (1e-15)
abnd(0)=1.

;get H density
reH=nenh(t,abund=abnd,eqfile=eqfile,elem=['H',elmZ],_extra=e) & nH=eden/reH

;wavelength range definitions.
if not keyword_set(wmn) then wmn=1. & wmn=(wmn>0.1)
if not keyword_set(wmx) then wmx=180. & wmx=(wmx>0.1)
if wmx eq wmn then begin
  message,'wavelength range incorrect',/info & return,0*TkeV
endif
if wmx lt wmn then begin
  wmn=wmn+wmx & wmx=wmn-wmx & wmn=wmn-wmx	;exchange
endif
;
if n_elements(nbin) ne 1 then nbin=1000L
if abs(nbin(0)) eq 1 then nbin=1000L & nbn=long(abs(nbin))
if nbin(0) lt 0 then logbin=1 else logbin=0

;	initialize
h=6.626176e-27	;[ergs s]
c=2.99792458e10	;[cm s^-1]
hc=h*c		;[ergs cm]
comfil='/tmp/cie_inp'		;command file
if keyword_set(cieinp) then comfil=string(cieinp)
ciefil='cie'			;CIE executable
dmpfil='tmp.spec'		;temporary ascii dump of spectrum
if keyword_set(ciespec) then dmpfil=string(ciespec)
rehsrch='electron/Hydrogen density:'	;ratio of n_e/n_H follows after this
lrehsrch=strlen(rehsrch)	;length of above string
ionbal=1		;1=A&R, 2=M&G, 3=Jacobs et al.
conly=2			;1=continuum+lines, 2=continuum, 3=lines
getemis=2		;1=[ph/m^2/s/A], 2=[ph/m^3/s/A]
broad=0			;0=none, 1=thermal broadening
chabun='y'		;y=change abundances, n=don't
bbflx='n'		;y=include diluted blackbody flux, n=no
dfact=0.1		;dilution factor
factor=fltarr(ncZ) & factor=abnd(cZ-1)/defabu(cZ-1)	;ratio of abundances

;	generate emissivity spectrum
emis=dblarr(nt,nbn) & rreH=fltarr(nt)
for it=0,nt-1 do begin			;{call CIE for each temperature
  kilroy; was here
  openw,ucom,comfil,/get_lun		;(create command file for CIE
    printf,ucom,strtrim(logbin,2)	;Wavelength grid: 0=lin, 1=log
    printf,ucom,strtrim(wmn,2),' ',strtrim(wmx,2)
					;start & stop wvls [Ang]
    printf,ucom,strtrim(nbn,2)		;number of steps
    printf,ucom,strtrim(ionbal,2)	;ion balance
    printf,ucom,strtrim(conly,2)	;continuum only
    printf,ucom,strtrim(getemis,2)	;spectrum or emissivity
    printf,ucom,strtrim(broad,2)	;broadening
    printf,ucom,strtrim(nH(it),2)	;H density
    printf,ucom,strtrim(TkeV(it),2)	;e temperature [keV]
    printf,ucom,strtrim(chabun,2)	;change abundances
    for iz=0,ncZ-1 do printf,ucom,strtrim(factor(iz),2)
					;ABUND/(Anders&Grevesse)
    printf,ucom,bbflx			;include diluted blackbody flux?
    if bbflx eq 'y' then printf,ucom,strtrim(dfact,2)	;dilution factor
    printf,ucom,dmpfil			;output file
  close,ucom,/all & free_lun,ucom	;close command file)

  spawn,'if ( -f '+dmpfil+' ) /bin/rm '+dmpfil
  spawn,ciedir+'/'+ciefil+' < '+comfil,scrdmp

  for i=0,n_elements(scrdmp)-1 do begin
    c1=scrdmp(i) & ireh=strpos(c1,rehsrch,0)
    if ireh ge 0 then begin
      rat_neh=strmid(c1,ireh+lrehsrch+1,14)
      rreh(it)=float(rat_neh)
    endif
  endfor
  print,it,logT(it),TkeV(it),rreh(it),reh(it)
  if rreh(it)/reh(it) gt 2 or rreh(it)/reh(it) lt 0.5 then $
	message,'n_e/n_H['+strtrim(logT(it),2)+']= '+strtrim(rreh(it),2)+$
	' (CIE) & '+strtrim(reh(it),2)+' (SCAR)',/info

  var=dblarr(2,nbn)
  openr,udmp,dmpfil,/get_lun		;(read in ASCII dump
    readf,udmp,var
    wvl=reform(var(0,*))		;wavelengths
    ems=reform(var(1,*))		;emissivities [ph/m^3/s/A]
  close,udmp,/all & free_lun,udmp	;close dump file)
  ems=ems*1e-6			;[../m^3..]->[../cm^3..]
  ems=ems/(rreh(it)*nH(it)^2)	;[ph/cm^3/..]->[ph cm^3/..]
	;NOTE -- doing this instead of straight division by eden(it)*nH(it)
	;because CIE says reH at low T is 1 (even when H is not ionized!)
	;which doesn't square with output of NENH.
  nrg=hc*1e8/wvl		;[ergs/ph]=[ergs cm]*[Ang/cm]/[Ang]
  ems=ems*nrg*1e23		;[ph..]->[1e-23 ergs..]
  emis(it,*)=ems		;[1e-23 ergs cm^3/s/A]

  ;	now get the bin boundaries
  ;(wvl=12.3985/(0.5*(E0+E1)), so convert it to keV first, then
  ;undo the averaging and then convert back to Ang)
  ;	NOTE: this introduces spurious wiggles in dW because of
  ;	machine inaccuracy.  better version when figured out.
  ww=[wvl,wmx] & ww=12.3985/ww
  for iw=nbn-1,0,-1 do ww(iw)=2*(ww(iw)-0.5*ww(iw+1))
  ww=12.3985/ww & ww(0)=wmn

  wvl=ww	;output

endfor					;IT=0,NT-1}

return,emis
end
