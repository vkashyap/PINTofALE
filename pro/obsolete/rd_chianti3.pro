function rd_chianti3,z,ion,delec=delec,chidir=chidir, _extra=e
;+
;function	rd_chianti3
;	one stop function to read in all relevant atomic parameters from
;	the CHIANTI database and return them in a comprehensive structure
;	(see structure variable CHELP for complete description).
;
;usage
;	line_database=rd_chianti3(Z,ion,delec=delec,chidir=chidir)
;
;parameters
;	Z	[INPUT; required] atomic number
;	ION	[INPUT] ionic state [if 0 or not given, assumes
;		fully ionized, i.e., =Z]
;		* Z and ION must be integer scalars
;
;keywords
;	delec	[INPUT] flag to specifies dielectronic recombination (1/0)
;	chidir	[INPUT] path to top directory of CHIANTI installation
;		[default: /data/fubar/SCAR/CHIANTI/]
;	_extra	[INPUT] junk -- ignore
;
;restrictions
;	* requires CHIANTI or something very similar
;	* requires subroutines READ_WGFA2, READ_ELVLC, READ_SPLUPS, which
;	  are part of CHIANTI
;
;history
;	written by vinay kashyap (Nov96), based on SYNTHETIC.PRO of CHIANTI
;	added keyword DELEC, modified default behavior of ION (VK; OctMM)
;	changed name to RD_CHIANTI3 because this won't work with v4
;	  anymore (VK; Jun02)
;-

message,'OBSOLETE: use RD_CHIANTI instead',/informational

chelp={CHIANTI_HELP,$
	NTRANS: '# transitions (LONG)',$
	NLEV: 'total # levels (LONG)',$
	LEV1: 'lower level indices (INTARR(NTRANS))',$
	LEV2: 'upper level indices (INTARR(NTRANS))',$
	WVL: 'transition wavelengths [Ang] (FLTARR(NTRANS))',$
	GF: 'weighted oscillator strengths (FLTARR(NTRANS))',$
	A: 'radiative transition probability [s^-1] (FLTARR(NTRANS))',$
	L1: 'level index (INTARR(NLEV))',$
	TERM: 'configuration index (STRARR(80,NLEV))',$
	CONF: 'configuration description (INTARR(NLEV))',$
	SS: '2S+1 (INTARR(NLEV))',$
	L: 'L (INTARR(NLEV))',$
	J: 'J (FLTARR(NLEV))',$
	ECMOBS: 'Observed Energy of level [cm^-1] (FLTARR(NLEV))',$
	ECMTHR: 'Theoretical Energy of level [cm^-1] (FLTARR(NLEV))',$
	ERYOBS: 'Observed Energy of level [Ryd] (FLTARR(NLEV))',$
	ERYTHR: 'Theoretical Energy of level [Ryd] (FLTARR(NLEV))',$
	TRANS: 'Transition type {cf.B&T,A&A1992,254,436} (FLTARR(NLEV,NLEV))',$
	GFU: 'weighted oscillator strengths (FLTARR(NLEV,NLEV))',$
	DERYD: 'energy difference between levels [Ryd] (FLTARR(NLEV,NLEV))',$
	C_UPS: 'Burgess & Tully scaling parameter (FLTARR(NLEV,NLEV))',$
	SPUPS: 'spline fits to scaled upsilons (FLTARR(NLEV,NLEV,5))',$
	WREF: 'List of references from .wgfa',$
	EREF: 'List of references from .elvlc',$
	UREF: 'List of references from .splups' }

;	usage
ok='ok' & np=n_params(0) & nZ=n_elements(Z) & nI=n_elements(ion)
sZ=size(Z) & nsZ=n_elements(sZ) & sI=size(ion) & nsI=n_elements(sI)
if np eq 0 then ok='Insufficient parameters' else $
 if nZ gt 1 then ok='Z cannot be an array' else $
  if nI gt 1 then ok='ION cannot be an array' else $
   if nZ ne nI then ok='Z and ION incompatible' else $
    if sZ[nsZ-2] gt 5 then ok='Z must be a number' else $
     if sI[nsI-2] gt 5 then ok='ION must be a number'
if ok ne 'ok' then begin
  print,'Usage: line_database = rd_chianti3(Z,ION_STATE,chidir=chidir)'
  print,'  reads CHIANTI database for given element and ionic state and'
  print,'  returns the info in a structure'
  if np ne 0 then message,ok,/info
  return,chelp
endif
;
if not keyword_set(ion) then ion=Z[0]		;list of ions to consider

;where is the CHIANTI database?
if not keyword_set(chidir) then chidir='/data/fubar/SCAR/CHIANTI/'

;initialize
atom = ['h','he','li','be','b','c','n','o','f','ne','na','mg','al','si',$
	'p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe','co','ni',$
	'cu','zn']

;	these are the input files
root=atom[z[0]-1]+'_'+strtrim(string(ion[0],'(i2)'),2)
if keyword_set(delec) then root=root+'d'
wfil=chidir+'/'+atom[z[0]-1]+'/'+root+'/'+root+'.wgfa'
efil=chidir+'/'+atom[z[0]-1]+'/'+root+'/'+root+'.elvlc'
ufil=chidir+'/'+atom[z[0]-1]+'/'+root+'/'+root+'.splups'

;	read wavelengths, oscillator strengths, and transition probs...
read_wgfa2,wfil,lev1,lev2,wvl,gf,a,wref
ntrans=n_elements(lev1) & nlev=max([lev1,lev2])	;how many lines/levels

;	read level configurations, angular momenta, transition energies...
read_elvlc,efil,l1,term,conf,ss,ll,jj,ecm,ery,ecmt,eryt,eref
;	if no observations available, fall back on theory
oo=where(ecm eq 0)
if max(oo) gt 0 then begin
  ecm[oo]=ecmt[oo] & ery[oo]=eryt[oo]
endif

;	read transition type, collissional osc.strengths, coeffs, upsilons
read_splups,ufil,trtype,gfu,deu,c_ups,spups,uref

;	place into structure for output

chline={NTRANS: ntrans, NLEV: nlev,$
	LEV1: lev1, LEV2: lev2, WVL: wvl, GF: gf, A: a,$
	L1: l1, TERM: term, CONF: conf, SS: ss, L: ll, J: jj,$
	ECMOBS: ecm, ECMTHR: ecmt, ERYOBS: ery, ERYTHR: eryt,$
	TRANS: trtype, GFU: gfu, DERYD: deu, C_UPS: c_ups, SPUPS: spups,$
	WREF: wref, EREF: eref, UREF: uref }

return,chline
end
