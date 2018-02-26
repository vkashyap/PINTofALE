function hydrogenic,Z,N,nmax=nmax,elimit=elimit,okeV=okeV,$
	lyman=lyman,balmer=balmer,paschen=paschen,brackett=brackett,$
	pfund=pfund,humphreys=humphreys, _extra=e
;+
;function	hydrogenic
;	return the theoretical line sequence for the hydrogenic
;	series for a given principal quantum number of the 1-electron
;	configuration of a given element.
;
;	warning: this routine does not take into account fine structure
;	corrections due to the velocity dependence of electron mass,
;	electron spin, and the Lamb shift due to radiation fields.
;
;syntax
;	Hwvl=hydrogenic(Z,N,nmax=nmax,elimit=elimit,/okeV,$
;	/lyman,/balmer,/paschen,/brackett,/pfund,/humphreys)
;
;parameters
;	Z	[INPUT; required] the atomic number of the element
;		for which the hydrogenic sequence is to be calculated
;		* may also be a symbol
;	N	[INPUT] the principal quantum number that defines the
;		series.  if not given, assumed to be 1
;		* if given, overrides the keywords below
;
;keywords
;	nmax	[INPUT] number of lines to include in the output
;		* default is NMAX=10
;	elimit	[OUTPUT] the series limit, in the same units as
;		the primary output
;	okeV	[INPUT] if set, returns the line energies in [keV]
;		* default is to return [Ang]
;	lyman	[INPUT] if set and N is not given, assumes N=1
;	balmer	[INPUT] if set and N is not given, assumes N=2
;	paschen	[INPUT] if set and N is not given, assumes N=3
;	brackett	[INPUT] if set and N is not given, assumes N=4
;	pfund	[INPUT] if set and N is not given, assumes N=5
;	humphreys	[INPUT] if set and N is not given, assumes N=6
;		* if more than one of the above keywords are specified, then
;		  the one that corresponds to the smallest N is adopted.
;
;subroutines
;	INICON
;
;history
;	vinay kashyap (Jun02)
;-

;	usage
ok='ok' & np=n_params() & nZ=n_elements(Z)
if np eq 0 then ok='Insufficient parameters' else $
 if nZ eq 0 then ok='Z: not defined' else $
  if nZ gt 1 then ok='Z: must be a scalar'
if ok ne 'ok' then begin
  print,'Usage: Hwvl=hydrogenic(Z,N,nmax=nmax,elimit=elimit,/okeV,$'
  print,'       /lyman,/balmer,/paschen,/brackett,/pfund,/humphreys)'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	what is the atomic number?
iZ=Z[0] & sZ=size(Z) & nsZ=n_elements(sZ)
if sZ[nsZ-2] eq 7 then symb2zion,Z[0],iZ,junk

;	what is the principal quantum number?
Nq=1
if n_elements(N) gt 0 then Nq=N[0] else $
 if keyword_set(lyman) then Nq=1 else $
  if keyword_set(balmer) then Nq=2 else $
   if keyword_set(paschen) then Nq=3 else $
    if keyword_set(brackett) then Nq=4 else $
     if keyword_set(pfund) then Nq=5 else $
      if keyword_set(humphreys) then Nq=6

;	how many lines in the output?
maxN=10 & if keyword_set(Nmax) then maxN=long(Nmax[0]) > 1

;	get the physical constants
inicon,fundae=fundae,amu=amu
;
hbar=fundae.H/2./!pi			;[ergs/s]
c=fundae.C				;[cm/s]
esu=fundae.ESU				;electron charge [esu]
eV=1e7*fundae.E				;[erg/eV]
mass_e=fundae.ME/fundae.AMU		;electron mass in [amu]
mass_n=(amu.(iZ-1))[0]			;mass of nucleus in [amu]
mu=mass_e*mass_n/(mass_e+mass_n)	;reduced mass [amu]
mu_gm=mu*fundae.AMU

;	so what is the Rydberg constant?
Ryd=10.^(alog10(mu_gm)+4.*alog10(esu)-2.*alog10(hbar))/2.	;[ergs]
Ryd_eV=Ryd/eV				;[eV]

;	get the energy levels
elimit=iZ^2*(Ryd_eV/1e3)/Nq^2		;[keV] by default
eserie=iZ^2*(Ryd_eV/1e3)/(findgen(maxN)+Nq+1)^2	;[keV]
elines=elimit-eserie				;[keV]

;	convert to [Ang] if necessary
if not keyword_set(okeV) then begin
  elimit=(fundae.KEVANG)/elimit
  elines=(fundae.KEVANG)/elines
endif

return,elines
end
