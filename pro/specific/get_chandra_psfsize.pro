function get_chandra_psfsize,theta,phi,energy=energy,eefrac=eefrac,$
	eeftbl=eeftbl,idx_eef=idx_eef,idx_nrg=idx_nrg,idx_phi=idx_phi,$
	verbose=verbose, _extra=e
;+
;function	get_chandra_psfsize
;	read in Diab's table that lists the PSF sizes and return
;	the radius [arcsec] that encloses a specified fraction
;	of the PSF
;
;syntax
;	psfradii=get_chandra_psfsize(theta,phi,energy=energy,eefrac=eefrac,$
;	eeftbl=eeftbl,idx_eef=idx_eef,idx_nrg=idx_nrg,idx_phi=idx_phi,$
;	verbose=verbose)
;
;parameters
;	theta	[INPUT; required] off-axis angles in [arcmin]
;	phi	[INPUT; optional] azimuthal angles in [degree]
;		* if phi is not given, assumed to be 0
;		* if N(PHI)=1, all PHI are taken to be PHI[0]
;		* if N(PHI).ne.N(THETA), only PHI[0] is used
;
;keywords
;	energy	[INPUT; default=1.5 keV] energy at which to determine
;		PSF size
;	eefrac	[INPUT; default=0.9] fraction of PSF enclosed
;		* note that in case of both ENERGY and EEFRAC,
;		  there is no interpolation -- entries closest
;		  to specified value will be chosen.  (users can
;		  always carry out the interpolations post facto)
;		* if <0, ABS(EEFRAC) is assumed
;		  if >1 and <100, assumed to be percentage
;		  if >100, then 1-1/EEFRAC is used
;	eeftbl	[INPUT; '/data/L3/fap/ecf/hrmaD1996-12-20hrci_ecf_N0002.fits.gz']
;		full path name to Diab Jerius' compilation of
;		enclosed energy fractions
;	idx_eef	[OUTPUT] index of EE fraction column used
;	idx_nrg	[OUTPUT] index of energy column used
;	idx_phi	[OUTPUT] index of PHI column used
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Apr2006)
;-

;	usage
ok='ok' & np=n_params() & nt=n_elements(theta)
if np eq 0 then ok='Insufficient parameters' else $
 if nt eq 0 then ok='THETA is undefined'
if ok ne 'ok' then begin
  print,'Usage: psfradii=get_chandra_psfsize(theta,phi,energy=energy,$'
  print,'       eefrac=eefrac,eeftbl=eeftbl,idx_eef=idx_eef,idx_nrg=idx_nrg,$'
  print,'       idx_phi=idx_phi,verbose=verbose)'
  print,'  read in table of enclosed energy fractions and return the
  print,'  size of the PSF'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
frac=0.9
if keyword_set(eefrac) then begin
  frac=abs(eefrac[0])
  if frac ge 1 and frac lt 100 then frac=frac/100.
  if frac ge 100 then frac=1.D - 1./frac
endif
;
nrg=1.5 & if keyword_set(energy) then nrg=float(energy[0])
;
eeftab='/data/L3/fap/ecf/hrmaD1996-12-20hrci_ecf_N0002.fits.gz'
;	note -- "ECF" seems to mean "EEF" here
if keyword_set(eeftbl) then begin
  sze=size(eeftbl) & nsze=n_elements(sze)
  if sze[nsze-2L] eq 7 then ecftab=eeftbl[0] else message,$
	'EEFTBL cannot be understood; using default',/informational
endif
;
th=[theta[*]] & ph=0*th
nph=n_elements(phi)
if nph eq nt then ph=[phi[*]] else begin
  if nph gt 0 then ph[*]=phi[0]
endelse

;	read in EEFTBL
if vv gt 10 then message,'Reading from '+eeftab,/informational
eef=mrdfits(eeftab,1,heef)
;	obtain RADIUS[ECF,THETA,PHI,ENERGY]
meef=n_elements(eef.ECF)
mth=n_elements(eef.THETA)
mph=n_elements(eef.PHI)
mnrg=n_elements(eef.ENERGY)

;	locate the appropriate columns
rad=eef.RADIUS
tmp=min(abs(frac-eef.ECF),idx_eef)
rad=reform(rad[idx_eef,*,*,*])
tmp=min(abs(nrg-eef.ENERGY),idx_nrg)
rad=reform(rad[*,*,idx_nrg])

;	how many unique PHIs?
phu=ph[uniq(ph,sort(ph))] & nphu=n_elements(phu)

;	output
psfradii=fltarr(nt) & idx_phi=lonarr(nt)

;	for each unique PHI, interpolate in THETA
for ip=0L,nphu-1L do begin
  op=where(ph eq phu[ip],mop)
  tmp=min(abs(phu[ip]-eef.PHI),iph)
  if vv gt 5 then print,'PHI = ',(eef.PHI)[iph]
  if mop eq 0 then message,'BUG!'
  idx_phi[op]=iph
  psfradii[op]=interpol(reform(rad[*,iph]),eef.THETA,th[op])>0
endfor

return,psfradii
end
