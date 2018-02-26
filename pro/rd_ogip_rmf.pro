function rd_ogip_rmf,rmfile,effar=effar,eqvar=eqvar,spcar=spcar,$
	rmfimg=rmfimg,fchcol=fchcol,shift1=shift1,verbose=verbose,$
	_extra=e
;+
;function	rd_ogip_rmf
;	reads in the response matrix and associated axes from OGIP-compliant
;	response file and returns everything in a structure of the form:
;
;	{NNRG, ELO, EHI, NCHAN, EMN, EMX, N_GRP, F_CHAN, N_CHAN, MATRIX, FIRSTCHAN}
;
;	where ELO and EHI refer to range of photon energies at which the
;	response is valid, EMN and EMX refer to bin boundaries for each
;	channel, N_GRP refers to number of groups of non-zero data, F_CHAN
;	refer to beginning indices of channels, N_CHAN are the number of
;	channels in each group, and the output response MATRIX excludes
;	zeros to save space.  The value of FIRSTCHAN indicates whether F_CHAN
;	indices are 0-based or 1-based.
;
;	http://legacy.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html
;	for a description of the format
;
;	WARNING: the EBOUNDS extension of RMFs, which defines EMN and EMX,
;	should never be used directly for scientific purposes.  They
;	contain the nominal energy boundaries for the PHA, PI, or pixel
;	bins that define the MATRIX.  They are never exact values, and
;	are useful only for plotting purposes.  It is the detector bins
;	and their indices that are the be-all and end-all of RMFs.
;
;syntax
;	rmstr=rd_ogip_rmf(rmfile,effar=effar,eqvar=eqvar,spcar=spcar,$
;	rmfimg=rmfimg,fchcol=fchcol,/shift1,verbose=verbose)
;
;parameters
;	rmfile	[INPUT; required] response matrix file
;
;keywords
;	effar	[OUTPUT] returns the effective area as a function of ELO
;		* if an RMF file is read instead of an RSP file, EFFAR are
;		  uniformly 1, unless it is ASCA in which case it is the QEs
;		* EFFAR are calculated ONLY if needed.  i.e., the keyword
;		  must be present in call (IDL5+) or defined on input (IDL5-)
;	eqvar	[OUTPUT] returns an equivalent effective area in each
;		detector channel, measuring the contribution to that
;		channel from all photon energies.
;	spcar	[OUTPUT] a variation on EQVAR, the convolution of a really
;		flat intensity spectrum [ph/.../keV] with the RSP.
;		In other words, what a flat spectrum would look like.
;		* BEWARE: EQVAR and SPCAR are _not_ true effective area.
;		  EFFAR are in _energy_ bins.  EQVAR and SPCAR are in
;		  _channel_ bins.
;		* CAVEAT VSVATOR!
;		* I will deny ever having said this, but if you _must_
;		  use one of these fake effective areas, use SPCAR.
;	rmfimg	[OUTPUT] returns the RMF matrix as an image of size
;		N(CHANNEL) x N(ENERGY)
;		* computed only if this keyword is present (must also
;		  be defined to be non-zero in IDL versions prior to 5.0)
;	shift1	[SEMI OBSOLETE] some RMFs have position indices that are
;		0-based and some are 1-based.  this keyword used to specify
;		which was what.  however, this information should be in the
;		FITS header in the keyword TLMIN{FCHCOL}, and that is where
;		the program takes its value from.
;		* BEWARE: if SHIFT1 is set _and_ TLMIN{FCHCOL} comes out as 0
;		  (usually if it is not defined, or perhaps it really is 0)
;		  then FIRSTCHAN is taken to be 1, and vice versa.
;		* in other words, if set, overrides whatever the file
;		  itself is saying
;	fchcol	[INPUT] column number of F_CHAN in RMFILE
;		* default is 4
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	requires IDLASTRO library routine MRDFITS
;
;history
;	vinay kashyap (Apr99; modified from RDRESP.PRO)
;	added keyword SHIFT1 (VK; JanMMI)
;	slight modification in defining IDL version (VK; FebMMI)
;	added keyword FCHCOL, overrode SHIFT1, added FIRSTCHAN to output
;	  structure (VK; Nov2001)
;	commented out the "L5_MATRIX" lines (VK; Dec2001)
;	added keywords EQVAR, SPCAR, and VERBOSE (VK; Jul02)
;	bug correction: rspfil variable changed to RMFILE (LL/VK;Aug04)
;       added ability to handle non-OGIP Astro-E2 XIS rmf file (LL/VK:Aug04)
;	tweaked to speed up calculation of effar (VK; Dec04)
;	added keyword RMFIMG (VK; May06)
;-

;	usage
ok='ok' & nrmf=n_elements(rmfile) & szr=size(rmfile) & nszr=n_elements(szr)
if nrmf eq 0 then ok='Name of Response Matrix File required' else $
 if nrmf gt 1 then ok='Cannot handle multiple RMFs!' else $
  if szr(nszr-2) ne 7 then ok='input must be filename'
if ok ne 'ok' then begin
  print,'Usage: rmstr=rd_ogip_rmf(rmfile,effar=effar,rmfimg=rmfimg,eqvar=eqvar,$'
  print,'       spcar=spcar,fchcol=fchcol,/shift1,verbose=verbose)'
  print,'  read OGIP-compliant RMF and return response matrix'
  if nrmf ne 0 then message,'	'+ok,/info
  return,{NNRG:0, ELO:0, EHI:0, NCHAN:0, EMN:0, EMX:0, N_GRP:0, F_CHAN:0,$
	N_CHAN:0, MATRIX:0, FIRSTCHAN:0}
endif

;	which column contains the first channel info?
fcol=4 & if n_elements(fchcol) ne 0 then fcol=fix(fchcol[0])
if fcol le 0 then begin
  message,'FCHCOL not understandable; ignoring and assuming =4',/info
  fcol=4
endif

;	verbosity
vv=0L & if keyword_set(verbose) then vv=long(verbose(0))>1

;	read from response matrix
if vv ge 5 then message,'reading extension 1 from file '+rmfile(0),/info
r1=mrdfits(rmfile(0),1,h1)
if vv ge 5 then message,'reading extension 2 from file '+rmfile(0),/info
r2=mrdfits(rmfile(0),2,h2)
x1=strtrim(sxpar(h1,'EXTNAME'),2) & x2=strtrim(sxpar(h2,'EXTNAME'),2)
x3=strtrim(sxpar(h2,'INSTRUME',2))
if x3 eq 'XIS' then begin
  message,'reading from non-OGIP compliant XIS RMF',/informational
  x2 = 'EBOUNDS'
endif

;	decode the structures r1 and r2
if x1 eq 'EBOUNDS' then begin 	;x2=SPECRESP MATRIX/MATRIX
  CHANNEL=r1.CHANNEL		;channel index
  EMN=r1.E_MIN			;channel bin minimum energy
  EMX=r1.E_MAX			;channel bin maximum energy
  ELO=r2.ENERG_LO		;minimum photon energy
  EHI=r2.ENERG_HI		;maximum photon energy
  N_GRP=r2.N_GRP		;number of groups of non-zero data
  F_CHAN=r2.F_CHAN		;beginning channel numbers of each group
  N_CHAN=r2.N_CHAN		;number of channels in each group
  MATRIX=r2.MATRIX		;non-zero elements of response
  ;L5_MATRIX=r2.L5_MATRIX	;<what?>
  FIRSTCHAN=sxpar(h2,'TLMIN'+strtrim(fcol,2))	;first channel index -- 0-based or 1-based
endif else if x2 eq 'EBOUNDS' then begin	;x1=SPECRESP MATRIX/MATRIX
    CHANNEL=r2.CHANNEL		;channel index
    EMN=r2.E_MIN		;channel bin minimum energy
    EMX=r2.E_MAX		;channel bin maximum energy
    ELO=r1.ENERG_LO		;minimum photon energy
    EHI=r1.ENERG_HI		;maximum photon energy
    N_GRP=r1.N_GRP		;number of groups of non-zero data
    F_CHAN=r1.F_CHAN		;beginning channel numbers of each group
    N_CHAN=r1.N_CHAN		;number of channels in each group
    MATRIX=r1.MATRIX		;non-zero elements of response
    ;L5_MATRIX=r1.L5_MATRIX	;number of non-zero elements of MATRIX in row
    FIRSTCHAN=sxpar(h1,'TLMIN'+strtrim(fcol,2))	;first channel index -- 0-based or 1-based
  endif else begin
    message,'	cannot understand input file: '+rmfil,/info
    return,0.
  endelse

if FIRSTCHAN eq 0 and keyword_set(shift1) then begin
  message,'FIRSTCHAN was 0, but SHIFT1 is forcing 1',/info
  FIRSTCHAN=1L
endif
if FIRSTCHAN eq 1 and n_elements(shift1) ne 0 then begin
  if shift1(0) eq 0 then begin
    message,'FIRSTCHAN was 1, but SHIFT1 is forcing 0',/info
    FIRSTCHAN=0L
  endif
endif

;	compatibility error checks
nchan=n_elements(channel) & nnrg=n_elements(elo) & szf=size(f_chan)
szn=size(n_chan) & szr=size(matrix)
ok=''		;assume everything is fine
if nchan ne n_elements(emn) then ok='CHANNEL x E_MIN' else $
 if nchan ne n_elements(emx) then ok='CHANNEL x E_MAX' else $
  if n_elements(emn) ne n_elements(emx) then ok='E_MIN x E_MAX'
if nnrg ne n_elements(ehi) then ok='ENERG_LO x ENERG_HI' else $
 if nnrg ne n_elements(n_grp) then ok='N_GRP x ENERG_LO' else $
  if n_elements(ehi) ne n_elements(n_grp) then ok='N_GRP x ENERG_HI' else $
   if nnrg ne szr(2) then ok='ENERG_LO x MATRIX'
if max(n_grp) gt 1 then begin
  if szf(1) ne max(n_grp) then ok='N_GRP x F_CHAN' else $
   if szf(2) ne nnrg then ok='ENERG_LO x F_CHAN' else $
    if szn(1) ne max(n_grp) then ok='N_GRP x N_CHAN' else $
     if szn(2) ne nnrg then ok='ENERG_LO x N_CHAN'
  tmp=intarr(nnrg) & for i=0,szn(1)-1 do tmp=tmp+n_chan(i,*)
  if max(tmp) ne szr(1) then ok='N_GRP x MATRIX'
endif
if ok ne '' then begin
  message,ok,/info & return,0.
endif

;	create output structure
rmstr=create_struct('NNRG',nnrg,'ELO',elo,'EHI',ehi,$
	'NCHAN',nchan,'EMN',emn,'EMX',emx,$
	'N_GRP',n_grp,'F_CHAN',f_chan,'N_CHAN',n_chan,$
	'MATRIX',matrix,'FIRSTCHAN',firstchan)

if keyword_set(shift1) then message,$
	'assuming that input RMF is 1-based, not 0-based',/info

;	figure out effective area
effar=total(matrix,1)
getea=0 & getimg=0 & idlver=float(strmid(!version.RELEASE,0,3))
if idlver gt 5 then begin
  if arg_present(eqvar) or arg_present(spcar) then getea=1
  if arg_present(rmfimg) then getimg=1
endif else begin
  if keyword_set(eqvar) or keyword_set(spcar) then getea=1
  if keyword_set(rmfimg) then getimg=1
endelse
if getimg eq 1 then rmfimg=fltarr(nchan,nnrg)
;
if getea eq 1 then begin
  if vv gt 5 then message,'extracting effective equivalent areas',/info
  ;effar=fltarr(nnrg)
  eqvar=fltarr(nchan) & spcar=fltarr(nchan)
  dbin=ehi-elo
  ;normspcar=fltarr(nchan) & numspcar=fltarr(nchan)
  cc='% completed: ' & bb='' & for i=1,5 do bb=bb+string(8b)
  print,cc+string(100.*0/float(nnrg),'(f5.1)'),bb,form='($,a18,a)'
  for inrg=0L,nnrg-1L do begin
    ngrp=n_grp(inrg) & jbeg=0
    for ig=0,ngrp-1 do begin
      if szn(0) gt 1 then begin
        ibeg=f_chan(ig,inrg) & iw=n_chan(ig,inrg)
      endif else begin
        ibeg=f_chan(inrg) & iw=n_chan(inrg)
      endelse
      if keyword_set(firstchan) then ibeg=ibeg-1	;IDL index correction
      if iw gt 0 then begin
	;effar(inrg)=effar(inrg)+total(matrix(jbeg:jbeg+iw-1,inrg))
	eqvar(ibeg:ibeg+iw-1)=eqvar(ibeg:ibeg+iw-1)+$
		matrix(jbeg:jbeg+iw-1,inrg)
	spcar(ibeg:ibeg+iw-1)=spcar(ibeg:ibeg+iw-1)+$
		matrix(jbeg:jbeg+iw-1,inrg)*dbin(inrg)
	;normspcar(ibeg:ibeg+iw-1)=normspcar(ibeg:ibeg+iw-1)+$
	;	dbin(inrg)
	;numspcar(ibeg:ibeg+iw-1)=numspcar(ibeg:ibeg+iw-1)+1
	if getimg gt 1 then rmfimg[ibeg:ibeg+iw-1,inrg]=matrix[jbeg:jbeg+iw-1,inrg]
      endif
      jbeg=jbeg+iw
    endfor
    if vv ge 10 and inrg eq 128L*long(inrg/128L) then print,$
	string(100.*inrg/float(nnrg),'(f5.1)'),bb,form='($,a5,a)'
  endfor
  ;ok=where(normspcar gt 0,mok)
  ;if mok gt 0 then spcar(ok)=spcar(ok)/normspcar(ok)
  dee=emx-emn & spcar=spcar/dee
  if vv ge 10 then print,'100.0'
  if vv gt 100 then stop,'HALTing.  type .CON to continue'
endif
;
if getea eq 0 and getimg eq 1 then begin
  for inrg=0L,nnrg-1L do begin
    ngrp=n_grp(inrg) & jbeg=0
    for ig=0,ngrp-1 do begin
      if szn(0) gt 1 then begin
        ibeg=f_chan(ig,inrg) & iw=n_chan(ig,inrg)
      endif else begin
        ibeg=f_chan(inrg) & iw=n_chan(inrg)
      endelse
      if keyword_set(firstchan) then ibeg=ibeg-1	;IDL index correction
      if iw gt 0 then rmfimg[ibeg:ibeg+iw-1,inrg]=matrix[jbeg:jbeg+iw-1,inrg]
      jbeg=jbeg+iw
    endfor
  endfor
endif

return,rmstr
end
