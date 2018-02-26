pro apedance,line,Z,abund=abund,apabref=apabref,verbose=verbose, _extra=e
;+
;procedure	apedance
;	remove or correct the abundance dependance of APED emissivities
;
;	APED emissivities as stored on disk include both ion balance
;	(Mazzotta et al.) and abundances (Anders & Grevesse).  This
;	can be inconvenient in some PINTofALE tasks since other line
;	emissivity databases in PINTofALE do not include either factor.
;	It is not possible to remove the effects of the included ion
;	balance, but removing the abundance dependance is simply a
;	matter of dividing the emissivity of each line by the abundance
;	appropriate to the element producing that line.
;
;syntax
;	apedance,line,Z,abund=abund,apabref=apabref,verbose=verbose
;
;parameters
;	line	[I/O] APED line emissivities
;		* emissivity structure out of RD_LINE()
;		* if 2D array, assumed to be an array of size (T,Wvl)
;	Z	[INPUT] atomic numbers
;		* required if LINE is an array and not a structure
;		* must match size of 2nd dimension in LINE
;
;keywords
;	abund	[INPUT] if given, updates LINE by multiplying by
;		the ratio of ABUND/(Anders & Grevesse)
;		* the default is to just divide by Anders & Grevesse
;		  abundances, effectively removing the abundance
;		  dependance from the emissivities
;	apabref	[INPUT] if set to a reference other than Anders & Grevesse,
;		it is assumed that the input LINE emissivities need to be
;		corrected according to the specified abundances.
;		* WARNING: do not set this unless you know exactly
;		 what you are doing
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to avoid crashing the program
;
;warning
;	no checks are made to verify that the emissivities have
;	have already had the abundances taken out, or that the
;	abundances are corrected according to some other baseline.
;	USE CAREFULLY!
;
;history
;	vinay kashyap (Apr2004)
;-

;	usage
ok='ok' & np=n_params() & nl=n_elements(line) & nZ=n_elements(Z)
nlt=n_tags(line) & szl=size(line)
if np eq 0 then ok='Insufficient parameters' else $
 if nl eq 0 then ok='LINE is undefined' else $
  if nlt ne 0 then begin	;(LINE is a structure
    tagnames=tag_names(line)
    i1=where(strupcase(tagnames) eq 'LINE_INT',mi1)
    i2=where(strupcase(tagnames) eq 'Z',mi2)
    if mi1 eq 0 or mi2 eq 0 then ok='structure in unknown format'
  endif else begin		;structure)(array
    nszl=n_elements(szl)
    if np eq 1 then ok='Z is required' else $
     if nZ eq 0 then ok='Z is undefined' else $
      if szl[0] ne 2 then ok='LINE must be a 2D array of size (T,Wvl)' else $
       if szl[2] ne nZ then ok='LINE and Z are incompatible'
  endelse			;array)
if ok ne 'ok' then begin
  print,'Usage: apedance,line,Z,abund=abund,apabref=apabref,verbose=verbose'
  print,'  remove or correct the abundance dependance of APED emissivities'
  if np gt 0 then message,ok,/informational
  return
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
if keyword_set(apabref) then begin
  if vv gt 0 then message,$
	'Assuming input emissivities have abundances from '+$
	apabref[0],/informational
  apdefab=getabund(apabref[0])
endif else begin
  if vv gt 0 then message,$
	'Assuming input emissivities have abundances from Anders & Grevesse',/informational
endelse
if n_elements(apdefab) lt 30 then begin
  apdefab=getabund('anders & grevesse')
endif
;
abnd=0*apdefab+1.
if n_elements(abund) ge 30 then begin
  if vv gt 0 then message,$
  	'resetting abundances of emissivities to input abundances',/informational
  abnd=abund
endif else begin
  if vv gt 0 then message,$
  	'stripping emissivities of abundances',/informational
endelse
;
if nlt ne 0 then begin
  lemis=line.LINE_INT
  lZ=long(line.Z)
endif else begin
  lemis=line
  lZ=long(Z)
endelse

;	correct for abundance
fZ=abnd[lZ-1L]/apdefab[lZ-1L]
szl=size(lemis)
for iT=0L,szl[1]-1L do lemis[iT,*]=lemis[iT,*]*fZ[*]

;	outputs
if nlt ne 0 then line.LINE_INT=lemis else line=lemis

return
end
