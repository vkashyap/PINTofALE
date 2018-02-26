function rdarf,arffil,ostr,units=units, _extra=e
;+
;function	rdarf
;	reads in the spectral response and associated axes from
;	OGIP-compliant ancillary response file
;
;syntax
;	arf=rdarf(arffil,ostr,units=units)
;
;parameters
;	arffil	[INPUT; required] name of ARF
;	ostr	[OUTPUT] structure containing the axes information for the
;		response: {NNRG, ELO, EHI}
;		* ELO and EHI refer to range of photon energies at which
;		  the response is valid
;
;keywords
;	units	[OUTPUT] contains the units of ARF; typically 'cm**2'
;	_extra	[INPUT] junk -- here only to prevent crashing the program
;
;restrictions
;	requires IDLASTRO library routine MRDFITS
;
;history
;	vinay kashyap (Jul97)
;	exits with 0 if ARFFIL not in fits format (VK; Apr05)
;-

;	usage
np=n_elements(arffil)
if np eq 0 then begin
  print,'Usage: arf=rdarf(arffil,ostr,units=units)'
  print,'  returns instrument spectral response curve'
  return,0.
endif

;	read from ancillary file
ar=mrdfits(arffil,1,har)
if n_tags(ar) eq 0 then begin
  if ar[0] eq 0 then begin
    message,arffil+' not in OGIP-compatible fits format',/information
    return,0.
  endif
endif
xnm=strtrim(sxpar(har,'EXTNAME'),2)

;	decode the structure
if xnm eq 'SPECRESP' then begin
  ELO=ar.ENERG_LO		;minimum photon energy
  EHI=ar.ENERG_HI		;maximum photon energy
  ARF=ar.SPECRESP		;spectral response
  coltyp=strtrim(sxpar(har,'TTYPE3'),2)
  if coltyp eq 'SPECRESP' then units=sxpar(har,'TUNIT3') else begin
    coltyp=strtrim(sxpar(har,'TTYPE2'),2)
    if coltyp eq 'SPECRESP' then units=sxpar(har,'TUNIT2') else begin
      coltyp=strtrim(sxpar(har,'TTYPE1'),2)
      if coltyp eq 'SPECRESP' then units=sxpar(har,'TUNIT1') else $
	message,'	UNITS keyword unreliable!',/info
    endelse
  endelse
endif else begin
  message,'	cannot understand input file: '+arffil,/info
  return,0.
endelse

;	create output structure
ostr=create_struct('NNRG',n_elements(elo),'ELO',elo,'EHI',ehi)

return,arf
end
